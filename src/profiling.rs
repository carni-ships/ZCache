//! Profiling utilities for CPU MSM implementations

use crate::{optimal_window_size, extract_window_bits};

/// Phase timing breakdown
#[derive(Default)]
pub struct PhaseBreakdown {
    pub total_ms: f64,
    pub scalar_cache_ms: f64,
    pub bucket_accumulation_ms: f64,
    pub power_factor_ms: f64,
    pub window_aggregation_ms: f64,
    pub bucket_additions: u64,
    pub bucket_doublings: u64,
    pub direct_muls: u64,
}

/// Profile a Pippenger implementation with phase timings
pub fn profile_pippenger<F>(bases: &[bls12_381::G1Affine], scalars: &[bls12_381::Scalar], func: F) -> (bls12_381::G1Projective, PhaseBreakdown)
where
    F: Fn(&[bls12_381::G1Affine], &[bls12_381::Scalar]) -> bls12_381::G1Projective,
{
    use std::time::Instant;
    use bls12_381::G1Projective;
    
    let start = Instant::now();
    
    let n = bases.len();
    if n == 0 {
        return (G1Projective::identity(), PhaseBreakdown::default());
    }
    
    let w = optimal_window_size(n);
    let num_windows = (255 + w - 1) / w;
    let bucket_count = 1usize << w;
    
    // Phase 1: Scalar cache
    let scalar_start = Instant::now();
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    let scalar_cache_ms = scalar_start.elapsed().as_secs_f64() * 1000.0;
    
    // Phase 2: Bucket accumulation
    let bucket_start = Instant::now();
    let mut window_buckets: Vec<Vec<G1Projective>> = Vec::with_capacity(num_windows);
    for _ in 0..num_windows {
        window_buckets.push(vec![G1Projective::identity(); bucket_count]);
    }
    
    let mut bucket_additions = 0u64;
    for (i, base) in bases.iter().enumerate() {
        let bytes = &scalar_bytes[i];
        
        for window_idx in 0..num_windows {
            let bit_pos = window_idx * w;
            let k = extract_window_bits(bytes, bit_pos, w);
            
            if k > 0 && k < bucket_count {
                window_buckets[window_idx][k] += *base;
                bucket_additions += 1;
            }
        }
    }
    let bucket_accumulation_ms = bucket_start.elapsed().as_secs_f64() * 1000.0;
    
    // Phase 3: Power factors
    let power_start = Instant::now();
    let mut scalar_2_pow_w = bls12_381::Scalar::one();
    for _ in 0..w {
        scalar_2_pow_w = scalar_2_pow_w.double();
    }
    
    let mut power_factors = vec![bls12_381::Scalar::one(); num_windows];
    for i in 1..num_windows {
        power_factors[i] = power_factors[i - 1] * scalar_2_pow_w;
    }
    let power_factor_ms = power_start.elapsed().as_secs_f64() * 1000.0;
    
    // Phase 4: Window aggregation
    let agg_start = Instant::now();
    let mut bucket_doublings = 0u64;
    let mut direct_muls = 0u64;
    
    let mut result = G1Projective::identity();
    
    for window_idx in 0..num_windows {
        let buckets = &window_buckets[window_idx];
        
        let mut window_result = G1Projective::identity();
        
        for k in 1..bucket_count {
            let bucket = buckets[k];
            if bool::from(bucket.is_identity()) {
                continue;
            }
            
            if k <= 16 {
                window_result += bucket * bls12_381::Scalar::from_raw([k as u64, 0, 0, 0]);
                direct_muls += 1;
            } else {
                let mut current = bucket;
                let mut remaining = k;
                while remaining > 0 {
                    if (remaining & 1) != 0 {
                        window_result += current;
                    }
                    remaining >>= 1;
                    if remaining > 0 {
                        current = current.double();
                        bucket_doublings += 1;
                    }
                }
            }
        }

        if !bool::from(window_result.is_identity()) {
            result += window_result * power_factors[window_idx];
        }
    }
    let window_aggregation_ms = agg_start.elapsed().as_secs_f64() * 1000.0;
    
    let total_ms = start.elapsed().as_secs_f64() * 1000.0;
    
    (
        result,
        PhaseBreakdown {
            total_ms,
            scalar_cache_ms,
            bucket_accumulation_ms,
            power_factor_ms,
            window_aggregation_ms,
            bucket_additions,
            bucket_doublings,
            direct_muls,
        },
    )
}

/// Print phase breakdown
pub fn print_phase_breakdown(n: usize, phase: &PhaseBreakdown) {
    println!("\n  Phase Breakdown for n={}:", n);
    println!("  ├── Scalar cache:        {:>8.3} ms ({:>5.1}%)", 
             phase.scalar_cache_ms, 100.0 * phase.scalar_cache_ms / phase.total_ms);
    println!("  ├── Bucket accumulation: {:>8.3} ms ({:>5.1}%)", 
             phase.bucket_accumulation_ms, 100.0 * phase.bucket_accumulation_ms / phase.total_ms);
    println!("  ├── Power factors:       {:>8.3} ms ({:>5.1}%)", 
             phase.power_factor_ms, 100.0 * phase.power_factor_ms / phase.total_ms);
    println!("  └── Window aggregation:  {:>8.3} ms ({:>5.1}%)", 
             phase.window_aggregation_ms, 100.0 * phase.window_aggregation_ms / phase.total_ms);
    println!("  └── Total:              {:>8.3} ms", phase.total_ms);
    println!("  Statistics:");
    println!("  ├── Bucket additions:   {:>10}", phase.bucket_additions);
    println!("  ├── Point doublings:    {:>10}", phase.bucket_doublings);
    println!("  ├── Direct muls (k≤16):{:>10}", phase.direct_muls);
    
    // Determine bottleneck
    let bottleneck = analyze_bottleneck(n, phase);
    println!("  └── BOTTLENECK: {}", bottleneck);
}

/// Analyze which phase is the bottleneck
pub fn analyze_bottleneck(n: usize, phase: &PhaseBreakdown) -> &'static str {
    let acc_pct = 100.0 * phase.bucket_accumulation_ms / phase.total_ms;
    let agg_pct = 100.0 * phase.window_aggregation_ms / phase.total_ms;
    
    if acc_pct > 50.0 {
        "bucket_accumulation"
    } else if agg_pct > 50.0 {
        "window_aggregation"
    } else if acc_pct > agg_pct {
        "bucket_accumulation"
    } else {
        "window_aggregation"
    }
}