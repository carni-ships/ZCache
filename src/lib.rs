//! CPU MSM Implementation - v25 (Clean Rewrite)
//!
//! Adaptive multi-scalar multiplication optimized for BLS12-381 G1.
//!
//! # Algorithm Selection (Empirically Tuned)
//!
//! | Input Size (n) | Algorithm | Reason                              |
//! |----------------|-----------|-------------------------------------|
//! | n <= 256       | Naive     | No bucket overhead, fastest at small |
//! | n > 256        | Pippenger | Parallel windows, amortized cost    |
//!
//! # Why We Beat/Lose to Bellman
//!
//! Bellman uses Pippenger with w=5-6 and point-first parallelization.
//! We use the same algorithm but with optimized:
//! - Bit extraction (single-pass per window)
//! - Bucket accumulation (identity skip)
//! - Window reduction (direct scalar multiplication)
//!
//! We lose at n=256-512 because our Pippenger has parallel overhead that
//! Bellman avoids (they use serial Pippenger for medium sizes).

use bls12_381::{G1Affine, G1Projective, Scalar};
use rayon::prelude::*;

// ============================================================================
// Constants
// ============================================================================

const SCALAR_BITS: usize = 255;
// Algorithm thresholds
// Parallel Pippenger beats naive at n >= 64 (due to parallelization)
const NAIVE_THRESHOLD: usize = 64;  // Naive for n <= 64

// ============================================================================
// Bit Extraction (Optimized for BLS12-381 scalars)
// ============================================================================

/// Extract `num_bits` bits starting at `start_bit` from a 32-byte scalar.
/// Uses little-endian byte order (BLS12-381 convention).
#[inline(always)]
fn extract_window_bits(bytes: &[u8; 32], start_bit: usize, num_bits: usize) -> usize {
    let mut result = 0usize;
    let mut bit_pos = start_bit;
    let mut byte_idx = start_bit / 8;
    let mut bits_extracted = 0;
    
    while bits_extracted < num_bits && byte_idx < 32 {
        let bits_in_byte = 8 - (bit_pos % 8);
        let bits_to_take = bits_in_byte.min(num_bits - bits_extracted);
        
        let mask = (1usize << bits_to_take) - 1;
        result |= ((bytes[byte_idx] >> (bit_pos % 8)) as usize & mask) << bits_extracted;
        
        bit_pos += bits_to_take;
        bits_extracted += bits_to_take;
        byte_idx += 1;
    }
    
    result
}

// ============================================================================
// Naive MSM (Optimal for n <= 256)
// ============================================================================

/// Direct scalar multiplication: result = Σ bases[i] * scalars[i]
/// No bucket overhead - fastest for small inputs.
#[inline]
pub fn naive_msm_stack(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    let mut result = G1Projective::identity();
    for i in 0..n {
        result += bases[i] * scalars[i];
    }
    result
}

#[inline]
pub fn naive_msm_small(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    naive_msm_stack(bases, scalars)
}

// ============================================================================
// Pippenger MSM (Optimal for n > 256)
// ============================================================================

/// Optimal window size based on input size
fn optimal_window_size(n: usize) -> usize {
    match n {
        n if n <= 256 => 4,
        n if n <= 1024 => 5,
        n if n <= 4096 => 6,
        n if n <= 16384 => 7,
        _ => 8,
    }
}

/// Pippenger serial - single-threaded window processing
pub fn pippenger_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= 32 { return naive_msm_stack(bases, scalars); }
    
    let w = optimal_window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    let bucket_count = 1usize << w;
    
    // Convert scalars to bytes once
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Precompute 2^(j*w) for each window
    let power_factors: Vec<Scalar> = (0..num_windows)
        .map(|j| {
            let mut pow = Scalar::one();
            for _ in 0..(j * w) {
                pow = pow.double();
            }
            pow
        })
        .collect();
    
    // Process windows sequentially
    let mut final_result = G1Projective::identity();
    
    for window_idx in 0..num_windows {
        let bit_pos = window_idx * w;
        let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
        
        // Accumulate points into buckets
        for i in 0..n {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, w);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Reduce: sum k * bucket[k] for k = 1..bucket_count
        let mut window_sum = G1Projective::identity();
        for k in 1..bucket_count {
            if !bool::from(buckets[k].is_identity()) {
                window_sum += buckets[k] * Scalar::from(k as u64);
            }
        }
        
        final_result += window_sum * power_factors[window_idx];
    }
    
    final_result
}

/// Pippenger parallel - parallelize over windows for large n only
pub fn pippenger_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= 64 { return naive_msm_stack(bases, scalars); }
    
    let w = optimal_window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    let bucket_count = 1usize << w;
    
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    let power_factors: Vec<Scalar> = (0..num_windows)
        .map(|j| {
            let mut pow = Scalar::one();
            for _ in 0..(j * w) {
                pow = pow.double();
            }
            pow
        })
        .collect();
    
    let window_results: Vec<G1Projective> = (0..num_windows)
        .into_par_iter()
        .map(|window_idx| {
            let bit_pos = window_idx * w;
            let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
            
            for i in 0..n {
                let k = extract_window_bits(&scalar_bytes[i], bit_pos, w);
                if k > 0 {
                    buckets[k] += bases[i];
                }
            }
            
            let mut window_sum = G1Projective::identity();
            for k in 1..bucket_count {
                if !bool::from(buckets[k].is_identity()) {
                    window_sum += buckets[k] * Scalar::from(k as u64);
                }
            }
            
            window_sum * power_factors[window_idx]
        })
        .collect();
    
    window_results.iter().fold(G1Projective::identity(), |acc, &r| acc + r)
}

pub fn pippenger_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    pippenger_msm_parallel(bases, scalars)
}

// ============================================================================
// Strauss MSM (Interleaved windows - cache-friendly for medium n)
// ============================================================================

/// Number of windows for Strauss with given window size
fn strauss_num_windows(_n: usize, w: usize) -> usize {
    (SCALAR_BITS + w - 1) / w
}

/// Strauss serial - all windows processed in single pass over bases
fn strauss_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // Use w=5 (32 buckets, 51 windows)
    let w = 5;
    let bucket_count = 1usize << w;
    let num_windows = strauss_num_windows(n, w);
    
    // One bucket array for all windows (Strauss key difference)
    let mut buckets: Vec<Vec<G1Projective>> = (0..num_windows)
        .map(|_| vec![G1Projective::identity(); bucket_count])
        .collect();
    
    // Single pass over bases - update all windows simultaneously
    for i in 0..n {
        let bytes = scalars[i].to_bytes();
        
        for window_idx in 0..num_windows {
            let k = extract_window_bits(&bytes, window_idx * w, w);
            if k > 0 {
                buckets[window_idx][k] += bases[i];
            }
        }
    }
    
    // Reduction: combine all windows
    let mut result = G1Projective::identity();
    
    for window_idx in 0..num_windows {
        let mut window_sum = G1Projective::identity();
        for k in 1..bucket_count {
            if !bool::from(buckets[window_idx][k].is_identity()) {
                window_sum += buckets[window_idx][k] * Scalar::from(k as u64);
            }
        }
        
        // Multiply by 2^(window_idx * w)
        let mut power = Scalar::one();
        for _ in 0..(window_idx * w) {
            power = power.double();
        }
        result += window_sum * power;
    }
    
    result
}

/// Strauss parallel - partition by points, not windows
pub fn strauss_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // For small n, serial is faster
    if n <= 1024 { return strauss_serial(bases, scalars); }
    
    // For large n, use Pippenger (better parallelization)
    pippenger_msm_parallel(bases, scalars)
}

pub fn strauss_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    strauss_msm_parallel(bases, scalars)
}

// ============================================================================
// Auto-select MSM (Main Entry Point)
// ============================================================================

/// Automatically select optimal algorithm based on input size.
pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    if n <= NAIVE_THRESHOLD {
        naive_msm_stack(bases, scalars)
    } else {
        pippenger_msm_parallel(bases, scalars)  // Use parallel for all larger inputs
    }
}

/// Force a specific algorithm (for benchmarking)
pub fn msm_with_algorithm(
    bases: &[G1Affine],
    scalars: &[Scalar],
    algorithm: Algorithm
) -> G1Projective {
    match algorithm {
        Algorithm::Naive => naive_msm_stack(bases, scalars),
        Algorithm::Strauss => strauss_msm_parallel(bases, scalars),
        Algorithm::Pippenger => pippenger_msm_parallel(bases, scalars),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    Naive,
    Strauss,
    Pippenger,
}

// ============================================================================
// Wrapper Functions (API Compatibility)
// ============================================================================

pub fn sliding_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    auto_msm(bases, scalars)
}

pub fn sliding_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    auto_msm(bases, scalars)
}

/// Get the optimal algorithm for a given input size
pub fn optimal_algorithm(n: usize) -> Algorithm {
    match n {
        n if n <= NAIVE_THRESHOLD => Algorithm::Naive,
        n if n <= NAIVE_THRESHOLD * 4 => Algorithm::Strauss,
        _ => Algorithm::Pippenger,
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bit_extraction() {
        let bytes = [0x12u8, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        
        // Little-endian: byte 0 is LSB
        assert_eq!(extract_window_bits(&bytes, 0, 8), 0x12);
        assert_eq!(extract_window_bits(&bytes, 8, 8), 0x34);
        assert_eq!(extract_window_bits(&bytes, 0, 16), 0x3412);
        assert_eq!(extract_window_bits(&bytes, 4, 8), 0x41);
    }
    
    #[test]
    fn test_correctness() {
        let g = G1Affine::generator();
        
        for n in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384] {
            let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
            let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
            
            let naive_result = naive_msm_stack(&bases, &scalars);
            let serial_result = pippenger_serial(&bases, &scalars);
            let parallel_result = pippenger_msm_parallel(&bases, &scalars);
            let strauss_result = strauss_msm_parallel(&bases, &scalars);
            let auto_result = auto_msm(&bases, &scalars);
            
            assert_eq!(naive_result, serial_result, "n={}: serial mismatch", n);
            assert_eq!(naive_result, parallel_result, "n={}: parallel mismatch", n);
            assert_eq!(naive_result, strauss_result, "n={}: strauss mismatch", n);
            assert_eq!(naive_result, auto_result, "n={}: auto mismatch", n);
        }
    }
    
    #[test]
    fn test_performance() {
        use std::time::Instant;
        
        let g = G1Affine::generator();
        
        println!("\n=== MSM Performance (release mode) ===");
        println!("|  n   |  Naive  | Pippenger | Strauss | Bellman |");
        
        for n in [64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384] {
            let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
            let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
            
            let bellman_ms = match n {
                64 => 20.0, 128 => 30.0, 256 => 35.0, 512 => 60.0,
                1024 => 90.0, 2048 => 120.0, 4096 => 180.0,
                8192 => 280.0, 16384 => 400.0, _ => 0.0,
            };
            
            let naive_time = {
                let start = Instant::now();
                for _ in 0..3 {
                    let _ = naive_msm_stack(&bases, &scalars);
                }
                start.elapsed().as_secs_f64() * 1000.0 / 3.0
            };
            
            let pippenger_time = {
                let start = Instant::now();
                for _ in 0..3 {
                    let _ = pippenger_msm_parallel(&bases, &scalars);
                }
                start.elapsed().as_secs_f64() * 1000.0 / 3.0
            };
            
            let strauss_time = {
                let start = Instant::now();
                for _ in 0..3 {
                    let _ = strauss_msm_parallel(&bases, &scalars);
                }
                start.elapsed().as_secs_f64() * 1000.0 / 3.0
            };
            
            println!("| {:5} | {:6.1}ms | {:8.1}ms | {:7.1}ms | {:6.1}ms |",
                     n, naive_time, pippenger_time, strauss_time, bellman_ms);
        }
    }
}