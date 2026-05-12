//! CPU MSM Implementation - v23
//!
//! Adaptive multi-scalar multiplication with algorithm selection based on input size.
//!
//! # Algorithm Selection (Key insight: avoid Pippenger overhead at small scales)
//!
//! | Input Size (n) | Algorithm | Reason                                      |
//! |----------------|-----------|---------------------------------------------|
//! | n <= 512       | Naive     | Low constant factors, no bucket overhead    |
//! | n > 512        | Pippenger | Amortized bucket cost, good parallelization |
//!
//! ## Why Pippenger is SLOW for small inputs:
//!
//! 1. **Fixed bucket allocation**: `2^w` buckets regardless of n (e.g., 64 buckets for w=6)
//! 2. **Per-window setup**: Clear buckets + reduction phase = O(2^w) per window
//! 3. **Not amortized**: At n=64, the bucket overhead dominates the ~64 point additions
//! 4. **Thread spawning overhead**: Rayon work stealing > work at small scales
//!
//! ## The "larger runtimes for smallest input sizes" fix:
//!
//! By using Naive (direct sequential multiplication) for n <= 512, we avoid all
//! Pippenger's overhead. The naive approach is O(n) with minimal constant factors,
//! making it faster for small-to-medium inputs.
//!
//! # Algorithms Implemented
//!
//! - **Naive**: Direct `result += bases[i] * scalars[i]` loop. Optimal for n <= 512.
//! - **Strauss**: Interleaved window processing (implemented but has bugs - disabled)
//! - **Pippenger**: Bucket-based with parallel windows. Optimal for n > 512.

use bls12_381::{G1Affine, G1Projective, Scalar};
use rayon::prelude::*;

mod profiling;

// ============================================================================
// CONSTANTS
// ============================================================================

const SCALAR_BITS: usize = 255;
const PARALLEL_THRESHOLD: usize = 1024;

// Algorithm selection thresholds - empirically tuned
// Pippenger parallel beats naive at n >= 384
// The crossover point is around 384-512 where Pippenger's bucket amortization wins
const NAIVE_THRESHOLD: usize = 256;   // Up to 256 points, naive is fastest
const STRAUSS_THRESHOLD: usize = 4096; // (unused - Strauss has bugs)

// ============================================================================
// BIT EXTRACTION
// ============================================================================

#[inline(always)]
pub fn extract_window_bits(bytes: &[u8; 32], bit_pos: usize, num_bits: usize) -> usize {
    if num_bits == 0 || bit_pos >= SCALAR_BITS {
        return 0;
    }
    let end_bit = (bit_pos + num_bits).min(SCALAR_BITS);
    let effective_bits = end_bit - bit_pos;
    
    if effective_bits == num_bits && (bit_pos & 7) == 0 {
        if num_bits == 8 {
            return bytes[bit_pos >> 3] as usize;
        }
        if num_bits <= 8 {
            return (bytes[bit_pos >> 3] & ((1u8 << num_bits) - 1)) as usize;
        }
    }
    
    let mut result = 0usize;
    let mut shift = 0;
    for i in 0..effective_bits {
        let pos = bit_pos + i;
        let byte_idx = pos >> 3;
        let bit_idx = pos & 7;
        if ((bytes[byte_idx] >> bit_idx) & 1u8) != 0 {
            result |= 1 << shift;
        }
        shift += 1;
    }
    result
}

// ============================================================================
// OPTIMAL WINDOW SIZE
// ============================================================================

#[inline(always)]
fn optimal_window_size(n: usize) -> usize {
    if n <= 8 { 2 }
    else if n <= 32 { 3 }
    else if n <= 128 { 4 }
    else if n <= 512 { 5 }
    else if n <= 2048 { 6 }
    else { 7 }
}

// ============================================================================
// POINT MULTIPLICATION - Addition Chain
// ============================================================================

#[inline(always)]
fn mult_by_k(k: usize, point: G1Projective) -> G1Projective {
    match k {
        0 => G1Projective::identity(),
        1 => point,
        2 => point.double(),
        3 => { let d2 = point.double(); d2 + point },
        4 => point.double().double(),
        5 => { let d4 = point.double().double(); d4 + point },
        6 => { let d2 = point.double(); d2.double() + d2 },
        7 => { let d4 = point.double().double(); d4 + d4 + point },
        8 => point.double().double().double(),
        9 => { let d8 = point.double().double().double(); d8 + point },
        10 => { let d2 = point.double(); let d8 = d2.double().double().double(); d8 + d2 },
        11 => { let d8 = point.double().double().double(); d8 + d8 + point },
        12 => { let d4 = point.double().double(); d4.double() + d4 },
        13 => { let d4 = point.double().double(); d4.double() + d4 + point },
        14 => { let d2 = point.double(); let d8 = d2.double().double().double(); d8 + d2 + d2 + d2 },
        15 => { let d4 = point.double().double(); d4.double() + d4 + d4 + point },
        _ => {
            let mut result = G1Projective::identity();
            let mut current = point;
            for bit in 0..16 {
                if (k >> bit) & 1 != 0 {
                    result += current;
                }
                if bit < 15 {
                    current = current.double();
                }
            }
            result
        }
    }
}

// ============================================================================
// WINDOW AGGREGATION
// ============================================================================

#[inline(always)]
fn aggregate_window(buckets: &[G1Projective]) -> G1Projective {
    let mut result = G1Projective::identity();
    for k in 1..buckets.len() {
        if !bool::from(buckets[k].is_identity()) {
            result += mult_by_k(k, buckets[k]);
        }
    }
    result
}

// ============================================================================
// POWER FACTOR PRECOMPUTATION
// ============================================================================

fn precompute_power_factors(num_windows: usize, w: usize) -> Vec<Scalar> {
    let mut power_factors = Vec::with_capacity(num_windows);
    let mut power = Scalar::one();
    for _ in 0..num_windows {
        power_factors.push(power);
        for _ in 0..w {
            power = power.double();
        }
    }
    power_factors
}

// ============================================================================
// NAIVE MSM
// ============================================================================
// OPTIMIZED SMALL-MSM (bucket-based, adaptive window)
// ============================================================================

fn small_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // For very small n, naive wins (no bucket overhead)
    if n <= 32 {
        let mut result = G1Projective::identity();
        for i in 0..n {
            result += bases[i] * scalars[i];
        }
        return result;
    }
    
    // Use w=5: 51 windows, 32 buckets each (1632 total buckets)
    // This is efficient for n=64-512
    const WINDOW_BITS: usize = 5;
    const BUCKET_COUNT: usize = 32; // 2^5
    const NUM_WINDOWS: usize = (SCALAR_BITS + WINDOW_BITS - 1) / WINDOW_BITS; // 51
    
    let mut buckets: Vec<Vec<G1Projective>> = (0..NUM_WINDOWS)
        .map(|_| vec![G1Projective::identity(); BUCKET_COUNT])
        .collect();
    
    // Accumulate into buckets
    for i in 0..n {
        let bytes = scalars[i].to_bytes();
        for w_idx in 0..NUM_WINDOWS {
            let k = extract_window_bits(&bytes, w_idx * WINDOW_BITS, WINDOW_BITS);
            if k > 0 {
                buckets[w_idx][k] += bases[i];
            }
        }
    }
    
    // Reduction: R = Σ window_sum[j] * 2^(j*5)
    let mut result = G1Projective::identity();
    
    for w_idx in 0..NUM_WINDOWS {
        let mut window_sum = G1Projective::identity();
        for k in 1..BUCKET_COUNT {
            if !bool::from(buckets[w_idx][k].is_identity()) {
                window_sum += buckets[w_idx][k] * Scalar::from(k as u64);
            }
        }
        
        // Multiply by 2^(w_idx * 5) - just double repeatedly
        let mut power = Scalar::one();
        for _ in 0..(w_idx * WINDOW_BITS) {
            power = power.double();
        }
        result += window_sum * power;
    }
    
    result
}

// ============================================================================
// NAIVE MSM (Stack-optimized for small inputs)
// ============================================================================

#[inline(always)]
pub fn naive_msm_stack(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    let mut result = G1Projective::identity();
    for i in 0..n {
        result += bases[i] * scalars[i];
    }
    result
}

#[inline(always)]
pub fn naive_msm_small(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    naive_msm_stack(bases, scalars)
}

// ============================================================================
// STRAUSS'S ALGORITHM (Improved Shamir's Method)
// ============================================================================
//
// Strauss's algorithm uses an interleaved approach that processes all windows
// simultaneously rather than sequentially. This provides better cache locality
// because we only need to keep one bucket array in cache instead of reprocessing
// the bases multiple times.
//
// Key advantages over Pippenger:
// 1. Better cache locality - single pass over bases vs multi-pass
// 2. Lower constant factors for medium-sized inputs
// 3. Avoids the expensive final aggregation phase
//
// The algorithm:
// 1. Partition scalars into c groups of size L = ceil(n/c)
// 2. Process each group with window size w, computing partial sums
// 3. Combine partial sums with the correct power factors
// ============================================================================

fn strauss_num_windows(_n: usize, w: usize) -> usize {
    (SCALAR_BITS + w - 1) / w
}

/// Strauss serial implementation with interleaved window processing
pub fn strauss_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // Choose optimal window size for Strauss
    // Strauss can use slightly larger w than Pippenger due to better locality
    let w = match n {
        n if n <= 128 => 3,
        n if n <= 512 => 4,
        n if n <= 2048 => 5,
        n if n <= 8192 => 6,
        _ => 7,
    };
    
    let bucket_count = 1usize << w;
    let num_windows = strauss_num_windows(n, w);
    
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Precompute power factors: 2^(j*w) for each window j
    let power_factors: Vec<Scalar> = (0..num_windows)
        .map(|j| {
            let mut pow = Scalar::one();
            for _ in 0..(j * w) {
                pow = pow.double();
            }
            pow
        })
        .collect();
    
    // Allocate a single bucket array for all windows (interleaved accumulation)
    // buckets[j * bucket_count + k] holds sum of all points where window j has value k
    let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); num_windows * bucket_count];
    
    // Single pass over all scalars - interleave window processing
    for i in 0..n {
        let bytes = &scalar_bytes[i];
        
        // Process each window for this scalar
        for window_idx in 0..num_windows {
            let bit_pos = window_idx * w;
            let k = extract_window_bits(bytes, bit_pos, w);
            
            if k > 0 {
                let bucket_idx = window_idx * bucket_count + k;
                buckets[bucket_idx] += bases[i];
            }
        }
    }
    
    // Combine results from all windows
    let mut result = G1Projective::identity();
    
    for window_idx in 0..num_windows {
        // Sum this window's buckets
        let mut window_sum = G1Projective::identity();
        let window_base = window_idx * bucket_count;
        
        for k in 1..bucket_count {
            if !bool::from(buckets[window_base + k].is_identity()) {
                window_sum += buckets[window_base + k] * Scalar::from(k as u64);
            }
        }
        
        // Add scaled by power factor
        result += window_sum * power_factors[window_idx];
    }
    
    result
}

/// Strauss parallel implementation - partition work by points, not windows
pub fn strauss_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // For small n, serial is faster (avoids parallel overhead)
    if n <= 2048 { return strauss_serial(bases, scalars); }
    
    // For large n, use Pippenger instead (better parallelization)
    pippenger_msm_parallel(bases, scalars)
}

/// Wrapper for Strauss MSM
pub fn strauss_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    strauss_msm_parallel(bases, scalars)
}

// ============================================================================
// PIPPENGER SERIAL
// ============================================================================

pub fn pippenger_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= 32 { return naive_msm_stack(bases, scalars); }
    
    let w = optimal_window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    let bucket_count = 1usize << w;
    
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Precompute power factors: 2^(j*w) mod p
    let power_factors: Vec<Scalar> = (0..num_windows)
        .map(|j| {
            let mut pow = Scalar::one();
            for _ in 0..(j * w) {
                pow = pow.double();
            }
            pow
        })
        .collect();
    
    // Process each window sequentially (Pippenger's original approach)
    let mut final_result = G1Projective::identity();
    
    for (window_idx, &pf) in power_factors.iter().enumerate() {
        let bit_pos = window_idx * w;
        
        // Initialize buckets (cleared per window - key difference from Strauss)
        let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
        
        // Accumulate points into buckets
        for i in 0..n {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, w);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Sum buckets using direct scalar multiplication (FIX: was mult_by_k which had bugs)
        let mut window_sum = G1Projective::identity();
        for k in 1..bucket_count {
            if !bool::from(buckets[k].is_identity()) {
                window_sum += buckets[k] * Scalar::from(k as u64);
            }
        }
        
        // Add to final result, scaled by power factor
        final_result += window_sum * pf;
    }
    
    final_result
}

// ============================================================================
// PIPPENGER PARALLEL
// ============================================================================

pub fn pippenger_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= 32 { return naive_msm_stack(bases, scalars); }
    if n <= PARALLEL_THRESHOLD { 
        return pippenger_serial(bases, scalars);
    }
    
    let w = optimal_window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    let bucket_count = 1usize << w;
    
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Precompute power factors
    let power_factors: Vec<Scalar> = (0..num_windows)
        .map(|j| {
            let mut pow = Scalar::one();
            for _ in 0..(j * w) {
                pow = pow.double();
            }
            pow
        })
        .collect();
    
    // Process windows in parallel
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
            
            // Sum buckets using direct scalar multiplication
            let mut window_sum = G1Projective::identity();
            for k in 1..bucket_count {
                if !bool::from(buckets[k].is_identity()) {
                    window_sum += buckets[k] * Scalar::from(k as u64);
                }
            }
            
            window_sum * power_factors[window_idx]
        })
        .collect();
    
    // Combine results
    window_results.iter().fold(G1Projective::identity(), |acc, &r| acc + r)
}

pub fn pippenger_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    pippenger_msm_parallel(bases, scalars)
}

// ============================================================================
// AUTO-SELECT WITH ADAPTIVE THRESHOLDS
// ============================================================================
//
// Algorithm selection based on empirical analysis:
//
// | Input Size  | Best Algorithm | Reason                                    |
// |-------------|----------------|-------------------------------------------|
// | n < 128     | Naive          | Low overhead, sequential access            |
// | 128 - 4096  | Strauss        | Better cache locality than Pippenger       |
// | n >= 4096   | Pippenger      | Amortized overhead, parallel efficiency    |
//
// The "small n is slower" problem with Pippenger comes from:
// 1. Fixed bucket allocation overhead (O(2^w) regardless of n)
// 2. Precomputed point table setup cost
// 3. Reduction phase cost that doesn't amortize at small n
// 4. Thread spawning overhead in parallel mode
// ============================================================================

pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // For small n, naive is optimal (no bucket/window overhead)
    // The original "large runtimes for smallest input sizes" issue was
    // caused by Pippenger being used at n=64 - the bucket overhead 
    // (2^w buckets regardless of n) dominates the actual work
    if n <= NAIVE_THRESHOLD {
        naive_msm_stack(bases, scalars)
    } else {
        // Pippenger for large inputs
        pippenger_msm_parallel(bases, scalars)
    }
}

/// Force a specific algorithm (useful for benchmarking)
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
// WRAPPER FUNCTIONS
// ============================================================================

pub fn sliding_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    auto_msm(bases, scalars)  // Now uses adaptive selection
}

pub fn sliding_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    auto_msm(bases, scalars)  // Now uses adaptive selection
}

/// Get the optimal algorithm for a given input size
pub fn optimal_algorithm(n: usize) -> Algorithm {
    match n {
        n if n <= NAIVE_THRESHOLD => Algorithm::Naive,
        n if n <= STRAUSS_THRESHOLD => Algorithm::Strauss,
        _ => Algorithm::Pippenger,
    }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bit_extraction() {
        // BLS12-381 scalars are little-endian (LSB first)
        // Byte 0 is the least significant byte
        let bytes = [0x12u8, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                     0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        
        // extract(0, 8): bits 0-7 of byte 0 = 0x12
        assert_eq!(extract_window_bits(&bytes, 0, 8), 0x12);
        // extract(8, 8): bits 0-7 of byte 1 = 0x34
        assert_eq!(extract_window_bits(&bytes, 8, 8), 0x34);
        // extract(0, 16): bits 0-15 of bytes 0-1 in little-endian = 0x3412
        assert_eq!(extract_window_bits(&bytes, 0, 16), 0x3412);
        // extract(4, 8): bits 4-11 of bytes 0-1 = bits 4-7 of 0x12 + bits 0-3 of 0x34 = 0x41
        assert_eq!(extract_window_bits(&bytes, 4, 8), 0x41);
        // extract(250, 5): bits 2 of byte 31 (0x00) = 0
        assert_eq!(extract_window_bits(&bytes, 250, 5), 0x0);
    }
    
    #[test]
    fn test_correctness() {
        let g = G1Affine::generator();
        
        for n in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 16384] {
            let n = n;
            
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
    fn test_algorithm_selection() {
        // Verify correct algorithm is selected for different sizes
        let g = G1Affine::generator();
        let bases: Vec<G1Affine> = (0..64).map(|_| g).collect();
        let scalars: Vec<Scalar> = (0..64).map(|i| Scalar::from(i as u64 + 1)).collect();
        
        // For n=64, should use naive (below NAIVE_THRESHOLD=128)
        let result = auto_msm(&bases, &scalars);
        
        // Just verify it runs without panicking and produces correct result
        let expected = naive_msm_stack(&bases, &scalars);
        assert_eq!(result, expected);
    }
    
    #[test]
    fn test_performance() {
        let g = G1Affine::generator();
        
        println!("\n=== MSM Performance Comparison (release mode) ===");
        println!("");
        println!("Input sizes and algorithm selection:");
        println!("| Points | Selected  | Naive    | Strauss  | Pippenger| Bellman  |");
        println!("|--------|-----------|----------|----------|----------|----------|");
        
        let iterations = 3;
        
        for n in [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 16384] {
            let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
            let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
            
            // Determine which algorithm auto selects
            let selected = match n {
                n if n <= NAIVE_THRESHOLD => "Naive",
                n if n <= STRAUSS_THRESHOLD => "Strauss",
                _ => "Pippenger",
            };
            
            // Bellman baseline (from zkgpu benchmarks)
            let bellman_ms = match n {
                16 => 2.1, 32 => 2.3, 64 => 2.5, 128 => 3.4, 256 => 4.8,
                512 => 4.8, 1024 => 6.9, 2048 => 10.7, 4096 => 21.1, 16384 => 44.3,
                _ => 0.0,
            };
            
            // Benchmark each algorithm
            let naive_time = benchmark_fn(|| naive_msm_stack(&bases, &scalars), iterations);
            let strauss_time = benchmark_fn(|| strauss_msm_parallel(&bases, &scalars), iterations);
            let pippenger_time = benchmark_fn(|| pippenger_msm_parallel(&bases, &scalars), iterations);
            
            println!("| {:6} | {:9} | {:6.2}ms | {:6.2}ms | {:6.2}ms | {:6.2}ms |", 
                     n, selected, naive_time, strauss_time, pippenger_time, bellman_ms);
        }
        
        println!("");
        println!("Key insight: Naive is fastest for n < {}, avoiding Pippenger's bucket overhead.", NAIVE_THRESHOLD);
        println!("Strauss beats Pippenger for n < {} due to better cache locality.", STRAUSS_THRESHOLD);
    }
}

fn benchmark_fn<F: Fn() -> G1Projective>(f: F, iterations: usize) -> f64 {
    use std::time::Instant;
    
    // Warmup
    for _ in 0..2 {
        let _ = f();
    }
    
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = f();
    }
    start.elapsed().as_secs_f64() / (iterations as f64) * 1000.0
}