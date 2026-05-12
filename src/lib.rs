//! CPU MSM Implementation - v45
//!
//! Point-parallelization strategy (same as Bellman):
//! 1. Split points across threads (each thread gets a range)
//! 2. Each thread processes ALL windows sequentially with thread-local buckets
//! 3. Each thread returns partial sums for each window
//! 4. Combine: sum thread results per window, then windows MSB to LSB
//!
//! Key: Thread-local storage means no mutex/race conditions.
//! The "window-sums" approach means we only return 1 value per window per thread.

use bls12_381::{G1Affine, G1Projective, Scalar};
use rayon::prelude::*;

const SCALAR_BITS: usize = 255;
const NAIVE_THRESHOLD: usize = 64;

// ============================================================================
// Public API
// ============================================================================

pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= NAIVE_THRESHOLD { return naive(bases, scalars); }
    parallel_msm(bases, scalars)
}

pub fn naive_msm_stack(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { naive(bases, scalars) }
pub fn bellman_style_multiexp(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn pippenger_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn pippenger_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn pippenger_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn strauss_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn strauss_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn optimized_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }

pub enum Algorithm { Naive, Pippenger }
pub fn msm_with_algorithm(bases: &[G1Affine], scalars: &[Scalar], alg: Algorithm) -> G1Projective {
    match alg { Algorithm::Naive => naive(bases, scalars), Algorithm::Pippenger => auto_msm(bases, scalars) }
}
pub fn optimal_algorithm(n: usize) -> Algorithm {
    if n <= NAIVE_THRESHOLD { Algorithm::Naive } else { Algorithm::Pippenger }
}

// ============================================================================
// Naive MSM
// ============================================================================

#[inline]
fn naive(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let mut result = G1Projective::identity();
    for i in 0..bases.len() {
        result += bases[i] * scalars[i];
    }
    result
}

// ============================================================================
// Bit Extraction
// ============================================================================

#[inline]
fn window_size(n: usize) -> usize {
    match n {
        n if n <= 256 => 4,
        n if n <= 1024 => 5,
        n if n <= 4096 => 6,
        _ => 7,
    }
}

#[inline]
fn extract_window_bits(bytes: &[u8; 32], start_bit: usize, num_bits: usize) -> usize {
    let mut result = 0usize;
    let mut bit_pos = start_bit;
    let mut byte_idx = bit_pos / 8;
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
// Parallel MSM - Point-Parallel with Thread-Local Buckets
// ============================================================================

/// Process a chunk of points through all windows
/// Returns partial sums for each window (one G1Projective per window)
fn process_chunk(
    bases: &[G1Affine],
    scalar_bytes: &[[u8; 32]],
    start: usize,
    end: usize,
    w: usize,
    num_windows: usize,
) -> Vec<G1Projective> {
    let bucket_count = 1usize << w;
    let mut window_sums = vec![G1Projective::identity(); num_windows];
    
    // Process each window
    for window_idx in 0..num_windows {
        let bit_pos = window_idx * w;
        let mut buckets = vec![G1Projective::identity(); bucket_count];
        
        // Accumulate this chunk's points into buckets
        for i in start..end {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, w);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Summation by parts: compute sum(k * bucket[k])
        let mut running_sum = G1Projective::identity();
        for k in (1..bucket_count).rev() {
            running_sum += buckets[k];
            window_sums[window_idx] += running_sum;
        }
    }
    
    window_sums
}

/// Point-parallel MSM
/// 
/// Strategy (same as Bellman's multicore approach):
/// 1. Split points into N chunks (N = number of threads)
/// 2. Each thread processes its chunk through ALL windows
/// 3. Each thread returns partial sums for each window
/// 4. Combine all partial sums
/// 5. Combine windows MSB to LSB
fn parallel_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    let w = window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    
    // Pre-convert scalars to bytes
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Clone bases for thread-safe access
    let bases_vec: Vec<G1Affine> = bases.iter().cloned().collect();
    
    // Split points across threads
    let num_threads = rayon::current_num_threads();
    let points_per_thread = (n + num_threads - 1) / num_threads;
    
    // Parallel over point chunks
    let chunk_results: Vec<Vec<G1Projective>> = (0..num_threads).into_par_iter().map(|thread_id| {
        let start = thread_id * points_per_thread;
        let end = (start + points_per_thread).min(n);
        
        if start >= n {
            return vec![G1Projective::identity(); num_windows];
        }
        
        process_chunk(&bases_vec, &scalar_bytes, start, end, w, num_windows)
    }).collect();
    
    // Combine chunk results: sum up values for each window
    let mut final_sums = vec![G1Projective::identity(); num_windows];
    for thread_result in chunk_results {
        for window_idx in 0..num_windows {
            final_sums[window_idx] += thread_result[window_idx];
        }
    }
    
    // Combine windows MSB to LSB
    let mut result = G1Projective::identity();
    for window_idx in (0..num_windows).rev() {
        for _ in 0..w {
            result = result.double();
        }
        result += final_sums[window_idx];
    }
    
    result
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
        assert_eq!(extract_window_bits(&bytes, 0, 8), 0x12);
        assert_eq!(extract_window_bits(&bytes, 8, 8), 0x34);
    }
    
    #[test]
    fn test_correctness() {
        let g = G1Affine::generator();
        for n in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384] {
            let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
            let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
            
            let naive_result = naive(&bases, &scalars);
            let parallel_result = auto_msm(&bases, &scalars);
            
            assert_eq!(naive_result, parallel_result, "n={}", n);
        }
    }
}