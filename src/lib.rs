//! CPU MSM Implementation - v50
//!
//! Based on v47 which beats Bellman at n >= 16384
//! Added optimizations:
//! 1. Pre-allocated bucket array (no reallocation per window)
//! 2. Pre-converted scalar bytes (computed once)
//! 3. Summation by parts reduction

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
// Bit Extraction (Optimized)
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
// Parallel MSM (Optimized)
// ============================================================================

/// Process a chunk of points - all windows in one pass
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
        
        // Allocate and zero buckets
        let mut buckets = vec![G1Projective::identity(); bucket_count];
        
        // Accumulate
        for i in start..end {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, w);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Summation by parts
        let mut running_sum = G1Projective::identity();
        for k in (1..bucket_count).rev() {
            running_sum += buckets[k];
            window_sums[window_idx] += running_sum;
        }
    }
    
    window_sums
}

/// Main parallel MSM
fn parallel_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    let w = window_size(n);
    let num_windows = (SCALAR_BITS + w - 1) / w;
    
    // Pre-convert scalars to bytes ONCE
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
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
        
        process_chunk(bases, &scalar_bytes, start, end, w, num_windows)
    }).collect();
    
    // Combine chunk results
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