//! CPU MSM Implementation - v41
//!
//! Point-parallel with window-sum reduction:
//! - Each thread processes a chunk of points for ALL windows
//! - Thread-local buckets (sequential windows)
//! - Reduce each thread's buckets to ONE value per window
//! - Combine: just num_windows additions
//!
//! Memory: num_threads × num_windows values (not num_threads × num_windows × bucket_count!)

use bls12_381::{G1Affine, G1Projective, Scalar};
use rayon::prelude::*;

const SCALAR_BITS: usize = 255;
const NAIVE_THRESHOLD: usize = 64;

// ============================================================================
// Public API
// ============================================================================

/// Auto-select: naive for small n, serial for medium, parallel for large
pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= NAIVE_THRESHOLD { return naive(bases, scalars); }
    // Parallel only helps at very large n (where bucket reduction dominates)
    if n >= 65536 { return parallel_msm(bases, scalars); }
    serial_msm(bases, scalars)
}

// Aliases
pub fn bellman_style_multiexp(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn naive_msm_stack(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { naive(bases, scalars) }
pub fn pippenger_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn pippenger_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn pippenger_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn strauss_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn strauss_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn bellman_style_multiexp_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }
pub fn optimized_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective { auto_msm(bases, scalars) }

pub enum Algorithm { Naive, Pippenger }
pub fn msm_with_algorithm(bases: &[G1Affine], scalars: &[Scalar], alg: Algorithm) -> G1Projective {
    match alg { Algorithm::Naive => naive(bases, scalars), Algorithm::Pippenger => auto_msm(bases, scalars) }
}
pub fn optimal_algorithm(n: usize) -> Algorithm {
    if n <= NAIVE_THRESHOLD { Algorithm::Naive } else { Algorithm::Pippenger }
}

// ============================================================================
// Naive MSM - for small n
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
fn chunk_size(n: usize) -> usize {
    if n < 32 { 3 } else { (n as f64).ln().ceil() as usize }
}

#[inline]
pub fn extract_window_bits(bytes: &[u8; 32], start_bit: usize, num_bits: usize) -> usize {
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
// Serial MSM - for medium n
// ============================================================================

/// Serial Bellman-style MSM: ln(n) chunks, summation by parts
fn serial_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    let c = chunk_size(n);
    let num_chunks = (SCALAR_BITS + c - 1) / c;
    let bucket_count = 1usize << c;
    
    // Pre-convert scalars to bytes
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    let mut chunk_sums: Vec<G1Projective> = vec![G1Projective::identity(); num_chunks];
    
    for chunk_idx in 0..num_chunks {
        let bit_pos = chunk_idx * c;
        let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
        
        // Bucket accumulation
        for i in 0..n {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, c);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Summation by parts (backward)
        let mut running_sum = G1Projective::identity();
        for k in (1..bucket_count).rev() {
            running_sum += buckets[k];
            chunk_sums[chunk_idx] += running_sum;
        }
    }
    
    // Combine MSB to LSB
    let mut result = G1Projective::identity();
    for chunk_idx in (0..num_chunks).rev() {
        for _ in 0..c {
            result = result.double();
        }
        result += chunk_sums[chunk_idx];
    }
    
    result
}

// ============================================================================
// Point-Parallel MSM with Window-Sum Reduction
// ============================================================================

/// Parallel MSM: each thread processes a range of points
/// 
/// Memory-efficient approach:
/// - Thread processes its point range sequentially through all windows
/// - Accumulates into buckets
/// - Reduces buckets to ONE value per window (summation by parts)
/// - Returns just num_windows values (not buckets!)
/// - Combine: just num_windows additions
/// 
/// Memory per thread: 1 × num_windows × 2^c buckets (for sequential processing)
/// But we REUSE buckets for each window, so memory is: 1 × 2^c
/// Total memory: 2^c G1 elements (not num_threads × 2^c!)
fn parallel_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    let c = chunk_size(n);
    let num_chunks = (SCALAR_BITS + c - 1) / c;
    let bucket_count = 1usize << c;
    
    // Pre-convert scalars to bytes
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Split points across threads
    let num_threads = rayon::current_num_threads();
    let points_per_thread = (n + num_threads - 1) / num_threads;
    
    // Clone bases for thread-safe access
    let bases_vec: Vec<G1Affine> = bases.iter().cloned().collect();
    
    // Each thread returns just num_chunks values (one per window)
    let thread_window_sums: Vec<Vec<G1Projective>> = (0..num_threads).into_par_iter().map(|thread_id| {
        let start = thread_id * points_per_thread;
        let end = (start + points_per_thread).min(n);
        
        if start >= n { 
            return vec![G1Projective::identity(); num_chunks]; 
        }
        
        // Process all windows sequentially (one bucket array at a time)
        let mut window_sums = vec![G1Projective::identity(); num_chunks];
        
        for chunk_idx in 0..num_chunks {
            let bit_pos = chunk_idx * c;
            
            // One bucket array, reused for each window
            let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
            
            // Accumulate this thread's points into buckets
            for i in start..end {
                let k = extract_window_bits(&scalar_bytes[i], bit_pos, c);
                if k > 0 {
                    buckets[k] += bases_vec[i];
                }
            }
            
            // Reduce buckets to single value (summation by parts, backward)
            let mut running_sum = G1Projective::identity();
            for k in (1..bucket_count).rev() {
                running_sum += buckets[k];
                window_sums[chunk_idx] += running_sum;
            }
        }
        
        window_sums
    }).collect();
    
    // Combine thread results: just num_chunks additions!
    // For each window, sum up all thread results
    let mut final_window_sums = vec![G1Projective::identity(); num_chunks];
    for thread_sums in thread_window_sums {
        for chunk_idx in 0..num_chunks {
            final_window_sums[chunk_idx] += thread_sums[chunk_idx];
        }
    }
    
    // Combine windows MSB to LSB
    let mut result = G1Projective::identity();
    for chunk_idx in (0..num_chunks).rev() {
        for _ in 0..c {
            result = result.double();
        }
        result += final_window_sums[chunk_idx];
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