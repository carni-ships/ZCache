//! CPU MSM Implementation - v30 (Parallel Window)
//!
//! Optimizations inspired by Bellman's multiexp implementation:
//! 1. ln(n) Chunk Size - Variable window size based on input size
//! 2. Summation by Parts - O(2^c) reduction
//! 3. Density Tracking - Skip empty windows
//! 4. Window-Parallel - Process each chunk in parallel using rayon

use bls12_381::{G1Affine, G1Projective, Scalar};
use rayon::prelude::*;

const SCALAR_BITS: usize = 255;
const NAIVE_THRESHOLD: usize = 64;

// ============================================================================
// Bit Extraction
// ============================================================================

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
// Helper Functions
// ============================================================================

fn optimal_chunk_size(n: usize) -> usize {
    if n < 32 { 3 } else { (n as f64).ln().ceil() as usize }
}

fn has_nonzero_bits(bytes: &[u8; 32], start_bit: usize, num_bits: usize) -> bool {
    for i in 0..num_bits {
        let bit_pos = start_bit + i;
        if bit_pos < 256 {
            let byte_idx = bit_pos / 8;
            let bit_idx = bit_pos % 8;
            if (bytes[byte_idx] >> bit_idx) & 1 != 0 {
                return true;
            }
        }
    }
    false
}

// ============================================================================
// Bellman-style multiexp (SERIAL)
// ============================================================================

/// Bellman-style multiexp using ln(n) chunks and summation by parts
pub fn bellman_style_multiexp(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= NAIVE_THRESHOLD { return naive_msm_stack(bases, scalars); }
    
    let c = optimal_chunk_size(n);
    let num_chunks = (SCALAR_BITS + c - 1) / c;
    let bucket_count = 1usize << c;
    
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    let mut chunk_results: Vec<G1Projective> = vec![G1Projective::identity(); num_chunks];
    
    for chunk_idx in 0..num_chunks {
        let bit_pos = chunk_idx * c;
        
        // Skip empty windows
        let mut has_nonzero = false;
        for bytes in &scalar_bytes {
            if has_nonzero_bits(bytes, bit_pos, c) {
                has_nonzero = true;
                break;
            }
        }
        if !has_nonzero { continue; }
        
        let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
        
        for i in 0..n {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, c);
            if k > 0 {
                buckets[k] += bases[i];
            }
        }
        
        // Summation by parts: O(2^c)
        let mut running_sum = G1Projective::identity();
        for bucket in buckets.into_iter().rev() {
            running_sum += bucket;
            chunk_results[chunk_idx] += running_sum;
        }
    }
    
    // Combine chunks: double c times for each chunk position
    let mut result = G1Projective::identity();
    for chunk_idx in 0..num_chunks {
        for _ in 0..c {
            result = result.double();
        }
        result += chunk_results[chunk_idx];
    }
    
    result
}

// ============================================================================
// Bellman-style multiexp (PARALLEL over windows)
// ============================================================================

/// Parallel version: process each window in parallel, then combine
pub fn bellman_style_multiexp_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    if n <= NAIVE_THRESHOLD { return naive_msm_stack(bases, scalars); }
    
    let c = optimal_chunk_size(n);
    let num_chunks = (SCALAR_BITS + c - 1) / c;
    let bucket_count = 1usize << c;
    
    let bases_vec = bases.to_vec();
    let scalar_bytes: Vec<[u8; 32]> = scalars.iter().map(|s| s.to_bytes()).collect();
    
    // Parallelize over chunks
    let chunk_results: Vec<G1Projective> = (0..num_chunks).into_par_iter().map(|chunk_idx| {
        let bit_pos = chunk_idx * c;
        
        // Skip empty windows
        let mut has_nonzero = false;
        for bytes in &scalar_bytes {
            if has_nonzero_bits(bytes, bit_pos, c) {
                has_nonzero = true;
                break;
            }
        }
        if !has_nonzero { return G1Projective::identity(); }
        
        let mut buckets: Vec<G1Projective> = vec![G1Projective::identity(); bucket_count];
        
        for i in 0..n {
            let k = extract_window_bits(&scalar_bytes[i], bit_pos, c);
            if k > 0 {
                buckets[k] += bases_vec[i];
            }
        }
        
        // Summation by parts
        let mut chunk_sum = G1Projective::identity();
        for bucket in buckets.into_iter().rev() {
            chunk_sum += bucket;
        }
        
        chunk_sum
    }).collect();
    
    // Combine chunks: double c times for each chunk position
    let mut result = G1Projective::identity();
    for chunk_idx in 0..num_chunks {
        for _ in 0..c {
            result = result.double();
        }
        result += chunk_results[chunk_idx];
    }
    
    result
}

// ============================================================================
// Naive MSM
// ============================================================================

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

// Aliases for compatibility
pub fn pippenger_serial(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

pub fn pippenger_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

pub fn pippenger_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

// ============================================================================
// Auto-select MSM
// ============================================================================

pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    let n = bases.len();
    if n == 0 { return G1Projective::identity(); }
    
    // Use parallel for large n, serial for small-medium
    if n <= NAIVE_THRESHOLD {
        naive_msm_stack(bases, scalars)
    } else if n <= 2048 {
        bellman_style_multiexp(bases, scalars)  // Serial for n <= 2048
    } else {
        bellman_style_multiexp_parallel(bases, scalars)  // Parallel for n > 2048
    }
}

pub fn sliding_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

pub fn sliding_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

pub fn strauss_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

pub fn strauss_msm_parallel(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    bellman_style_multiexp(bases, scalars)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm { Naive, Strauss, Pippenger }

pub fn msm_with_algorithm(bases: &[G1Affine], scalars: &[Scalar], algorithm: Algorithm) -> G1Projective {
    match algorithm {
        Algorithm::Naive => naive_msm_stack(bases, scalars),
        Algorithm::Strauss | Algorithm::Pippenger => bellman_style_multiexp(bases, scalars),
    }
}

pub fn optimal_algorithm(n: usize) -> Algorithm {
    if n <= NAIVE_THRESHOLD { Algorithm::Naive } else { Algorithm::Pippenger }
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
            let bellman_result = bellman_style_multiexp(&bases, &scalars);
            let auto_result = auto_msm(&bases, &scalars);
            
            assert_eq!(naive_result, bellman_result, "n={}: bellman mismatch", n);
            assert_eq!(naive_result, auto_result, "n={}: auto mismatch", n);
        }
    }
}