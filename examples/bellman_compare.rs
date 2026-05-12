//! Real Bellman vs Our Implementation Benchmark
//! 
//! Run with: cargo run --release --example bellman_compare

use std::time::Instant;
use std::sync::Arc;
use bls12_381::{G1Affine, G1Projective, Scalar};
use cpu_msm_optimized::{auto_msm, bellman_style_multiexp, naive_msm_stack, pippenger_msm_parallel};

use bellman::multiexp::{multiexp, FullDensity};
use bellman::multicore::Worker;

fn main() {
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║          Bellman vs CPU-MSM-Optimized Benchmark (REAL)                    ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝\n");
    
    let sizes = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384];
    
    println!("┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐");
    println!("│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │");
    println!("├─────────┼────────────────┼────────────────┼────────────┼─────────────┤");
    
    let worker = Worker::new();
    
    for &n in &sizes {
        let g = G1Affine::generator();
        let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
        let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
        
        // Benchmark Bellman (real)
        let bellman_time = {
            let start = Instant::now();
            for _ in 0..3 {
                let bases_arc = Arc::new(bases.clone());
                let scalars_exp: Vec<_> = scalars.iter()
                    .map(|s| bellman::multiexp::Exponent::from(s))
                    .collect();
                let scalars_arc = Arc::new(scalars_exp);
                let _: Result<G1Projective, _> = multiexp(&worker, (bases_arc.clone(), 0), FullDensity, scalars_arc).wait();
            }
            start.elapsed().as_secs_f64() * 1000.0 / 3.0
        };
        
        // Benchmark our optimized (Bellman-style ln(n) chunks)
        let optimized_time = {
            let start = Instant::now();
            for _ in 0..3 {
                let _ = auto_msm(&bases, &scalars);
            }
            start.elapsed().as_secs_f64() * 1000.0 / 3.0
        };
        
        let speedup = bellman_time / optimized_time;
        let winner = if speedup > 1.0 { "OPT" } else { "BNM" };
        
        println!("│ {:>7} │ {:>14.2}ms │ {:>14.2}ms │ {:>10.2}x │ {:>11} │",
                 n, bellman_time, optimized_time, speedup, winner);
    }
    
    println!("└─────────┴────────────────┴────────────────┴────────────┴─────────────┘\n");
    
    // Algorithm breakdown
    println!("=== Algorithm Breakdown ===\n");
    
    let sizes = [64, 256, 1024, 4096, 16384];
    println!("┌─────────┬────────────┬────────────┬────────────┬────────────┐");
    println!("│   n     │    Naive   │  Bellman   │  Ours ln(n) │  Winner    │");
    println!("├─────────┼────────────┼────────────┼────────────┼────────────┤");
    
    for &n in &sizes {
        let g = G1Affine::generator();
        let bases: Vec<G1Affine> = (0..n).map(|_| g).collect();
        let scalars: Vec<Scalar> = (0..n).map(|i| Scalar::from(i as u64 + 1)).collect();
        
        let naive_time = {
            let start = Instant::now();
            for _ in 0..3 { let _ = naive_msm_stack(&bases, &scalars); }
            start.elapsed().as_secs_f64() * 1000.0 / 3.0
        };
        
        let bellman_time = {
            let start = Instant::now();
            for _ in 0..3 {
                let bases_arc = Arc::new(bases.clone());
                let scalars_exp: Vec<_> = scalars.iter()
                    .map(|s| bellman::multiexp::Exponent::from(s))
                    .collect();
                let scalars_arc = Arc::new(scalars_exp);
                let _: Result<G1Projective, _> = multiexp(&worker, (bases_arc.clone(), 0), FullDensity, scalars_arc).wait();
            }
            start.elapsed().as_secs_f64() * 1000.0 / 3.0
        };
        
        let ours_time = {
            let start = Instant::now();
            for _ in 0..3 { let _ = bellman_style_multiexp(&bases, &scalars); }
            start.elapsed().as_secs_f64() * 1000.0 / 3.0
        };
        
        let winner = if bellman_time < ours_time && bellman_time < naive_time {
            "Bellman"
        } else if ours_time < naive_time {
            "Ours"
        } else {
            "Naive"
        };
        
        println!("│ {:>7} │ {:>10.2}ms │ {:>10.2}ms │ {:>10.2}ms │ {:>10} │",
                 n, naive_time, bellman_time, ours_time, winner);
    }
    
    println!("└─────────┴────────────┴────────────┴────────────┴────────────┘");
    println!("\nNote: Bellman times are REAL (using actual bellman crate)");
    println!("      Our times use ln(n) chunks + summation by parts");
}