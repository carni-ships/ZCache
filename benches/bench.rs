//! CPU MSM Performance Benchmark

use bls12_381::{G1Affine, G1Projective, Scalar};
use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use std::time::Instant;
use cpu_msm_optimized::{
    auto_msm, glv_msm, glv_msm_parallel, naive_msm, pippenger_msm, pippenger_msm_parallel,
};

fn generate_points(n: usize) -> (Vec<G1Affine>, Vec<Scalar>) {
    use rand::RngCore;
    let mut rng = rand::thread_rng();
    let g = G1Projective::generator();

    let bases: Vec<G1Affine> = (0..n)
        .map(|_| {
            let mut bytes = [0u8; 64];
            rng.fill_bytes(&mut bytes);
            let s = Scalar::from_bytes_wide(&bytes);
            (g * s).into()
        })
        .collect();

    let scalars: Vec<Scalar> = (0..n)
        .map(|_| {
            let mut bytes = [0u8; 64];
            rng.fill_bytes(&mut bytes);
            Scalar::from_bytes_wide(&bytes)
        })
        .collect();

    (bases, scalars)
}

fn bench_naive(c: &mut Criterion) {
    let mut group = c.benchmark_group("naive");
    for &n in &[32, 64, 128] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(naive_msm(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_pippenger(c: &mut Criterion) {
    let mut group = c.benchmark_group("pippenger");
    for &n in &[64, 256, 512, 1024, 2048, 4096, 8192] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(pippenger_msm(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_pippenger_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("pippenger-parallel");
    for &n in &[512, 1024, 2048, 4096, 8192, 16384] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(pippenger_msm_parallel(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_glv(c: &mut Criterion) {
    let mut group = c.benchmark_group("glv");
    for &n in &[128, 256, 512, 1024, 2048, 4096, 8192] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(glv_msm(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_glv_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("glv-parallel");
    for &n in &[1024, 2048, 4096, 8192, 16384] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(glv_msm_parallel(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_auto(c: &mut Criterion) {
    let mut group = c.benchmark_group("auto");
    for &n in &[32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("n", n), |b| {
            b.iter(|| black_box(auto_msm(&bases, &scalars)))
        });
    }
    group.finish();
}

fn bench_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("throughput");
    for n in [1024, 4096, 16384] {
        let (bases, scalars) = generate_points(n);
        group.bench_function(BenchmarkId::new("pts_per_sec", n), |b| {
            b.iter(|| {
                let start = Instant::now();
                let result = auto_msm(&bases, &scalars);
                let elapsed = start.elapsed();
                let rate = n as f64 / elapsed.as_secs_f64();
                black_box((result, rate))
            })
        });
    }
    group.finish();
}

fn bench_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("comparison");
    let n = 4096;
    let (bases, scalars) = generate_points(n);

    group.bench_function(BenchmarkId::new("naive", n), |b| {
        b.iter(|| black_box(naive_msm(&bases, &scalars)));
    });
    group.bench_function(BenchmarkId::new("pippenger", n), |b| {
        b.iter(|| black_box(pippenger_msm(&bases, &scalars)));
    });
    group.bench_function(BenchmarkId::new("pippenger-parallel", n), |b| {
        b.iter(|| black_box(pippenger_msm_parallel(&bases, &scalars)));
    });
    group.bench_function(BenchmarkId::new("glv", n), |b| {
        b.iter(|| black_box(glv_msm(&bases, &scalars)));
    });
    group.bench_function(BenchmarkId::new("glv-parallel", n), |b| {
        b.iter(|| black_box(glv_msm_parallel(&bases, &scalars)));
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_naive,
    bench_pippenger,
    bench_pippenger_parallel,
    bench_glv,
    bench_glv_parallel,
    bench_auto,
    bench_throughput,
    bench_comparison,
);
criterion_main!(benches);