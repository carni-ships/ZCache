# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance Results (Real Benchmark vs Bellman, up to 2^20)

```
╔══════════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized Benchmark (REAL)                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      32 │           0.53ms │           7.51ms │       0.07x │         BNM │
│      64 │           0.62ms │          14.95ms │       0.04x │         BNM │
│     128 │           0.57ms │           0.32ms │       1.80x │         OPT │
│     256 │           0.91ms │           0.50ms │       1.81x │         OPT │
│     512 │           1.32ms │           0.93ms │       1.42x │         OPT │
│    1024 │           1.68ms │           1.46ms │       1.15x │         OPT │
│    2048 │           3.19ms │           2.75ms │       1.16x │         OPT │
│    4096 │           4.94ms │           5.41ms │       0.91x │         BNM │
│    8192 │           9.50ms │          10.69ms │       0.89x │         BNM │
│   16384 │          14.74ms │          19.54ms │       0.75x │         BNM │
│   32768 │          32.22ms │          38.79ms │       0.83x │         BNM │
│   65536 │          59.41ms │          77.15ms │       0.77x │         BNM │
│  131072 │         111.53ms │         147.73ms │       0.75x │         BNM │
│  262144 │         227.49ms │         300.86ms │       0.76x │         BNM │
│  524288 │         498.07ms │         595.25ms │       0.84x │         BNM │
│ 1048576 │         861.25ms │        1159.68ms │       0.74x │         BNM │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Key Results

| Range | Winner | Speedup |
|-------|--------|---------|
| n = 128 - 2048 | **Ours** | **1.15x - 1.81x faster** |
| n = 32-64 | Bellman (naive region) | - |
| n ≥ 4096 | Bellman (parallel) | 1.2x - 1.4x slower |

### At Scale (2^20 = 1,048,576 points)

| Implementation | Time | 
|----------------|------|
| Bellman (parallel) | 861ms |
| Ours (serial) | 1160ms |
| Ratio | 0.74x (Bellman wins) |

## Algorithm Selection

| Input Size (n) | Algorithm | Reason |
|----------------|-----------|--------|
| n ≤ 64 | Naive | No bucket overhead, O(n) direct multiplication |
| n > 64 | Bellman-style serial | ln(n) chunks + summation by parts |

### Implementation

```rust
pub fn auto_msm(bases, scalars) {
    if n <= 64 { naive(bases, scalars) }      // O(n), no buckets
    else { bellman_style_multiexp(bases, scalars) }  // ln(n) chunks
}
```

## Why We Win at n=128-2048

1. **ln(n) Chunk Size**: Optimal bucket count based on input size
2. **Summation by Parts**: O(2^c) reduction instead of O(c × 2^c)
3. **Density Tracking**: Skip windows with no non-zero bits
4. **Pure Rust**: No C/assembly, optimized LLVM codegen

## Why Bellman Wins at n≥4096

1. **Multicore Parallelization**: Bellman uses parallel over points
2. **Worker Thread Pool**: Optimized work distribution
3. **Production Hardened**: Years of Zcash optimization

## Benchmark Methodology

Real benchmark using the actual `bellman` crate (not estimates):
- Generate random bases and scalars
- Run 5 iterations, take average
- Warmup iterations to ensure JIT-free measurement

## Running Benchmarks

```bash
# Run real Bellman comparison (n up to 2^20)
cargo run --release --example bellman_compare

# Run tests
cargo test
```

## Architecture

```
bellman_style_multiexp(bases, scalars):
├─ c = ceil(ln(n)) chunk size (4-8 based on n)
├─ For each chunk (51-85 windows):
│  ├─ Skip if no non-zero bits (density tracking)
│  ├─ Allocate 2^c buckets
│  ├─ Accumulate bases[i] into bucket[k]
│  └─ Summation by parts: O(2^c)
└─ Combine chunks: double c times for each, then add

Time: O(n × c + 2^c × num_chunks)
Space: O(2^c) buckets
```

## Comparison: Algorithm Breakdown

```
┌─────────┬────────────┬────────────┬────────────┬────────────┐
│   n     │    Naive   │  Bellman   │  Ours ln(n) │  Winner    │
├─────────┼────────────┼────────────┼────────────┼────────────┤
│      64 │      14.78ms │       0.60ms │      14.81ms │    Bellman │
│     256 │      58.69ms │       0.85ms │       0.51ms │       Ours │
│    1024 │     235.85ms │       1.70ms │       1.44ms │       Ours │
│    4096 │     951.61ms │       4.92ms │       5.43ms │    Bellman │
│   16384 │    3791.85ms │      16.49ms │      19.43ms │    Bellman │
└─────────┴────────────┴────────────┴────────────┴────────────┘
```

## Key Files

- `src/lib.rs` - Main MSM implementation
- `examples/bellman_compare.rs` - Real Bellman comparison (up to 2^20)
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## References

- Pippenger, N. (1979). On the evaluation of powers and logarithms
- Bellman Multiexp: https://github.com/zkcrypto/bellman
- Zcash Sapling Protocol Specification

## License

MIT or Apache 2.0