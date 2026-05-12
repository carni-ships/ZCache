# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance Results (Real Benchmark vs Bellman)

```
╔══════════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized Benchmark (REAL)                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      32 │           0.92ms │          11.79ms │       0.08x │         BNM │
│      64 │           0.60ms │          15.36ms │       0.04x │         BNM │
│     128 │           0.53ms │           0.31ms │       1.74x │         OPT │
│     256 │           0.88ms │           0.46ms │       1.93x │         OPT │
│     512 │           1.12ms │           0.81ms │       1.38x │         OPT │
│    1024 │           1.62ms │           1.37ms │       1.19x │         OPT │
│    2048 │           3.05ms │           2.63ms │       1.16x │         OPT │
│    4096 │           4.57ms │           4.98ms │       0.92x │         BNM │
│    8192 │           7.42ms │           9.83ms │       0.75x │         BNM │
│   16384 │          13.47ms │          18.85ms │       0.71x │         BNM │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Key Results

| Points | Bellman (real) | Ours ln(n) | Winner | Speedup |
|--------|----------------|------------|--------|---------|
| 32 | 0.92ms | 11.79ms | Bellman | - |
| 64 | 0.60ms | 15.36ms | Bellman | - |
| **128** | 0.53ms | 0.31ms | **Ours** | **1.74x** |
| **256** | 0.88ms | 0.46ms | **Ours** | **1.93x** |
| **512** | 1.12ms | 0.81ms | **Ours** | **1.38x** |
| **1024** | 1.62ms | 1.37ms | **Ours** | **1.19x** |
| **2048** | 3.05ms | 2.63ms | **Ours** | **1.16x** |
| 4096 | 4.57ms | 4.98ms | Bellman | - |
| 8192 | 7.42ms | 9.83ms | Bellman | - |
| 16384 | 13.47ms | 18.85ms | Bellman | - |

## Algorithm Selection

The implementation uses adaptive algorithm selection based on input size:

| Input Size (n) | Algorithm | Reason |
|----------------|-----------|--------|
| n ≤ 64 | Naive | No bucket overhead, fastest at small sizes |
| n > 64 | Bellman-style ln(n) | Optimal chunk size for medium/large inputs |

### Why We Beat Bellman at n=128-2048

1. **ln(n) Chunk Size**: Bellman uses `ceil(ln(n))` for optimal bucket count
2. **Summation by Parts**: O(2^c) reduction instead of O(c * 2^c)
3. **Density Tracking**: Skip windows where no scalar has non-zero bits

### Why Bellman Beats Us at n>4096

1. **Better Parallelization**: Bellman uses multicore Worker with better task distribution
2. **Batch Processing**: Bellman processes bases in batches with better cache locality
3. **Optimization Maturity**: Bellman is battle-tested Zcash code

## Implemented Optimizations

1. **Adaptive Algorithm Selection** - Naive for n≤64, Bellman-style ln(n) for n>64
2. **ln(n) Chunk Size** - Optimal window size based on input size
3. **Summation by Parts** - O(2^c) reduction instead of O(c * 2^c)
4. **Density Tracking** - Skip empty windows
5. **Batch Scalar Conversion** - Convert scalars to bytes once

## Running Benchmarks

```bash
# Run real Bellman comparison
cargo run --release --example bellman_compare

# Run tests
cargo test

# Run performance tests
cargo test test_performance -- --nocapture
```

## Key Files

- `src/lib.rs` - Main MSM implementation (naive + Bellman-style)
- `examples/bellman_compare.rs` - Real Bellman comparison benchmark
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## Architecture

```
auto_msm(bases, scalars):
├─ n ≤ 64:     naive_msm_stack(bases, scalars)
│              → O(n) direct scalar multiplication
│              → No bucket allocation overhead
│
└─ n > 64:     bellman_style_multiexp(bases, scalars)
               → c = ceil(ln(n)) chunk size
               → 2^c buckets per chunk
               → Summation by parts reduction
               → Double c times for each chunk position

Bellman-style Formula:
result = Σ (2^c)^i * (Σ (2^(c-1) * bucket + ... + bucket[0]))
```

## References

- Pippenger, N. (1979). On the evaluation of powers and logarithms.
- Bellman Multiexp: https://github.com/zkcrypto/bellman
- Zcash Sapling Protocol Specification

## License

MIT or Apache 2.0 (same as Zcash)