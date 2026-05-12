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
│      32 │           0.50ms │           7.50ms │       0.07x │         BNM │
│      64 │           0.68ms │          15.07ms │       0.04x │         BNM │
│     128 │           0.66ms │           0.32ms │       2.07x │         OPT │
│     256 │           0.98ms │           0.51ms │       1.92x │         OPT │
│     512 │           1.53ms │           0.89ms │       1.72x │         OPT │
│    1024 │           1.73ms │           1.44ms │       1.20x │         OPT │
│    2048 │           3.84ms │           2.75ms │       1.39x │         OPT │
│    4096 │           5.59ms │           2.58ms │       2.17x │         OPT │
│    8192 │          10.43ms │           7.89ms │       1.32x │         OPT │
│   16384 │          16.92ms │           9.03ms │       1.87x │         OPT │
│   32768 │          32.84ms │          18.42ms │       1.78x │         OPT │
│   65536 │          64.63ms │          35.25ms │       1.83x │         OPT │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Key Results

| Range | Winner | Speedup |
|-------|--------|---------|
| n = 32-64 | Bellman | naive region |
| **n = 128 - 65536** | **Ours** | **1.2x - 2.2x faster** |

### At Scale

| Points | Bellman (ms) | Ours (ms) | Speedup |
|--------|--------------|-----------|---------|
| 2^10 = 1024 | 1.73 | 1.44 | **1.20x** |
| 2^12 = 4096 | 5.59 | 2.58 | **2.17x** |
| 2^14 = 16384 | 16.92 | 9.03 | **1.87x** |
| 2^16 = 65536 | 64.63 | 35.25 | **1.83x** |

## Algorithm Selection

| Input Size (n) | Algorithm | Reason |
|----------------|-----------|--------|
| n ≤ 64 | Naive | No bucket overhead, O(n) direct multiplication |
| n ≤ 2048 | Bellman-style (serial) | ln(n) chunks, avoids thread overhead |
| n > 2048 | Bellman-style (parallel) | Window-parallel with Rayon |

### Implementation

```rust
pub fn auto_msm(bases, scalars) {
    if n <= 64 { naive(bases, scalars) }           // O(n)
    else if n <= 2048 { bellman_style_serial() }    // ln(n) chunks
    else { bellman_style_parallel() }               // parallel over windows
}
```

## Optimizations

1. **ln(n) Chunk Size**: Optimal window size based on input size
2. **Summation by Parts**: O(2^c) reduction instead of O(c × 2^c)
3. **Density Tracking**: Skip windows with no non-zero bits
4. **Window-Parallel**: Process each chunk in parallel using Rayon
5. **Adaptive Selection**: Use serial for small-medium n, parallel for large n

## Why We Beat Bellman

1. **Window-Parallelization**: Process 51-85 independent windows in parallel
2. **ln(n) Chunk Size**: Optimal bucket count (4-8 bits per window)
3. **Summation by Parts**: Efficient O(2^c) reduction
4. **Pure Rust**: No C/assembly dependencies, optimized LLVM codegen

## Running Benchmarks

```bash
# Run real Bellman comparison (n up to 65536)
cargo run --release --example bellman_compare

# Run tests
cargo test

# Run performance tests
cargo test test_performance -- --nocapture
```

## Architecture

```
bellman_style_multiexp_parallel(bases, scalars):
├─ n > 2048: use parallel version
│
├─ c = ceil(ln(n)) chunk size (4-8 based on n)
├─ num_chunks = 255 / c (51-85 windows)
│
├─ Parallelize over chunks using rayon:
│  └─ For each chunk_idx in parallel:
│     ├─ Allocate 2^c buckets
│     ├─ Accumulate bases[i] into bucket[k]
│     └─ Summation by parts: O(2^c)
│
└─ Combine chunks: for each, double c times then add

Time: O((n × c + 2^c) × num_chunks / num_threads)
```

## Key Files

- `src/lib.rs` - Main MSM implementation
- `examples/bellman_compare.rs` - Real Bellman comparison
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## References

- Pippenger, N. (1979). On the evaluation of powers and logarithms
- Bellman Multiexp: https://github.com/zkcrypto/bellman
- Zcash Sapling Protocol Specification

## License

MIT or Apache 2.0