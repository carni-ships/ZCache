# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance vs Bellman (Real Benchmark)

```
╔══════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized (v50)                           ║
║          Point-Parallel MSM with Rayon                               ║
╚══════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      64 │           0.61ms │          14.85ms │       0.04x │       BNM │
│     128 │           0.58ms │           2.78ms │       0.21x │       BNM │
│     256 │           0.83ms │           3.04ms │       0.27x │       BNM │
│     512 │           1.38ms │           3.94ms │       0.35x │       BNM │
│    1024 │           1.86ms │           4.33ms │       0.43x │       BNM │
│    2048 │           3.39ms │           6.74ms │       0.50x │       BNM │
│    4096 │           4.75ms │           7.69ms │       0.62x │       BNM │
│    8192 │           9.12ms │          13.92ms │       0.66x │       BNM │
│   16384 │          16.32ms │          14.55ms │       1.12x │  **OPT**  │
│   32768 │          31.08ms │          16.99ms │       1.83x │  **OPT**  │
│   65536 │          61.40ms │          24.36ms │       2.52x │  **OPT**  │
│  131072 │         112.89ms │          41.79ms │       2.70x │  **OPT**  │
│  262144 │         225.50ms │          68.61ms │       3.29x │  **OPT**  │
│  524288 │         446.98ms │         136.35ms │       3.28x │  **OPT**  │
│ 1048576 │         857.13ms │         252.96ms │       3.38x │  **OPT**  │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Algorithm Overview

### Point-Parallelization Strategy

```
Thread 0: points[0..n/threads), ALL windows → returns num_windows partial sums
Thread 1: points[n/threads..2n/threads), ALL windows → returns num_windows partial sums
...
Thread N-1: points[(N-1)n/threads..n), ALL windows → returns num_windows partial sums

Combine: sum thread results per window, then windows MSB to LSB
```

Each thread has its own buckets (no race conditions), and returns ONE value per window.

### Pippenger Algorithm

1. **Window decomposition**: Split scalar into c-bit windows (w = 4-7)
2. **Bucket accumulation**: bases[i] → buckets[k] where k = bits_w(scalars[i])
3. **Summation by parts**: sum(k × bucket[k]) = Σ Σ bucket[j] for j >= k
4. **Window combination**: Process windows MSB to LSB with doubling

### Adaptive Window Size

| n Range | Window Size (w) | Bucket Count (2^w) |
|---------|----------------|-------------------|
| n <= 256 | 4 | 16 |
| n <= 1024 | 5 | 32 |
| n <= 4096 | 6 | 64 |
| n > 4096 | 7 | 128 |

## Why We Win at Large N

- **Better parallelization** - each thread processes ALL windows for its points
- **Thread-local buckets** - no mutex contention
- **Efficient bucket reduction** - summation by parts: O(2^c) instead of O(c × 2^c)
- **Scales better** - asymptotics favor point-parallelization at large n

## Why Bellman Wins at Small-Medium N

- **8-core parallelization from the start** - even for n=64
- **Work-stealing thread pool** - efficient load balancing
- **Optimized C code** - years of production optimization

## Optimization Attempts

| Optimization | Effort | Result |
|-------------|--------|--------|
| NAF Encoding | Medium | ❌ Had bugs |
| AVX2 Bit Extraction | Medium | ❌ No gain (point ops still scalar) |
| AVX2 Point Operations | Very High | ❌ Not implemented (requires field intrinsics) |
| Window-Parallel | Low | ❌ Slow (too few parallel tasks) |
| **Point-Parallel** | Medium | ✅ **Success** |

## Running

```bash
# Compare against real Bellman
cargo run --release --example bellman_compare

# Run tests
cargo test --release

# Run benchmarks
cargo run --release --example benchmark
```

## Files

- `src/lib.rs` - Main implementation (~210 lines)
- `examples/bellman_compare.rs` - Real Bellman comparison
- `examples/benchmark.rs` - Benchmark comparison
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## History

| Version | Key Change |
|---------|------------|
| v50 | Point-parallel with thread-local buckets, beats Bellman at n >= 16384 |
| v48-49 | SIMD attempts (not fully working) |
| v47 | Bottleneck analysis: 95% time in bucket accumulation |
| v45 | Point-parallelization (thread-local buckets) |
| v41 | Hybrid serial/parallel |
| v39 | Honest assessment - serial loses to Bellman |
| v38 | Incorrect "3.5x faster" (measurement error) |

## License

MIT or Apache 2.0