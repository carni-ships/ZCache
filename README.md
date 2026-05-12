# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance vs Bellman (Real Benchmark)

```
╔══════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized (v45)                           ║
║          Point-Parallel MSM with Rayon                               ║
╚══════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      64 │           0.61ms │          14.85ms │       0.04x │       BNM │
│     128 │           0.59ms │           2.83ms │       0.21x │       BNM │
│     256 │           0.84ms │           2.95ms │       0.28x │       BNM │
│     512 │           1.38ms │           4.18ms │       0.33x │       BNM │
│    1024 │           1.88ms │           4.44ms │       0.42x │       BNM │
│    2048 │           3.44ms │           6.83ms │       0.50x │       BNM │
│    4096 │           4.89ms │           7.73ms │       0.63x │       BNM │
│    8192 │           9.29ms │          14.58ms │       0.64x │       BNM │
│   16384 │          16.35ms │          14.75ms │       1.11x │  **OPT**  │
│   32768 │          33.10ms │          19.04ms │       1.74x │  **OPT**  │
│   65536 │          64.35ms │          26.71ms │       2.41x │  **OPT**  │
│  131072 │         114.70ms │          44.59ms │       2.57x │  **OPT**  │
│  262144 │         228.77ms │          72.12ms │       3.17x │  **OPT**  │
│  524288 │         451.62ms │         153.73ms │       2.94x │  **OPT**  │
│ 1048576 │         859.53ms │         261.59ms │       3.29x │  **OPT**  │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Algorithm Breakdown

| n | Naive | Bellman | Ours | Winner |
|---|-------|---------|------|--------|
| 64 | 14.77ms | 0.61ms | 14.85ms | Bellman |
| 256 | 59.80ms | 1.05ms | 3.10ms | Bellman |
| 1024 | 244.12ms | 1.71ms | 4.35ms | Bellman |
| 4096 | 948.49ms | 4.44ms | 7.58ms | Bellman |
| 16384 | 3815.53ms | 16.62ms | 14.38ms | **Ours** |

## Key Insights

### Why Bellman Wins at Medium Sizes (n < 16384)

- **8-core parallelization** kicks in immediately (even for n=64)
- **Work-stealing thread pool** - efficient load balancing
- **Point-parallelization** - parallelizes over points, not windows

### Why We Win at Large Sizes (n >= 16384)

- **Better point-parallelization** - each thread processes ALL windows for its points
- **Thread-local buckets** - no mutex contention
- **Efficient bucket reduction** - summation by parts: O(2^c) instead of O(c × 2^c)

### Parallelization Strategy

Our point-parallel approach:
```
Thread 0: points[0..n/threads), ALL windows → returns num_windows partial sums
Thread 1: points[n/threads..2n/threads), ALL windows → returns num_windows partial sums
...
Thread N-1: points[(N-1)n/threads..n), ALL windows → returns num_windows partial sums

Combine: sum thread results per window, then windows MSB to LSB
```

Each thread has its own buckets (no race conditions), and returns ONE value per window.

## Implementation

### Key Components

1. **Bit Extraction**: Extract k-bit chunks from 32-byte scalars
2. **Bucket Accumulation**: bases[i] → buckets[k] where k = bits_w(scalars[i])
3. **Summation by Parts**: sum(k × bucket[k]) = Σ Σ bucket[j] for j >= k
4. **Window Combination**: Process windows MSB to LSB with doubling

### Adaptive Window Size

| n Range | Window Size (w) | Bucket Count (2^w) |
|---------|----------------|-------------------|
| n <= 256 | 4 | 16 |
| n <= 1024 | 5 | 32 |
| n <= 4096 | 6 | 64 |
| n > 4096 | 7 | 128 |

## Running

```bash
# Compare against real Bellman
cargo run --release --example bellman_compare

# Run tests (verifies correctness against naive)
cargo test --release

# Run benchmarks
cargo run --release --example benchmark
```

## Files

- `src/lib.rs` - Main implementation (220 lines)
- `examples/bellman_compare.rs` - Real Bellman comparison
- `examples/benchmark.rs` - Benchmark comparison
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## Performance Analysis

| Region | Algorithm | Why |
|--------|-----------|-----|
| n <= 64 | Naive | Bucket overhead not worth it |
| 64 < n < 16384 | Bellman (parallel) | Their 8-core parallelization effective |
| n >= 16384 | Ours (parallel) | Our point-parallelization scales better |

## History

- v45: Point-parallel with thread-local buckets, beats Bellman at n >= 16384
- v43-44: SIMD attempts (not fully implemented)
- v42: Point-parallel with window-sum reduction
- v41: Hybrid serial/parallel
- v40: Point-parallel with window-sum reduction
- v39: Honest assessment - serial loses to Bellman
- v38: Incorrect "3.5x faster" (measurement error)

## License

MIT or Apache 2.0