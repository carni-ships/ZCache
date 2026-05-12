# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Current Performance vs Bellman (Real Benchmark)

```
╔══════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized (v41)                           ║
╚══════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      64 │           0.53ms │          14.92ms │       0.04x │        BNM │
│     128 │           0.54ms │           1.85ms │       0.29x │        BNM │
│     256 │           0.85ms │           3.21ms │       0.26x │        BNM │
│     512 │           1.15ms │           5.60ms │       0.21x │        BNM │
│    1024 │           1.66ms │           6.15ms │       0.27x │        BNM │
│    4096 │           4.41ms │          19.70ms │       0.22x │        BNM │
│   16384 │          15.11ms │          44.43ms │       0.34x │        BNM │
│   65536 │          61.86ms │         187.87ms │       0.33x │        BNM │
│ 1048576 │         844.66ms │         725.81ms │       1.16x │  **OPT**   │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Analysis: Why We Lose

### Bellman's Advantage

| Factor | Bellman | Our Implementation |
|--------|---------|-------------------|
| Cores | 8-core parallel | Serial (only parallel at n≥65536) |
| Scheduling | Work-stealing | Rayon thread pool |
| Memory | Batched, cache-friendly | Sequential |

### Our Serial Algorithm

The bucket method (ln(n) chunks, summation by parts) is correct and 10-20x faster than naive for medium n:
- n=64: naive=14.6ms, bucket=14.5ms (same, overhead dominates)
- n=256: naive=58.4ms, bucket=3.3ms (17x faster!)
- n=1024: naive=236ms, bucket=6.2ms (38x faster!)

But Bellman's **8-core parallelization** gives them 4-25x speedup over our serial code.

### What We Tried

1. **Window-parallelization** - Parallel over windows (19-51 tasks) - too few for 8 cores
2. **Point-parallelization** - Parallel over points - memory explosion (2.5M buckets)
3. **Window-sum reduction** - Thread returns one value per window - still too slow
4. **Hybrid (serial + parallel for n≥65536)** - Parallel only helps at 1M+ points

### Why Parallelization Didn't Help

The fundamental issue: each parallel task still does O(n) work, and the thread coordination overhead (Vec allocation, thread spawning, result combining) dominates at small-medium n.

## Algorithm Implementation

### Serial Bellman-style MSM (ln(n) chunks)

```
serial_msm(bases, scalars):
├─ c = ceil(ln(n)) chunk size (3-15 based on n)
├─ num_chunks = 255/c
├─ For each chunk:
│  ├─ Allocate 2^c buckets
│  ├─ Accumulate bases[i] into bucket[k]
│  └─ Summation by parts: sum(k * bucket[k])
└─ Combine: MSB to LSB with doubling
```

### Parallel MSM (hybrid, for n≥65536)

```
parallel_msm(bases, scalars):
├─ Split points across threads (each thread = points/n_threads)
├─ Each thread processes ALL windows sequentially:
│  ├─ Allocate 2^c buckets (reuse for each window)
│  ├─ Accumulate thread's points into buckets
│  └─ Reduce to one value per window
├─ Return: num_chunks values (one per window)
└─ Combine: just num_chunks additions + MSB-to-LSB doubling
```

## To Beat Bellman

**Priority order:**

1. **SIMD bit extraction** - AVX2 to extract bits from 8 scalars simultaneously
2. **Better parallelization** - Process points in cache-friendly batches
3. **Atomic bucket accumulation** - Shared buckets with atomics (avoids per-thread bucket allocation)
4. **Assembly intrinsics** - x86_64 optimized point operations

The serial implementation has reached its limit. Further improvements require **proper SIMD + parallelization working together**.

## Running

```bash
# Compare against real Bellman
cargo run --release --example bellman_compare

# Run tests
cargo test --release
```

## Files

- `src/lib.rs` - Main implementation (261 lines)
- `examples/bellman_compare.rs` - Real Bellman comparison
- `Cargo.toml` - Dependencies

## History

- v41: Hybrid serial/parallel (parallel only at n≥65536)
- v40: Point-parallel with window-sum reduction
- v39: Clean serial implementation
- v38: Summation by parts bug fixed
- v37: Lessons learned on parallelization attempts

## License

MIT or Apache 2.0