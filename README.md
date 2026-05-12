# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Current Performance vs Bellman (Real Benchmark)

```
╔══════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized (Current State)                    ║
╚══════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│     128 │           0.64ms│           1.98ms│       0.32x│        BNM  │
│     256 │           0.85ms│           3.39ms│       0.25x│        BNM  │
│     512 │           1.39ms│           5.68ms│       0.24x│        BNM  │
│    1024 │           1.94ms│           6.23ms│       0.31x│        BNM  │
│    4096 │           4.96ms│          20.50ms│       0.24x│        BNM  │
│   16384 │          16.54ms│          45.76ms│       0.36x│        BNM  │
│   65536 │          61.51ms│         167.94ms│       0.37x│        BNM  │
│ 1048576 │         875.75ms│        1422.14ms│       0.62x│        BNM  │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Analysis: Why We Lose

| Factor | Bellman | Our Implementation |
|--------|---------|-------------------|
| Parallelization | 8-core Worker pool | Serial (single-threaded) |
| Chunk handling | Point-parallel | Serial per-window |
| Memory access | Batched, cache-friendly | Sequential |

**Key insight**: Our serial implementation loses to Bellman's parallel implementation at all sizes.

## Implementation

### Algorithm: ln(n) Chunks + Summation by Parts

```
optimized_msm(bases, scalars):
├─ c = ceil(ln(n)) chunk size
├─ num_chunks = 255/c
├─ bucket_count = 2^c
├─ For each chunk:
│  ├─ Accumulate bases[i] into bucket[k]
│  └─ Summation by parts: sum(k * bucket[k])
└─ Combine: MSB to LSB with doubling

Time: O(n × num_chunks) for accumulation + O(num_chunks × 2^c) for reduction
```

### Summation by Parts Formula

```
sum_{k=1}^{m} k * bucket[k] = sum_{k=1}^{m} running_sum[k]

where running_sum[k] = bucket[k] + bucket[k+1] + ... + bucket[m]

This is computed by iterating backward and accumulating:
running_sum += bucket[k]
result += running_sum
```

## To Improve: Add Proper Parallelization

Priority order:
1. **Point-parallel with shared buckets** - Each thread processes a range of points, accumulate into shared bucket arrays
2. **Window-batch parallelization** - Process 4-8 windows at a time in parallel
3. **SIMD bit extraction** - Extract bits from multiple scalars simultaneously

## Files

- `src/lib.rs` - Main implementation (160 lines)
- `examples/bellman_compare.rs` - Real Bellman comparison
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## Running

```bash
# Compare against real Bellman
cargo run --release --example bellman_compare

# Run tests
cargo test --release
```

## History

- v39: Clean implementation, verified correctness
- v38: Reported 3.5x speedup (later found to be measurement error)
- v37: Lessons learned on parallelization attempts
- v30-v35: Window-parallelization (didn't scale)

## License

MIT or Apache 2.0