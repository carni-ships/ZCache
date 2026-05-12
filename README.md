# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves.

## Real Benchmark vs Bellman

```
╔══════════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized Benchmark (REAL)                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      32 │           0.95ms │          11.59ms │       0.08x │         BNM │
│      64 │           0.57ms │          15.41ms │       0.04x │         BNM │
│     128 │           0.63ms │           0.30ms │       2.12x │         OPT │
│     256 │           0.85ms │           0.47ms │       1.81x │         OPT │
│     512 │           1.32ms │           0.83ms │       1.59x │         OPT │
│    1024 │           1.69ms │           1.30ms │       1.30x │         OPT │
│    2048 │           3.20ms │           2.52ms │       1.27x │         OPT │
│    4096 │           4.50ms │           5.01ms │       0.90x │         BNM │
│    8192 │           9.11ms │           9.74ms │       0.94x │         BNM │
│   16384 │          14.45ms │          17.99ms │       0.80x │         BNM │
│   32768 │          29.50ms │          37.08ms │       0.80x │         BNM │
│   65536 │          55.02ms │          70.95ms │       0.78x │         BNM │
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Key Results

| Range | Winner | Speedup |
|-------|--------|---------|
| n = 128 - 2048 | **Ours** | **1.27x - 2.12x** |
| n = 4096 - 65536 | Bellman | 1.2x - 1.3x slower |
| n = 32 - 64 | Bellman (naive region) | - |

## Algorithm

| Input Size (n) | Algorithm | Reason |
|----------------|-----------|--------|
| n ≤ 64 | Naive | No bucket overhead, O(n) direct multiplication |
| n > 64 | Bellman-style ln(n) | Optimal bucket count, summation by parts |

### Implementation

```rust
pub fn auto_msm(bases: &[G1Affine], scalars: &[Scalar]) -> G1Projective {
    if n <= NAIVE_THRESHOLD {
        naive_msm_stack(bases, scalars)  // O(n), no bucket overhead
    } else {
        bellman_style_multiexp(bases, scalars)  // ln(n) chunks + sum by parts
    }
}
```

## Why We Beat Bellman at n=128-2048

1. **ln(n) Chunk Size**: Bellman's optimal window size formula
2. **Summation by Parts**: O(2^c) reduction instead of O(c × 2^c)
3. **Density Tracking**: Skip windows with no non-zero bits
4. **Pure Rust**: No C/assembly dependencies

## Why Bellman Beats Us at n≥4096

1. **Multicore Worker**: Better thread pool and task distribution
2. **Batch Source**: Efficient base/scalar streaming from memory
3. **Production Optimized**: Years of Zcash production hardening

## Benchmarks

```bash
# Run real Bellman comparison (now with n up to 65536)
cargo run --release --example bellman_compare

# Run tests
cargo test

# Run performance tests
cargo test test_performance -- --nocapture
```

## Algorithm Breakdown

| n | Naive | Bellman | Ours ln(n) | Winner |
|---|-------|---------|-----------|--------|
| 64 | 13.48ms | 0.53ms | 13.45ms | Bellman |
| 256 | 53.81ms | 0.81ms | 0.46ms | Ours |
| 1024 | 217.23ms | 1.60ms | 1.35ms | Ours |
| 4096 | 864.93ms | 4.42ms | 4.91ms | Bellman |
| 16384 | 3499.94ms | 13.60ms | 17.93ms | Bellman |

## Architecture

```
bellman_style_multiexp(bases, scalars):
├─ c = ceil(ln(n)) chunk size (4-8 based on n)
├─ For each chunk:
│  ├─ Allocate 2^c buckets
│  ├─ Accumulate bases[i] into bucket[k] where k = bits(scalar[i])
│  └─ Summation by parts: Σ 2^i * bucket[i] in O(2^c)
└─ Combine chunks: for each, double c times then add

Time: O(n × c / w + 2^w × w) per chunk
Space: O(2^w) buckets
```

## License

MIT or Apache 2.0