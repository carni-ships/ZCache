# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance Results vs Bellman (Real Benchmark)

```
╔══════════════════════════════════════════════════════════════════════════╗
║          Bellman vs CPU-MSM-Optimized Benchmark (Real)                  ║
╚══════════════════════════════════════════════════════════════════════════╝

┌─────────┬────────────────┬────────────────┬────────────┬─────────────┐
│   n     │    Bellman     │   Optimized    │  Speedup   │  Winner     │
├─────────┼────────────────┼────────────────┼────────────┼─────────────┤
│      32 │           0.93ms│           7.58ms│       0.12x│          BNM│
│      64 │           0.62ms│          14.89ms│       0.04x│          BNM│
│     128 │           0.54ms│           0.33ms│       1.65x│         OPT│
│     256 │           0.89ms│           0.50ms│       1.77x│         OPT│
│     512 │           1.28ms│           0.90ms│       1.43x│         OPT│
│    1024 │           1.66ms│           1.42ms│       1.17x│         OPT│
│    2048 │           3.30ms│           2.73ms│       1.21x│         OPT│
│    4096 │           4.87ms│           5.42ms│       0.90x│          BNM│
│    8192 │           9.51ms│          10.67ms│       0.89x│          BNM│
│   16384 │          14.78ms│          19.44ms│       0.76x│          BNM│
│   32768 │          32.15ms│          38.68ms│       0.83x│          BNM│
│   65536 │          59.53ms│          77.13ms│       0.77x│          BNM│
│  131072 │         111.47ms│         147.61ms│       0.76x│          BNM│
│  262144 │         227.47ms│         300.85ms│       0.76x│          BNM│
│  524288 │         498.12ms│         595.21ms│       0.84x│          BNM│
│ 1048576 │         861.32ms│        1159.70ms│       0.74x│          BNM│
└─────────┴────────────────┴────────────────┴────────────┴─────────────┘
```

## Key Results

| Range | Winner | Speedup |
|-------|--------|---------|
| n = 128 - 2048 | **Ours** | **1.17x - 1.77x faster** |
| n = 32-64 | Bellman (naive region) | - |
| n ≥ 4096 | Bellman | 1.2x - 1.4x slower |

### At Scale (2^16 = 65,536 points)

| Implementation | Time | 
|----------------|------|
| Bellman (parallel) | 59.5ms |
| Ours (serial) | 77.1ms |
| Ratio | 0.77x (Bellman wins) |

## Algorithm Design

### Serial Implementation (ln(n) Chunks)

```
bellman_style_multiexp(bases, scalars):
├─ c = ceil(ln(n)) chunk size (3-15 based on n)
├─ num_chunks = 255/c (85 down to 17 windows)
├─ For each chunk (serial):
│  ├─ Allocate 2^c buckets
│  ├─ Accumulate bases[i] into bucket[k]
│  └─ Weighted sum: sum(k * bucket[k])
└─ Combine: Process from MSB to LSB (doubling + adding)

Time: O(n × c + num_chunks × 2^c)
```

## Why We Win at n=128-2048

1. **ln(n) Chunk Size**: Optimal bucket count vs fixed windows
2. **Direct Weighted Sum**: O(2^c) instead of O(c × 2^c)
3. **No Parallelization Overhead**: Serial avoids thread coordination costs

## Why Bellman Wins at n≥4096

Bellman uses **point-parallelization** with:
- Multicore worker pool with work-stealing
- Better load balancing across threads
- Batched memory access patterns
- Production-hardened over years

Our serial implementation is memory-efficient but lacks parallel speedup.

## Parallelization Attempts (Lessons Learned)

### 1. Window-Parallelization
- Process each window in parallel using Rayon
- Problem: Fixed number of windows (17-85) vs 8+ cores
- Result: Poor CPU utilization at large n

### 2. Point-Parallelization
- Each thread processes a range of points
- Problem: Memory explosion (num_threads × num_chunks × bucket_count buckets)
- For n=1M: 8 × 19 × 16384 = 2.6M buckets = 125MB just for buckets!
- Result: Memory bandwidth saturation

### 3. Window-Batching
- Process 4-8 windows at a time in parallel
- Problem: Complex scaling logic was error-prone
- Result: Correctness bugs in chunk combination

### Current Strategy
- Pure serial for all n (memory-efficient, correct)
- Future work: Hybrid approach with proper memory management

## Benchmark Methodology

Real benchmark using the actual `bellman` crate:
```bash
cargo run --release --example bellman_compare
```

- Generate random bases and scalars
- Run 5 iterations, take average
- Warmup iterations ensure JIT-free measurement

## Key Files

- `src/lib.rs` - Main MSM implementation (190 lines, pure serial)
- `examples/bellman_compare.rs` - Real Bellman comparison benchmark
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## Running Benchmarks

```bash
# Run Bellman comparison
cargo run --release --example bellman_compare

# Run tests
cargo test --release

# Run specific test with output
cargo test --release test_correctness -- --nocapture
```

## License

MIT or Apache 2.0