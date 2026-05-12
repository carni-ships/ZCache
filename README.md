# CPU MSM Optimization

Optimized Multi-Scalar Multiplication (MSM) for BLS12-381 G1 curves, written in pure Rust.

## Performance vs Bellman

| n | Bellman | Ours | Winner | Speedup |
|---|---------|------|--------|---------|
| 64 | 0.61ms | 14.85ms | Bellman | 0.04x |
| 128 | 0.58ms | 2.78ms | Bellman | 0.21x |
| 256 | 0.83ms | 3.04ms | Bellman | 0.27x |
| 512 | 1.38ms | 3.94ms | Bellman | 0.35x |
| 1024 | 1.86ms | 4.33ms | Bellman | 0.43x |
| 2048 | 3.39ms | 6.74ms | Bellman | 0.50x |
| 4096 | 4.75ms | 7.69ms | Bellman | 0.62x |
| 8192 | 9.12ms | 13.92ms | Bellman | 0.66x |
| 16384 | 16.32ms | 14.55ms | **Ours** | 1.12x |
| 32768 | 31.08ms | 16.99ms | **Ours** | 1.83x |
| 65536 | 61.40ms | 24.36ms | **Ours** | 2.52x |
| 131072 | 112.89ms | 41.79ms | **Ours** | 2.70x |
| 262144 | 225.50ms | 68.61ms | **Ours** | 3.29x |
| 524288 | 446.98ms | 136.35ms | **Ours** | 3.28x |
| 1048576 | 857.13ms | 252.96ms | **Ours** | 3.38x |

- **n <= 8192**: Bellman wins (8-core parallelization effective)
- **n >= 16384**: Ours wins (point-parallelization scales better)

## Algorithm

### Point-Parallelization

```
Thread 0: points[0..n/threads), ALL windows -> returns num_windows partial sums
Thread 1: points[n/threads..2n/threads), ALL windows -> returns num_windows partial sums
...
Thread N-1: points[(N-1)n/threads..n), ALL windows -> returns num_windows partial sums

Combine: sum thread results per window, then windows MSB to LSB
```

Each thread has its own buckets (no mutex/race conditions).

### Pippenger Algorithm

1. Window decomposition: Split scalar into c-bit windows (w = 4-7)
2. Bucket accumulation: bases[i] -> buckets[k] where k = bits_w(scalars[i])
3. Summation by parts: sum(k * bucket[k]) = Sum_j>=k bucket[j]
4. Window combination: Process windows MSB to LSB with doubling

### Adaptive Window Size

| n Range | w | Buckets (2^w) |
|---------|---|---------------|
| n <= 256 | 4 | 16 |
| n <= 1024 | 5 | 32 |
| n <= 4096 | 6 | 64 |
| n > 4096 | 7 | 128 |

## Why We Win at Large N

- Better parallelization (each thread processes ALL windows for its points)
- Thread-local buckets (no mutex contention)
- Summation by parts: O(2^c) reduction

## Why Bellman Wins at Small N

- 8-core parallelization from the start
- Optimized C code (years of production optimization)

## Optimization Attempts

| Optimization | Result |
|-------------|--------|
| NAF Encoding | Buggy |
| AVX2 Bit Extraction | No gain |
| AVX2 Point Operations | Not implemented |
| Window-Parallel | Too few tasks |
| **Point-Parallel** | **Success** |

## Usage

```bash
# Compare against real Bellman
cargo run --release --example bellman_compare

# Run tests
cargo test --release
```

## Files

- `src/lib.rs` - Main implementation (~210 lines)
- `examples/bellman_compare.rs` - Bellman comparison
- `Cargo.toml` - Dependencies (bls12_381, bellman, rayon)

## History

- v50: Point-parallel, beats Bellman at n >= 16384
- v48-49: SIMD attempts (not working)
- v47: Bottleneck analysis: 95% in bucket accumulation
- v45: Point-parallelization with thread-local buckets
- v39: Honest assessment
- v38: Incorrect "3.5x faster" (measurement error)

## License

MIT or Apache 2.0