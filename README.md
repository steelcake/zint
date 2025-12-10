# zint

Integer compression library for Zig based on FastLanes

## Features

- Zint uses FastLanes internally so it implements fully vectorized delta encoding and bitpacking.
- It is safe to decompress untrusted input.
- Supports signed integers.
- Supports 128 and 256 bit integers.
- Supports delta, FrameOfReference and regular bitpacking.

## Example usage

```zig
const zint = @import("zint");
const Zint = zint.Zint;
const Ctx = zint.Ctx;

const ctx = try Ctx.init(std.heap.page_allocator);
defer ctx.deinit(std.heap.page_allocator);

const compress_buf = try std.heap.page_allocator.alloc(u8, Z.bitpack_compress_bound(input.len));

const compressed_size = Z.bitpack_compress(ctx, input, compress_buf);

const compressed = compress_buf[0..compressed_size];

const output = try std.heap.page_allocator.alloc(T, input.len);
const consumed_size = try Z.decompress(ctx, compressed, output);

std.debug.assert(compressed_size == consumed_size);

std.debug.assert(std.mem.eql(T, output, input));
```

# FastLanes

Zint is based on the FastLanes integer compression.

The fastlanes implementation on this repo is ported from rust implementation of fastlanes by spiraldb.

- [paper](https://www.vldb.org/pvldb/vol16/p2132-afroozeh.pdf)
- [original implementation](https://github.com/cwida/FastLanes) 
- [rust implementation by spiraldb](https://github.com/spiraldb/fastlanes)

# Benchmark Results

Some example results are located in [./benchmark_results](./benchmark_results).

The benchmark measures performance of different Zint variants on some synthetic datasets.

Run with `make bench`.

My interpretation of the benchmarks is that using the right Zint variant for right kind of data is massively faster than lz4 or zstd.
Bitpacking and FrameOfReference+BitPacking almost match memcopy performance but Delta+Bitpacking is slower than memcopy.

In terms of compression ratio, the right variant of Zint is better or on par with zstd(level 1) while lz4 tends to be worse.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.


