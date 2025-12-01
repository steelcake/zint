# zint

Integer compression library for Zig based on FastLanes

## Features

- High performance. Zint uses FastLanes internally so it implements fully vectorized delta encoding and bitpacking.
- Safe decompress API. It is safe to decompress untrusted input.
- Compresses each block (1024) elements inside the input dynamically. Automatically applies delta encoding if it is beneficial.

## Example usage

```zig
const Z = @require("zint").Zint(T);

const compress_buf = try std.heap.page_allocator.alloc(u8, Z.compress_bound(input.len));

const compressed_size = Z.compress(input, compress_buf);

const compressed = compress_buf[0..compressed_size];

const output = try std.heap.page_allocator.alloc(T, input.len);
try Z.decompress(compressed, output);

std.debug.assert(std.mem.eql(T, output, input));
```

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


