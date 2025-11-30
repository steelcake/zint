const std = @import("std");

const FastLanes = @import("fastlanes.zig").FastLanes;

pub const Error = error{Invalid};

const Header = packed struct(u8) {
    encoding: enum(u2) {
        noencoding,
        bitpacked,
        deltabitpacked,
        dictbitpacked,
    },
    /// only useful if `encoding != .noencoding`
    num_bits: u6,
};

pub fn Zint(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("Zint only supports u8, u16, u32 and u64"),
    }

    const FL = FastLanes(T);

    return struct {
        fn needed_num_bits(data: [1024]T) u8 {
            var m = data[0];
            for (0..1024) |idx| {
                m = @max(data[idx], m);
            }

            if (m == 0) {
                return 0;
            }

            return @as(u8, std.math.log2_int(T, m)) + 1;
        }

        fn make_dict(input: [1024]T) struct { [1024]T, usize  } {
            var dict: [1024]T = undefined;
            var len: usize = 0;

            for (0..1024) |i| {
                for (0..len) |j| {
                    if (dict[j] == input[i]) {
                        break;
                    }
                } else {
                    dict[len] = input[i];
                    len += 1;
                }
            }

            return .{ dict, len };
        }

        fn dict_encode(dict: [1024]T, input: [1024]T) [1024]u16 {
            var indices: [1024]u16 = undefined;
            for (0..1024) |i| {
                for (0..1024) |j| {
                    if (dict[j] == input[i]) {
                        indices[i] = @intCast(j);
                        break;
                    }
                } else unreachable;
            }
            return indices;
        }

        fn dict_decode(dict: []const align(1) T, indices: [1024]u16) [1024]T {
            var out: [1024]T = undefined;
            for (0..1024) |i| {
                out[i] = dict[indices[i]];
            }
            return out;
        }

        const PACK_BOUND = 1 + @sizeOf(T) * 1024;

        fn pack(input: [1024]T, noalias out: []u8) usize {
            std.debug.assert(out.len >= PACK_BOUND);

            const normal_nbits = needed_num_bits(input);
            const normal_size = @as(usize, normal_nbits) * 1024;

            const dict, const dict_len = make_dict(input);
            const dict_idx_nbits = @as(u8, std.math.log2_int(usize, dict_len)) + 1;

            // +2 is for the dictionary length encoded as u16
            const dict_encoded_size = dict_len * FL.N_BITS + 1024 * @as(usize, dict_idx_nbits) + 2;

            const transposed = FL.transpose(input);
            const bases: [FL.N_LANES]T = transposed[0..FL.N_LANES].*;
            const delta = FL.delta(transposed, bases);

            const delta_nbits = needed_num_bits(delta);
            const delta_encoded_size = FL.N_BITS * FL.N_LANES + @as(usize, delta_nbits) * 1024;

            if (dict_len > 0 and dict_encoded_size < delta_encoded_size and dict_encoded_size < normal_size) {
                // dict encode and bitpack values

                const dicted = dict_encode(dict, input);

                const DictedFL = FastLanes(u16);

                // write the dicted values right after header
                inline for (0..DictedFL.N_BITS) |nb| {
                    if (nb == dict_idx_nbits) {
                        out[0] = @bitCast(Header {
                            .encoding = .dictbitpacked,
                            .num_bits = nb,
                        });
                        const packed_data = DictedFL.Packer(nb).pack(dicted);
                        const packed_bytes_t = [packed_data.len * @sizeOf(u16)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[1 .. 1 + packed_bytes.len].* = packed_bytes;
                        break;
                    }
                } else unreachable;

                const dict_out_offset = 1 + @as(usize, dict_idx_nbits) * 1024 / 8;
                const dict_len_o: u16 = @intCast(dict_len);
                @as([][2]u8, @ptrCast(out[dict_out_offset..dict_out_offset+2]))[0] = @as([2]u8, @bitCast(dict_len_o));

                // write the dictionary last
                const dict_bytes: []const u8 = @ptrCast(dict[0..dict_len]);
                @memcpy(out[dict_out_offset+2..dict_out_offset+2+dict_bytes.len], dict_bytes);

                return dict_out_offset + dict_bytes.len + 2;
            } else if (delta_encoded_size < dict_encoded_size and delta_encoded_size < normal_size) {
                // bitpack delta encoded values
                const bases_bytes_t = [bases.len * @sizeOf(T)]u8;
                const bases_bytes: bases_bytes_t = @bitCast(bases);
                out[1 .. 1 + bases_bytes.len].* = bases_bytes;

                inline for (0..FL.N_BITS) |nb| {
                    if (nb == delta_nbits) {
                        out[0] = @bitCast(Header{
                            .num_bits = nb,
                            .encoding = .deltabitpacked,
                        });

                        const packed_data = FL.Packer(nb).pack(delta);
                        const packed_bytes_t = [packed_data.len * @sizeOf(T)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[1 + bases_bytes.len .. 1 + bases_bytes.len + packed_bytes.len].* = packed_bytes;
                        return 1 + bases_bytes.len + packed_bytes.len;
                    }
                }

                unreachable;
            } else if (normal_nbits < FL.N_BITS) {
                // normal bitpacking
                inline for (0..FL.N_BITS) |nb| {
                    if (nb == normal_nbits) {
                        out[0] = @bitCast(Header{
                            .num_bits = nb,
                            .encoding = .bitpacked,
                        });
                        const packed_data = FL.Packer(nb).pack(input);
                        const packed_bytes_t = [packed_data.len * @sizeOf(T)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[1 .. 1 + packed_bytes.len].* = packed_bytes;
                        return 1 + packed_bytes.len;
                    }
                }

                unreachable;
            } else {
                // can't use compression, just copy the bytes over 

                out[0] = @bitCast(Header {
                    .encoding = .noencoding,
                    // not used when encoding is `noencoding`
                    // putting FL.N_BITS here would overflow since FL.N_BITS can be 64
                    // and this field is a u6
                    .num_bits = 0,
                });

                out[1..1 + @sizeOf(T) * 1024].* = @bitCast(input);

                return 1 + @sizeOf(T) * 1024;
            }
        }

        fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
            if (input.len == 0) {
                return Error.Invalid;
            }
            const head: Header = @bitCast(input[0]);

            switch (head.encoding) {
                .dictbitpacked => {
                    const DictedFL = FastLanes(u16);
                    // first read the dict indices
                    
                    const dict_indices_len = @as(usize, head.num_bits) * 1024 / 8;
                    // +2 is for dict_len encoded as u16
                    if (input.len < 1 + dict_indices_len + 2) {
                        return Error.Invalid;
                    }

                    const dict_indices = inline for (0..DictedFL.N_BITS + 1) |nb| {
                        if (nb == head.num_bits) {
                            const Packer = DictedFL.Packer(nb);
                            const packed_b_len = Packer.PACKED_LEN * @sizeOf(u16);
                            const packed_bytes: [packed_b_len]u8 = input[1 .. 1 + packed_b_len].*;
                            const packed_data: [Packer.PACKED_LEN]u16 = @bitCast(packed_bytes);

                            break Packer.unpack(packed_data);
                        }
                    } else unreachable;
                    
                    const dict_len: usize = @as(u16, @bitCast(
                        @as([]const [2]u8, @ptrCast(input[1+dict_indices_len..1+dict_indices_len+2]))[0]
                    ));

                    std.debug.assert(dict_len > 0);

                    if (input.len < 1 + dict_indices_len + 2 + dict_len * @sizeOf(T)) {
                        return Error.Invalid;
                    }

                    const dict_offset = 1 + dict_indices_len + 2;
                    const dict: []const align(1) T = @ptrCast(input[dict_offset..dict_offset+dict_len*@sizeOf(T)]);

                    const data = dict_decode(dict, dict_indices);
                    
                    return .{ data, dict_offset + dict_len * @sizeOf(T) };
                },
                .deltabitpacked => {
                    const expected_input_len = 1 + @sizeOf(T) * FL.N_LANES + @as(usize, head.num_bits) * 1024 / 8;
                    if (input.len < expected_input_len) {
                        return Error.Invalid;
                    }

                    const bases_bytes: [@sizeOf(T) * FL.N_LANES]u8 = input[1 .. 1 + @sizeOf(T) * FL.N_LANES].*;
                    const bases: [FL.N_LANES]T = @bitCast(bases_bytes);
                    inline for (0..FL.N_BITS + 1) |nb| {
                        if (nb == head.num_bits) {
                            const Packer = FL.Packer(nb);
                            const packed_b_len = Packer.PACKED_LEN * @sizeOf(T);
                            const packed_bytes: [packed_b_len]u8 = input[1 + bases_bytes.len .. 1 + bases_bytes.len + packed_b_len].*;
                            const packed_data: [Packer.PACKED_LEN]T = @bitCast(packed_bytes);

                            const delta = Packer.unpack(packed_data);
                            const transposed = FL.undelta(delta, bases);
                            const data = FL.untranspose(transposed);

                            return .{ data, expected_input_len };
                        }
                    }

                    unreachable;
                },
                .bitpacked => {
                    const expected_input_len = 1 + @as(usize, head.num_bits) * 1024 / 8;
                    if (input.len < expected_input_len) {
                        return Error.Invalid;
                    }

                    inline for (0..FL.N_BITS + 1) |nb| {
                        if (nb == head.num_bits) {
                            const Packer = FL.Packer(nb);
                            const packed_b_len = Packer.PACKED_LEN * @sizeOf(T);
                            const packed_bytes: [packed_b_len]u8 = input[1 .. 1 + packed_b_len].*;
                            const packed_data: [Packer.PACKED_LEN]T = @bitCast(packed_bytes);

                            const data = Packer.unpack(packed_data);

                            return .{ data, expected_input_len };
                        }
                    }
                    
                    unreachable;
                },
                .noencoding => {
                    if (input.len < 1 + 1024 * @sizeOf(T)) {
                        return Error.Invalid;
                    }

                    const data: [1024]T = @bitCast(input[1..1+@sizeOf(T)*1024].*);

                    return .{ data, 1 + @sizeOf(T) * 1024 };
                },
            } 
        }

        /// Calculate the maximum number of output bytes needed to compress `len` integers.
        pub fn compress_bound(len: usize) usize {
            const n_blocks = (len + 1024 - 1) / 1024;
            return n_blocks * PACK_BOUND;
        }

        /// Compress the given integers into the output buffer.
        ///
        /// Output buffer size needs to be at least `compress_bound(input.len)`
        pub fn compress(noalias input: []const T, noalias out: []u8) usize {
            std.debug.assert(compress_bound(input.len) <= out.len);

            var offset: usize = 0;

            const n_whole_blocks = input.len / 1024;
            const whole_blocks: []const [1024]T = @ptrCast(input[0 .. n_whole_blocks * 1024]);
            for (whole_blocks) |block| {
                offset += pack(block, out[offset..]);
            }

            const remaining = input[n_whole_blocks * 1024 ..];
            if (remaining.len > 0) {
                var final_block: [1024]T = undefined;
                var idx: usize = 0;
                for (remaining) |v| {
                    final_block[idx] = v;
                    idx += 1;
                }

                // repeat the last element
                const last_elem = final_block[idx - 1];
                for (idx..1024) |i| {
                    final_block[i] = last_elem;
                }

                offset += pack(final_block, out[offset..]);
            }

            return offset;
        }

        /// Decompress `out.len` integers from the `input` into the `out` buffer.
        ///
        /// This function will error if the input is invalid or doesn't deserialize into exactly
        /// `input.len` integers.
        pub fn decompress(noalias input: []const u8, noalias out: []T) Error!void {
            const n_whole_blocks = out.len / 1024;
            const whole_blocks: [][1024]T = @ptrCast(out[0 .. n_whole_blocks * 1024]);

            var offset: usize = 0;

            for (whole_blocks) |*b| {
                b.*, const n_read = try unpack(input[offset..]);
                offset += n_read;
            }

            const remaining = out[n_whole_blocks * 1024 ..];
            if (remaining.len > 0) {
                const final_block, const n_read = try unpack(input[offset..]);
                offset += n_read;
                @memcpy(remaining, final_block[0..remaining.len]);
            }

            if (offset != input.len) {
                return Error.Invalid;
            }
        }
    };
}

test "zint pack" {
    var data: [1024]u32 = undefined;
    for (0..1024) |i| {
        data[i] = @as(u32, @intCast(i)) + 1023;
    }

    const Z = Zint(u32);

    var packed_d: [1 + 1024 * @sizeOf(u32)]u8 = undefined;
    const packed_size = Z.pack(data, &packed_d);

    const unpacked_d, const consumed = try Z.unpack(&packed_d);

    std.debug.assert(consumed == packed_size);

    for (0..1024) |i| {
        std.debug.assert(unpacked_d[i] == data[i]);
    }
}

test "zint compress" {
    const T = u64;
    const Z = Zint(T);

    const len = 512312;
    const input = try std.heap.page_allocator.alloc(T, len);
    var prng = std.Random.DefaultPrng.init(69);
    const rand = prng.random();
    input[0] = rand.int(T);
    for (1..len) |i| {
        input[i] = input[i - 1] +% rand.int(u8);
    }

    const compress_buf = try std.heap.page_allocator.alloc(u8, Z.compress_bound(len));

    const compressed_size = Z.compress(input, compress_buf);

    const compressed = compress_buf[0..compressed_size];

    const output = try std.heap.page_allocator.alloc(T, len);
    try Z.decompress(compressed, output);

    std.debug.assert(std.mem.eql(T, output, input));
}
