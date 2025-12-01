const std = @import("std");

const FastLanes = @import("fastlanes.zig").FastLanes;

pub const Error = error{Invalid};

const Header = packed struct(u16) {
    is_delta: u1,
    num_bits: u7,
    is_zigzag: u1,
    _padding: u7 = 0,
};

const HeaderSize = @sizeOf(Header);

pub fn Zint(comptime T: type) type {
    const is_signed = switch (T) {
        u8, u16, u32, u64 => false,
        i8, i16, i32, i64 => true,
        else => @compileError("Zint only supports u8, u16, u32 and u64"),
    };

    const U = std.meta.Int(.unsigned, @typeInfo(T).int.bits);

    const FL = FastLanes(U);

    const DELTA_BASES_N_BITS = FL.N_BITS * FL.N_LANES;

    return struct {
        fn zigzag_encode(data: [1024]T) struct { [1024]U, bool } {
            if (!is_signed) {
                return .{ data, false };
            }

            for (0..1024) |i| {
                if (data[i] < 0) {
                    break;
                }
                return .{ @bitCast(data), false };
            }

            var out: [1024]U = undefined;
            for (0..1024) |i| {
                const v = data[i];
                out[i] = @bitCast(std.math.shr(T, v, (FL.N_BITS - 1)) ^ std.math.shl(T, v, 1));
            }

            return .{ out, true };
        }

        pub fn negate(x: U) U {
            return ~x +% 1;
        }

        fn zigzag_decode(data: [1024]U) [1024]T {
            if (!is_signed) {
                return data;
            }

            var out: [1024]U = undefined;
            for (0..1024) |i| {
                const v = data[i];
                out[i] = std.math.shr(U, v, 1) ^ negate(v & 1);
            }

            return @bitCast(out);
        }

        fn needed_num_bits(data: [1024]U) u7 {
            var m = data[0];
            for (0..1024) |idx| {
                m = @max(data[idx], m);
            }

            if (m == 0) {
                return 0;
            }

            return @as(u7, std.math.log2_int(U, m)) + 1;
        }

        const PACK_BOUND = HeaderSize + @sizeOf(U) * 1024;

        fn pack(input_v: [1024]T, noalias out: []u8) usize {
            const input, const is_zigzag = zigzag_encode(input_v);

            std.debug.assert(out.len >= PACK_BOUND);

            const normal_nbits = needed_num_bits(input);

            const transposed = FL.transpose(input);
            const bases: [FL.N_LANES]U = transposed[0..FL.N_LANES].*;
            const delta = FL.delta(transposed, bases);

            const delta_nbits = needed_num_bits(delta);

            const delta_is_better = @as(usize, delta_nbits) * 1024 + DELTA_BASES_N_BITS < @as(usize, normal_nbits) * 1024;

            if (delta_is_better) {
                const bases_bytes_t = [bases.len * @sizeOf(U)]u8;
                const bases_bytes: bases_bytes_t = @bitCast(bases);
                out[HeaderSize .. HeaderSize + bases_bytes.len].* = bases_bytes;

                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == delta_nbits) {
                        out[0..2].* = @bitCast(Header{
                            .num_bits = nb,
                            .is_delta = 1,
                            .is_zigzag = @intFromBool(is_zigzag),
                        });
                        out[1] = @intFromBool(is_zigzag);
                        const packed_data = FL.Packer(nb).pack(delta);
                        const packed_bytes_t = [packed_data.len * @sizeOf(U)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[HeaderSize + bases_bytes.len .. HeaderSize + bases_bytes.len + packed_bytes.len].* = packed_bytes;
                        return HeaderSize + bases_bytes.len + packed_bytes.len;
                    }
                }

                unreachable;
            } else {
                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == normal_nbits) {
                        out[0..2].* = @bitCast(Header{
                            .num_bits = nb,
                            .is_delta = 0,
                            .is_zigzag = @intFromBool(is_zigzag),
                        });
                        const packed_data = FL.Packer(nb).pack(input);
                        const packed_bytes_t = [packed_data.len * @sizeOf(U)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[HeaderSize .. HeaderSize + packed_bytes.len].* = packed_bytes;
                        return HeaderSize + packed_bytes.len;
                    }
                }

                unreachable;
            }
        }

        fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
            if (input.len < HeaderSize) {
                return Error.Invalid;
            }
            const head: Header = @bitCast(input[0..2].*);

            const is_delta: bool = @bitCast(head.is_delta);

            if (is_delta) {
                const expected_input_len = HeaderSize + @sizeOf(U) * FL.N_LANES + @as(usize, head.num_bits) * 1024 / 8;
                if (input.len < expected_input_len) {
                    return Error.Invalid;
                }

                const bases_bytes: [@sizeOf(U) * FL.N_LANES]u8 = input[HeaderSize .. HeaderSize + @sizeOf(U) * FL.N_LANES].*;
                const bases: [FL.N_LANES]U = @bitCast(bases_bytes);
                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == head.num_bits) {
                        const Packer = FL.Packer(nb);
                        const packed_b_len = Packer.PACKED_LEN * @sizeOf(U);
                        const packed_bytes: [packed_b_len]u8 = input[HeaderSize + bases_bytes.len .. HeaderSize + bases_bytes.len + packed_b_len].*;
                        const packed_data: [Packer.PACKED_LEN]U = @bitCast(packed_bytes);

                        const delta = Packer.unpack(packed_data);
                        const transposed = FL.undelta(delta, bases);
                        const d = FL.untranspose(transposed);
                        const data = if (head.is_zigzag == 1) zigzag_decode(d) else @as([1024]T, @bitCast(d));

                        return .{ data, expected_input_len };
                    }
                }

                unreachable;
            } else {
                const expected_input_len = HeaderSize + @as(usize, head.num_bits) * 1024 / 8;
                if (input.len < expected_input_len) {
                    return Error.Invalid;
                }

                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == head.num_bits) {
                        const Packer = FL.Packer(nb);
                        const packed_b_len = Packer.PACKED_LEN * @sizeOf(T);
                        const packed_bytes: [packed_b_len]u8 = input[HeaderSize .. HeaderSize + packed_b_len].*;
                        const packed_data: [Packer.PACKED_LEN]U = @bitCast(packed_bytes);

                        const d = Packer.unpack(packed_data);
                        const data = if (head.is_zigzag == 1) zigzag_decode(d) else @as([1024]T, @bitCast(d));

                        return .{ data, expected_input_len };
                    }
                }

                unreachable;
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

    var packed_d: [HeaderSize + 1024 * @sizeOf(u32)]u8 = undefined;
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
