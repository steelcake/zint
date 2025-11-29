const std = @import("std");

const FastLanes = @import("fastlanes.zig").FastLanes;

pub const Error = error{Invalid};

const Header = packed struct(u8) {
    is_delta: u1,
    num_bits: u7,
};

pub fn Zint(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("Zint only supports u8, u16, u32 and u64"),
    }

    const FL = FastLanes(T);

    const DELTA_BASES_N_BITS = FL.N_BITS * FL.N_LANES;

    return struct {
        fn needed_num_bits(data: [1024]T) u7 {
            var m = data[0];
            for (0..1024) |idx| {
                m = @max(data[idx], m);
            }

            if (m == 0) {
                return 0;
            }

            return std.math.log2_int_ceil(T, m) + 1;
        }

        const PACK_BOUND = 1 + @sizeOf(T) * 1024;

        fn pack(input: [1024]T, noalias out: []u8) usize {
            std.debug.assert(out.len >= PACK_BOUND);

            const normal_nbits = needed_num_bits(input);

            const transposed = FL.transpose(input);
            const bases: [FL.N_LANES]T = transposed[0..FL.N_LANES].*;
            const delta = FL.delta(transposed, bases);

            const delta_nbits = needed_num_bits(delta);

            const delta_is_better = @as(usize, delta_nbits) * 1024 + DELTA_BASES_N_BITS < @as(usize, normal_nbits) * 1024;

            if (delta_is_better) {
                const bases_bytes_t = [bases.len * @sizeOf(T)]u8;
                const bases_bytes: bases_bytes_t = @bitCast(bases);
                out[1 .. 1 + bases_bytes.len].* = bases_bytes;

                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == delta_nbits) {
                        out[0] = @bitCast(Header{
                            .num_bits = nb,
                            .is_delta = 1,
                        });
                        const packed_data = FL.Packer(nb).pack(delta);
                        const packed_bytes_t = [packed_data.len * @sizeOf(T)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[1 + bases_bytes.len .. 1 + bases_bytes.len + packed_bytes.len].* = packed_bytes;
                        return 1 + bases_bytes.len + packed_bytes.len;
                    }
                }

                unreachable;
            } else {
                inline for (0..FL.N_BITS + 1) |nb| {
                    if (nb == normal_nbits) {
                        out[0] = @bitCast(Header{
                            .num_bits = nb,
                            .is_delta = 0,
                        });
                        const packed_data = FL.Packer(nb).pack(input);
                        const packed_bytes_t = [packed_data.len * @sizeOf(T)]u8;
                        const packed_bytes: packed_bytes_t = @bitCast(packed_data);
                        out[1 .. 1 + packed_bytes.len].* = packed_bytes;
                        return 1 + packed_bytes.len;
                    }
                }

                unreachable;
            }
        }

        fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
            if (input.len == 0) {
                return Error.Invalid;
            }
            const head: Header = @bitCast(input[0]);

            const is_delta: bool = @bitCast(head.is_delta);

            if (is_delta) {
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
            } else {
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
            }
        }

        /// Calculate the maximum number of output bytes needed to compress `len` integers.
        pub fn compress_bound(len: usize) usize {
            const n_blocks = (len + 1024 - 1) / 1024;
            return n_blocks * PACK_BOUND;
        }

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
