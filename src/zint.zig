const std = @import("std");

const FastLanes = @import("fastlanes.zig").FastLanes;

pub const Error = error{Invalid};

pub fn Zint(comptime T: type) type {
    const FL = FastLanes(T);

    const Header = packed struct(u8) {
        is_delta: u1,
        num_bits: u7,
    };

    const DELTA_BASES_N_BITS = FL.N_BITS * FL.N_LANES;

    return struct {
        fn needed_num_bits(data: [1024]T) u7 {
            var m = @abs(data[0]);
            for (0..1024) |idx| {
                m = @max(data[idx], @abs(m));
            }
            return std.math.log2_int_ceil(T, m);
        }

        pub fn pack(input: [1024]T, noalias out: []u8) usize {
            std.debug.assert(out.len >= 1 + @sizeOf(T) * 1024);

            const normal_nbits = needed_num_bits(input);

            const transposed = FL.transpose(input);
            const bases: [FL.N_LANES]T = transposed[0..FL.N_LANES].*;
            const delta = FL.delta(input, bases);

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
                        return 1 + packed_bytes.len;
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

        // pub fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
        //     if (input.len == 0) {
        //         return Error.Invalid;
        //     }
        //     const num_bits = input[0];

        //     inline for (0..FL.N_BITS) |nb| {
        //         if (nb == num_bits) {
        //             const input_end = 1 + FL.packed_len(nb) * @sizeOf(T);
        //             if (input.len < input_end) {
        //                 return Error.Invalid;
        //             }
        //             const packed_input: [FL.packed_len(nb)]T = @bitCast(input[1..input_end]);
        //             return .{ FL.unpack(nb, packed_input), input_end };
        //         }
        //     }

        //     std.debug.assert(num_bits == FL.N_BITS);

        //     const input_end = 1 + 1024 * @sizeOf(T);
        //     if (input.len < input_end) {
        //         return Error.Invalid;
        //     }
        //     return .{ @bitCast(input[1..input_end]), input_end };
        // }
    };
}

test "zint pack" {
    var data: [1024]u64 = undefined;
    for (0..1024) |i| {
        data[i] = i;
    }

    const Z = Zint(u64);

    var packed_d: [1 + 1024 * @sizeOf(u64)]u8 = undefined;
    _ = Z.pack(data, &packed_d);

    // const unpacked_d, const consumed = Z.unpack(&packed_d);

    // std.debug.assert(consumed == packed_len);

    // for (0..1024) |i| {
    //     std.debug.assert(unpacked_d[i] == data[i]);
    // }
}
