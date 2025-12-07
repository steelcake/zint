pub fn ScalarBitpack(comptime T: type) type {
    switch (T) {
        u8, u16, u32, u64 => {},
        else => @compileError("unsupported type"),
    }

    return struct {
        const N_BITS = @sizeOf(T) * 8;

        pub const Error = error{
            InvalidInput,
        };

        /// Bit packs and writes integers from the input into the output.
        ///
        /// Returns the number of elements written to the output.
        pub fn bitpack(input: []const T, reference: T, output: []align(1) T, bit_width: u8) Error!usize {
            if (bit_width > N_BITS) {
                return Error.InvalidInput;
            }

            if (bit_width == N_BITS) {
                if (input.len > output.len) {
                    return Error.InvalidInput;
                }

                if (reference == 0) {
                    @memcpy(output[0..input.len], input);
                } else {
                    for (0..input.len) |i| {
                        output[i] = input[i] -% reference;
                    }
                }

                return input.len;
            }

            if (bit_width == 0 or input.len == 0) {
                return 0;
            }

            const total_bits: u64 = @as(u64, input.len) * @as(u64, bit_width);
            const total_len = (total_bits + N_BITS - 1) / N_BITS;
            if (output.len < total_len) {
                return Error.InvalidInput;
            }

            var out_idx: usize = 0;
            var buffer: T = 0;
            var bits_in_buffer: u8 = 0;

            const mask = (@as(T, 1) << @intCast(bit_width)) - 1;

            for (input) |raw_val| {
                const val = (raw_val -% reference) & mask;

                buffer |= (val << @intCast(bits_in_buffer));

                const space_left = N_BITS - bits_in_buffer;

                if (bit_width >= space_left) {
                    output[out_idx] = buffer;
                    out_idx += 1;

                    if (bit_width > space_left) {
                        buffer = val >> @intCast(space_left);
                        bits_in_buffer = bit_width - space_left;
                    } else {
                        buffer = 0;
                        bits_in_buffer = 0;
                    }
                } else {
                    bits_in_buffer += bit_width;
                }
            }

            if (bits_in_buffer > 0) {
                output[out_idx] = buffer;
                out_idx += 1;
            }

            if (out_idx != total_len) unreachable;

            return total_len;
        }

        /// Reads packed integers from input into the output.
        ///
        /// Returns the number of elements consumed from the input.
        pub fn bitunpack(input: []align(1) const T, reference: T, output: []T, bit_width: u8) Error!usize {
            if (bit_width > N_BITS) return Error.InvalidInput;

            if (bit_width == N_BITS) {
                if (input.len < output.len) return Error.InvalidInput;

                if (reference == 0) {
                    @memcpy(output, input[0..output.len]);
                } else {
                    for (0..output.len) |i| {
                        output[i] = input[i] +% reference;
                    }
                }

                return output.len;
            }

            if (output.len == 0) {
                return 0;
            }

            if (bit_width == 0) {
                @memset(output, reference);
                return 0;
            }

            const total_bits_needed = @as(u64, output.len) * @as(u64, bit_width);
            const input_len_needed = (total_bits_needed + N_BITS - 1) / N_BITS;
            if (input.len < input_len_needed) {
                return Error.InvalidInput;
            }

            var in_idx: usize = 0;
            var buffer: T = 0;
            var bits_in_buffer: u8 = 0;

            const mask = (@as(T, 1) << @intCast(bit_width)) - 1;

            if (in_idx < input.len) {
                buffer = input[in_idx];
                in_idx += 1;
                bits_in_buffer = N_BITS;
            }

            for (output) |*out_val| {
                if (bits_in_buffer >= bit_width) {
                    out_val.* = (buffer & mask) +% reference;

                    buffer = buffer >> @intCast(bit_width);
                    bits_in_buffer -= bit_width;
                } else {
                    const low_bits = buffer;
                    const bits_taken = bits_in_buffer;
                    const bits_needed = bit_width - bits_taken;

                    if (in_idx >= input.len) unreachable;
                    const next_word = input[in_idx];
                    in_idx += 1;

                    const high_mask = (@as(T, 1) << @intCast(bits_needed)) - 1;
                    const high_bits = next_word & high_mask;

                    out_val.* = (low_bits | (high_bits << @intCast(bits_taken))) +% reference;

                    buffer = next_word >> @intCast(bits_needed);
                    bits_in_buffer = N_BITS - bits_needed;
                }
            }

            if (in_idx != input_len_needed) unreachable;

            return input_len_needed;
        }
    };
}

fn TestScalarBitpack(comptime T: type) type {
    return struct {
        const std = @import("std");
        const SBP = ScalarBitpack(T);

        const Context = struct {
            in: []T,
            packed_out: []T,
            roundtrip: []T,
        };

        fn read_input(in: []T, input_raw: []const u8) ?[]const T {
            if (input_raw.len < 1) return null;

            const has_zeroes = input_raw[0] % 2 == 0;

            const input_raw_len = (input_raw.len - 1) / @sizeOf(T) * @sizeOf(T) + 1;
            const input: []align(1) const T = @ptrCast(input_raw[1..input_raw_len]);

            const len = @min(in.len, input.len);
            @memcpy(in[0..len], input[0..len]);

            if (!has_zeroes and len > 0) {
                const replacement = in[0] +| 1;
                for (0..len) |i| {
                    if (in[i] == 0) {
                        in[i] = replacement;
                    }
                }
            }

            return in[0..len];
        }

        fn fuzz_scalar_for(ctx: Context, input: []const u8) anyerror!void {
            const in = read_input(ctx.in, input) orelse return;

            const min, const max = if (in.len > 0) std.mem.minMax(T, in) else .{ 0, 0 };

            const width = (@sizeOf(T) * 8) - @clz(max - min);
            std.debug.assert(width <= @bitSizeOf(T));

            const packed_len = SBP.bitpack(in, min, ctx.packed_out, width) catch unreachable;

            const consumed = SBP.bitunpack(ctx.packed_out, min, ctx.roundtrip[0..in.len], width) catch unreachable;

            std.debug.assert(consumed == packed_len);

            try std.testing.expectEqualSlices(T, in, ctx.roundtrip[0..in.len]);
        }

        test "fuzz scalar FOR bitpacking" {
            const max_len = 12333;
            const ctx = Context{
                .roundtrip = std.heap.page_allocator.alloc(T, max_len) catch unreachable,
                .in = std.heap.page_allocator.alloc(T, max_len) catch unreachable,
                .packed_out = std.heap.page_allocator.alloc(T, max_len) catch unreachable,
            };
            try std.testing.fuzz(ctx, fuzz_scalar_for, .{});
        }
    };
}

// Instantiate tests for all supported types
test "scalar bitpack fuzz u8" {
    _ = TestScalarBitpack(u8);
}
test "scalar bitpack fuzz u16" {
    _ = TestScalarBitpack(u16);
}
test "scalar bitpack fuzz u32" {
    _ = TestScalarBitpack(u32);
}
test "scalar bitpack fuzz u64" {
    _ = TestScalarBitpack(u64);
}
