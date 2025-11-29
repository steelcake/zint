const std = @import("std");

const fastlanes = @import("fastlanes.zig");
const FastLanes = fastlanes.FastLanes;

// pub const Error = error{Invalid};

// fn Zint(comptime T: type) type {
//     const FL = FastLanes(T);

//     const Header = packed struct(u8) {
//         is_delta: u1,
//         _padding: u1,
//         num_bits: u6,
//     };

//     return struct {
//         fn max(data: [1024]T) T {
//             var m = @abs(data[0]);
//             for (0..1024) |idx| {
//                 m = @max(data[idx], @abs(m));
//             }
//             return m;
//         }

//         pub fn pack(input: [1024]T, noalias out: []u8) usize {
//             std.debug.assert(out.len >= 1 + @sizeOf(T) * 1024);
//             const num_bits: u8 = std.math.log2_int_ceil(T, max);

//             out[0] = num_bits;
//             inline for (0..FL.N_BITS) |nb| {
//                 if (nb == num_bits) {
//                     const packed_input = FL.pack(nb, input);
//                     @as([][FL.packed_len(nb)]T, @ptrCast(out[1..1+nb*1024/8]))[0] = packed_input;
//                     return 1 + nb * 1024 / 8;
//                 }
//             }

//             std.debug.assert(num_bits == FL.N_BITS);
//             @memcpy(out[1..1+@sizeOf(T) * 1024], input);

//             return 1 + @sizeOf(T) * 1024;
//         }

//         pub fn unpack(noalias input: []const u8) Error!struct { [1024]T, usize } {
//             if (input.len == 0) {
//                 return Error.Invalid;
//             }
//             const num_bits = input[0];

//             inline for (0..FL.N_BITS) |nb| {
//                 if (nb == num_bits) {
//                     const input_end = 1 + FL.packed_len(nb) * @sizeOf(T);
//                     if (input.len < input_end) {
//                         return Error.Invalid;
//                     }
//                     const packed_input: [FL.packed_len(nb)]T = @bitCast(input[1..input_end]);
//                     return .{ FL.unpack(nb, packed_input), input_end };
//                 }
//             }

//             std.debug.assert(num_bits == FL.N_BITS);

//             const input_end = 1 + 1024 * @sizeOf(T);
//             if (input.len < input_end) {
//                 return Error.Invalid;
//             }
//             return .{ @bitCast(input[1..input_end]), input_end };
//         }
//     };
// }

// test "zint pack" {
//     var data: [1024]u64 = undefined;
//     for (0..1024) |i| {
//         data[i] = i;
//     }

//     const Z = Zint(u64);

//     var packed_d: [1024 * @sizeOf(u64)]u8 = undefined;
//     const packed_len = Z.pack(data, &packed_d);

//     const unpacked_d, const consumed = Z.unpack(&packed_d);

//     std.debug.assert(consumed == packed_len);

//     for (0..1024) |i| {
//         std.debug.assert(unpacked_d[i] == data[i]);
//     }
// }
//

test {
    _ = fastlanes;
}
