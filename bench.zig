const std = @import("std");
const Prng = std.Random.DefaultPrng;
const Random = std.Random;
const page_allocator = std.heap.page_allocator;
const Allocator = std.mem.Allocator;
const FixedBufferAllocator = std.heap.FixedBufferAllocator;

const zint = @import("zint");
const Zint = zint.Zint;

const TYPES = .{
    u8, i8, u16, i16, u32, i32, u64, i64, u128, i128, u256, i256,
};

const LENGTHS: []const u32 = &.{
    1023, 1024, 1025, 123321, 1 << 18, (1 << 18) + 1023,
};

const BUFFER_SIZE = 1 << 34;

const N_RUNS = 100;

pub fn main() anyerror!void {
    const mem = try page_allocator.alloc(u8, BUFFER_SIZE);
    defer page_allocator.free(mem);

    var fb_alloc = FixedBufferAllocator.init(mem);
    const alloc = fb_alloc.allocator();

    var ctx = try zint.Ctx.init(alloc);
    defer ctx.deinit(alloc);

    const ALIGN = comptime std.mem.Alignment.fromByteUnits(64);

    const input_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 30);
    defer alloc.free(input_buf);
    const output_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 30);
    defer alloc.free(output_buf);
    const compressed_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 31);

    @memset(input_buf, 69);
    @memset(output_buf, 69);
    @memset(compressed_buf, 69);

    inline for (TYPES) |T| {
        const input_b: []T = @ptrCast(input_buf);
        const output_b: []T = @ptrCast(output_buf);

        inline for (DATASETS) |ds| {
            ds(T).fill_input(input_b);

            for (LENGTHS) |len| {
                std.debug.print(
                    "T={s}/DATASET={s}/LENGTH={}\n",
                    .{
                        @typeName(T),
                        ds(T).name(),
                        len,
                    },
                );

                inline for (ALGOS) |algo| {
                    const alg = algo(T).init();
                    defer alg.deinit();

                    const res = try bench_one(T, alg, len, input_b, compressed_buf, output_b);

                    std.debug.print(
                        "ALGO={s}\tcompress={d:.3}GB/s\tdecompress={d:.3}GB/s\tratio={d}\n",
                        .{
                            algo(T).name(),
                            res.compress_gb_s,
                            res.decompress_gb_s,
                            res.ratio,
                        },
                    );
                }

                std.debug.print("###########################\n", .{});
            }
        }
    }
}

const DATASETS = .{Width7};

fn Width7(comptime T: type) type {
    return struct {
        fn fill_input(input: []T) void {
            var prng = Prng.init(0);
            const rand = prng.random();
            for (0..input.len) |i| {
                input[i] = rand.int(u7);
            }
        }

        fn name() []const u8 {
            return "width7";
        }
    };
}

const ALGOS = .{ ZintBitpack, MemCopy };

fn MemCopy(comptime T: type) type {
    return struct {
        const Self = @This();

        fn init() Self {
            return .{};
        }

        fn deinit(self: Self) void {
            _ = self;
        }

        fn compress_bound(self: Self, len: u32) u32 {
            _ = self;
            return len * @sizeOf(T);
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            _ = self;
            @memcpy(output[0 .. @sizeOf(T) * input.len], @as([]const u8, @ptrCast(input)));
            return @sizeOf(T) * input.len;
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            _ = self;
            @memcpy(@as([]u8, @ptrCast(output)), input[0 .. @sizeOf(T) * output.len]);
            return @sizeOf(T) * output.len;
        }

        fn name() []const u8 {
            return "memcopy";
        }
    };
}

fn ZintBitpack(comptime T: type) type {
    const Z = Zint(T);

    return struct {
        const Self = @This();

        ctx: zint.Ctx,

        fn init() Self {
            return .{
                .ctx = zint.Ctx.init(std.heap.page_allocator) catch unreachable,
            };
        }

        fn deinit(self: Self) void {
            self.ctx.deinit(std.heap.page_allocator);
        }

        fn compress_bound(self: Self, len: u32) u32 {
            _ = self;
            return Z.bitpack_compress_bound(len);
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            return try Z.bitpack_compress(self.ctx, input, output);
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            return try Z.bitpack_decompress(self.ctx, input, output);
        }

        fn name() []const u8 {
            return "zint_bitpack";
        }
    };
}

const Result = struct {
    compress_gb_s: f64,
    decompress_gb_s: f64,
    ratio: usize,
};

fn rand_slice(src: anytype, rand: Random, len: usize) @TypeOf(src) {
    const offset = rand.int(usize) % (src.len + 1 - len);
    return src[offset .. offset + len];
}

fn bench_one(
    comptime T: type,
    algo: anytype,
    len: u32,
    input_buf: []const T,
    compressed_buf: []u8,
    output_buf: []T,
) anyerror!Result {
    const bound = algo.compress_bound(len);

    var prng = Prng.init(69);
    const rand = prng.random();

    var compress_time_ns: u64 = 0;
    var compressed_len: usize = 0;

    const n_compressed = compressed_buf.len / bound;

    var compressed_indices: [N_RUNS]usize = undefined;
    var input_slices: [N_RUNS][]const T = undefined;

    for (0..N_RUNS) |run_idx| {
        const input = rand_slice(input_buf, rand, len);
        input_slices[run_idx] = input;

        const compressed_idx = while (true) {
            const idx = rand.int(usize) % n_compressed;

            const found = for (0..run_idx) |i| {
                if (compressed_indices[i] == idx) {
                    break true;
                }
            } else false;

            if (!found) {
                compressed_indices[run_idx] = idx;
                break idx;
            }
        };
        const compressed = compressed_buf[bound * compressed_idx .. bound * (compressed_idx + 1)];

        var t = timer();

        compressed_len += try algo.compress(input, compressed);

        compress_time_ns += t.lap();
    }

    var decompress_time_ns: u64 = 0;

    for (0..N_RUNS) |run_idx| {
        const compressed_idx = compressed_indices[run_idx];
        const compressed = compressed_buf[bound * compressed_idx .. bound * (compressed_idx + 1)];

        const output = rand_slice(output_buf, rand, len);

        var t = timer();

        _ = try algo.decompress(compressed, output);

        decompress_time_ns += t.lap();

        std.debug.assert(std.mem.eql(T, input_slices[run_idx], output));
    }

    const total_size = @as(usize, len) * @sizeOf(T) * N_RUNS;

    const ratio = compressed_len * 100 / total_size;

    const compress_gb_s = @as(f64, @floatFromInt(total_size)) / @as(f64, @floatFromInt(compress_time_ns));
    const decompress_gb_s = @as(f64, @floatFromInt(total_size)) / @as(f64, @floatFromInt(decompress_time_ns));

    return .{
        .compress_gb_s = compress_gb_s,
        .decompress_gb_s = decompress_gb_s,
        .ratio = ratio,
    };
}

fn timer() std.time.Timer {
    return std.time.Timer.start() catch unreachable;
}
