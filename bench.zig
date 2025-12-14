const std = @import("std");
const Prng = std.Random.DefaultPrng;
const Random = std.Random;
const page_allocator = std.heap.page_allocator;
const Allocator = std.mem.Allocator;
const FixedBufferAllocator = std.heap.FixedBufferAllocator;

const zint = @import("zint");
const Zint = zint.Zint;

const TYPES = .{
    // u8,
    // i8,
    // u16,
    // i16,
    u32,
    // i32,
    // u64,
    // i64,
    // u128,
    // i128,
    // u256,
    // i256,
};

const LENGTHS: []const u32 = &.{
    // 10,
    // 69,
    // 1023,
    // 1024,
    // 1025,
    // 123321,
    1 << 18,
        // (1 << 18) + 1023,
};

const WIDTHS = .{
    7,
    15,
    // 32,
    // 33,
};

const DATASETS = .{
    // Width,
    DeltaWidth,
    // FrameWidth,
};

const ALGOS = .{
    MemCopy,
    Lz4,
    Zstd,
    // ZintBitpack,
    ZintForpack,
    ZintDeltapack,
};

const BUFFER_SIZE = 1 << 34;

const N_RUNS = 500;

pub fn main() anyerror!void {
    const mem = try page_allocator.alloc(u8, BUFFER_SIZE);
    defer page_allocator.free(mem);

    var fb_alloc = FixedBufferAllocator.init(mem);
    const alloc = fb_alloc.allocator();

    var ctx = try zint.Ctx.init(alloc);
    defer ctx.deinit(alloc);

    const ALIGN = comptime std.mem.Alignment.fromByteUnits(64);

    const input_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 31);
    defer alloc.free(input_buf);
    const output_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 31);
    defer alloc.free(output_buf);
    const compressed_buf = try alloc.alignedAlloc(u8, ALIGN, 1 << 31);

    @memset(input_buf, 69);
    @memset(output_buf, 69);
    @memset(compressed_buf, 69);

    var prng = Prng.init(69);
    const rand = prng.random();

    @setEvalBranchQuota(4096);

    inline for (TYPES) |T| {
        const input_b: []T = @ptrCast(input_buf);
        const output_b: []T = @ptrCast(output_buf);

        inline for (WIDTHS) |W| {
            if (W > @bitSizeOf(T)) continue;

            inline for (DATASETS) |ds| {
                ds(T, W).fill_input(input_b);

                for (LENGTHS) |len| {
                    std.debug.print(
                        "T={s}/DATASET={s}/LENGTH={}\n",
                        .{
                            @typeName(T),
                            ds(T, W).name(),
                            len,
                        },
                    );

                    inline for (ALGOS) |algo| {
                        const alg = algo(T).init();
                        defer alg.deinit();

                        const res = try bench_one(T, alg, rand, len, input_b, compressed_buf, output_b);

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
}

fn Width(comptime T: type, comptime W: comptime_int) type {
    return struct {
        fn fill_input(input: []T) void {
            var prng = Prng.init(0);
            const rand = prng.random();
            for (0..input.len) |i| {
                input[i] = if (@bitSizeOf(T) > W) mask(T, W, rand.int(T)) else rand.int(T);
            }
        }

        fn name() []const u8 {
            return std.fmt.comptimePrint("width{}", .{W});
        }
    };
}

fn mask(comptime T: type, comptime W: comptime_int, v: T) T {
    const m: T = (1 << W) - 1;

    if (@typeInfo(T).int.signedness == .signed) {
        return (v << (@bitSizeOf(T) - W)) >> (@bitSizeOf(T) - W);
    } else {
        return v & m;
    }
}

fn DeltaWidth(comptime T: type, comptime W: comptime_int) type {
    return struct {
        fn fill_input(input: []T) void {
            var prng = Prng.init(0);
            const rand = prng.random();
            var val = rand.int(T);
            for (0..input.len) |i| {
                input[i] = val;
                val +%= if (@bitSizeOf(T) > W) @intCast(@abs(mask(T, W, rand.int(T))) / 2) else rand.int(T);
            }
        }

        fn name() []const u8 {
            return std.fmt.comptimePrint("delta_width{}", .{W});
        }
    };
}

fn FrameWidth(comptime T: type, comptime W: comptime_int) type {
    return struct {
        fn fill_input(input: []T) void {
            var prng = Prng.init(0);
            const rand = prng.random();
            var val = rand.int(T);
            for (0..input.len) |i| {
                input[i] = val;
                val +%= if (@bitSizeOf(T) > W) mask(T, W, rand.int(T)) else rand.int(T);
            }
        }

        fn name() []const u8 {
            return std.fmt.comptimePrint("frame_width{}", .{W});
        }
    };
}

const ZSTD_LEVEL = 1;

fn Zstd(comptime T: type) type {
    return struct {
        const Self = @This();

        zstd_cctx: *sys.ZSTD_CCtx,
        zstd_cstate: []align(8) u8,

        zstd_dctx: *sys.ZSTD_DCtx,
        zstd_dstate: []align(8) u8,

        fn init() Self {
            const alloc = page_allocator;

            const zstd_dctx_cap = sys.ZSTD_estimateDCtxSize();

            const zstd_dstate = alloc.alignedAlloc(
                u8,
                std.mem.Alignment.fromByteUnits(8),
                @intCast(zstd_dctx_cap),
            ) catch @panic("unreachable");

            const zstd_dctx = sys.ZSTD_initStaticDCtx(zstd_dstate.ptr, zstd_dstate.len) orelse {
                @panic("failed to init zstd decompression ctx");
            };

            const zstd_cctx_cap = sys.ZSTD_estimateCCtxSize(ZSTD_LEVEL);

            const zstd_cstate = alloc.alignedAlloc(
                u8,
                std.mem.Alignment.fromByteUnits(8),
                @intCast(zstd_cctx_cap),
            ) catch @panic("unreachable");

            const zstd_cctx = sys.ZSTD_initStaticCCtx(zstd_cstate.ptr, zstd_cstate.len) orelse {
                @panic("failed to init zstd compression ctx");
            };

            return .{
                .zstd_dctx = zstd_dctx,
                .zstd_cctx = zstd_cctx,
                .zstd_cstate = zstd_cstate,
                .zstd_dstate = zstd_dstate,
            };
        }

        fn deinit(self: Self) void {
            page_allocator.free(self.zstd_dstate);
            page_allocator.free(self.zstd_cstate);
        }

        fn compress_bound(self: Self, len: u32) u32 {
            _ = self;
            return @intCast(sys.ZSTD_compressBound(@intCast(len * @sizeOf(T))));
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            return try zstd_compress(self.zstd_cctx, @as([]const u8, @ptrCast(input)), output, ZSTD_LEVEL);
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            return try zstd_decompress(self.zstd_dctx, input, @as([]u8, @ptrCast(output)));
        }

        fn name() []const u8 {
            return "zstd(level 1)";
        }
    };
}

fn Lz4(comptime T: type) type {
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
            return @intCast(sys.LZ4_compressBound(@intCast(len * @sizeOf(T))));
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            _ = self;
            return try lz4_compress(@as([]const u8, @ptrCast(input)), output);
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            _ = self;
            return try lz4_decompress(input, @as([]u8, @ptrCast(output)));
        }

        fn name() []const u8 {
            return "lz4";
        }
    };
}

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
                .ctx = zint.Ctx.init(std.heap.page_allocator) catch @panic("unreachable"),
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

fn ZintForpack(comptime T: type) type {
    const Z = Zint(T);

    return struct {
        const Self = @This();

        ctx: zint.Ctx,

        fn init() Self {
            return .{
                .ctx = zint.Ctx.init(std.heap.page_allocator) catch @panic("unreachable"),
            };
        }

        fn deinit(self: Self) void {
            self.ctx.deinit(std.heap.page_allocator);
        }

        fn compress_bound(self: Self, len: u32) u32 {
            _ = self;
            return Z.forpack_compress_bound(len);
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            return try Z.forpack_compress(self.ctx, input, output);
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            return try Z.forpack_decompress(self.ctx, input, output);
        }

        fn name() []const u8 {
            return "zint_forpack";
        }
    };
}

fn ZintDeltapack(comptime T: type) type {
    const Z = Zint(T);

    return struct {
        const Self = @This();

        ctx: zint.Ctx,

        fn init() Self {
            return .{
                .ctx = zint.Ctx.init(std.heap.page_allocator) catch @panic("unreachable"),
            };
        }

        fn deinit(self: Self) void {
            self.ctx.deinit(std.heap.page_allocator);
        }

        fn compress_bound(self: Self, len: u32) u32 {
            _ = self;
            return Z.deltapack_compress_bound(len);
        }

        fn compress(self: Self, noalias input: []const T, noalias output: []u8) anyerror!usize {
            return try Z.deltapack_compress(self.ctx, input, output);
        }

        fn decompress(self: Self, noalias input: []const u8, noalias output: []T) anyerror!usize {
            return try Z.deltapack_decompress(self.ctx, input, output);
        }

        fn name() []const u8 {
            return "zint_deltapack";
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
    rand: Random,
    len: u32,
    input_buf: []const T,
    compressed_buf: []u8,
    output_buf: []T,
) anyerror!Result {
    const bound = algo.compress_bound(len);

    var compress_time_ns: u64 = 0;
    var compressed_len: usize = 0;

    const n_compressed = compressed_buf.len / bound;

    var compressed_indices: [N_RUNS]usize = undefined;
    var compressed_sizes: [N_RUNS]usize = undefined;
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

        const cl = try algo.compress(input, compressed);
        compressed_len += cl;

        compressed_sizes[run_idx] = cl;

        compress_time_ns += t.lap();
    }

    var decompress_time_ns: u64 = 0;

    for (0..N_RUNS) |run_idx| {
        const compressed_idx = compressed_indices[run_idx];
        const compressed = compressed_buf[bound * compressed_idx .. bound * compressed_idx + compressed_sizes[run_idx]];

        const output = rand_slice(output_buf, rand, len);

        var t = timer();

        _ = try algo.decompress(compressed, output);

        decompress_time_ns += t.lap();

        if (!std.mem.eql(T, input_slices[run_idx], output)) {
            @panic("mismatch");
        }
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
    return std.time.Timer.start() catch @panic("unreachable");
}

const sys = @cImport({
    @cDefine("ZSTD_STATIC_LINKING_ONLY", "1");
    @cInclude("zstd.h");
    @cInclude("lz4.h");
    @cInclude("lz4hc.h");
});

const CompressError = error{
    CompressFail,
};

const DecompressError = error{
    DecompressFail,
};

fn lz4_compress(src: []const u8, dst: []u8) CompressError!usize {
    const lz4_size = sys.LZ4_compress_default(
        @ptrCast(src.ptr),
        @ptrCast(dst.ptr),
        @intCast(src.len),
        @intCast(dst.len),
    );
    if (lz4_size != 0) {
        return @intCast(lz4_size);
    } else {
        return CompressError.CompressFail;
    }
}

fn zstd_compress(ctx: *sys.ZSTD_CCtx, src: []const u8, dst: []u8, level: i32) CompressError!usize {
    const res = sys.ZSTD_compressCCtx(ctx, dst.ptr, dst.len, src.ptr, src.len, level);
    if (sys.ZSTD_isError(res) == 0) {
        return res;
    } else {
        return CompressError.CompressFail;
    }
}

fn zstd_decompress(ctx: *sys.ZSTD_DCtx, src: []const u8, dst: []u8) DecompressError!usize {
    const res = sys.ZSTD_decompressDCtx(ctx, dst.ptr, dst.len, src.ptr, src.len);
    if (sys.ZSTD_isError(res) != 0) {
        return DecompressError.DecompressFail;
    }

    return @intCast(res);
}

fn lz4_decompress(src: []const u8, dst: []u8) DecompressError!usize {
    const res = sys.LZ4_decompress_safe(
        @ptrCast(src.ptr),
        @ptrCast(dst.ptr),
        @intCast(src.len),
        @intCast(dst.len),
    );
    if (res < 0) {
        return DecompressError.DecompressFail;
    }

    return @intCast(res);
}
