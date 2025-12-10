const std = @import("std");
const Prng = std.Random.DefaultPrng;
const page_allocator = std.heap.page_allocator;
const Allocator = std.mem.Allocator;
const posix = std.posix;
const linux = std.os.linux;
const FixedBufferAllocator = std.heap.FixedBufferAllocator;

const zint = @import("zint");
const Zint = zint.Zint;
const Ctx = zint.Ctx;

const TYPES = .{
    u8, u16, u32, u64, u128, u256,
    i8, i16, i32, i64, i128, i256,
};

const LENGTHS: []const u32 = &.{
    1023, 1024, 1025, 123321, 1 << 25, (1 << 25) + 1023,
};

pub fn main() anyerror!void {
    const mem = alloc_thp(1 << 33).?;
    defer posix.munmap(mem);

    var fb_alloc = FixedBufferAllocator.init(mem);
    const alloc = fb_alloc.allocator();

    inline for (TYPES) |T| {
        for (LENGTHS) |len| {
            fb_alloc.reset();

            var ctx = try Ctx.init(alloc);
            defer ctx.deinit(alloc);

            const input = try alloc.alloc(T, len);
            defer alloc.free(input);

            const compressed = try alloc.alloc(u8, Zint(T).bitpack_compress_bound(len));
            defer alloc.free(compressed);

            const output = try alloc.alloc(T, len);
            defer alloc.free(output);

            var prng = Prng.init(0);
            const rand = prng.random();
            for (0..len) |i| {
                input[i] = rand.int(u7);
            }

            // _ = try bench_one(T, ctx, input, compressed, output);
            const res = try bench_one(T, ctx, input, compressed, output);

            std.debug.print("T={s}\tlen={}\tcompress={d:.3}GB/s\tdecompress={d:.3}GB/s\tratio={d}\n", .{
                @typeName(T),
                len,
                res.compress_gb_s,
                res.decompress_gb_s,
                res.ratio,
            });
        }
    }
}

const Result = struct {
    compress_gb_s: f64,
    decompress_gb_s: f64,
    ratio: usize,
};

fn bench_one(
    comptime T: type,
    ctx: Ctx,
    input: []const T,
    compressed: []u8,
    output: []T,
) anyerror!Result {
    const Z = Zint(T);

    var t = Timer.start();

    const compressed_size = try Z.bitpack_compress(ctx, input, compressed);

    const size = input.len * @sizeOf(T);

    const compress_gb_s = t.gb_s(size);

    const cs = try Z.bitpack_decompress(ctx, compressed[0..compressed_size], output);
    std.debug.assert(cs == compressed_size);

    const decompress_gb_s = t.gb_s(size);

    std.debug.assert(std.mem.eql(T, input, output));

    return .{
        .compress_gb_s = compress_gb_s,
        .decompress_gb_s = decompress_gb_s,
        .ratio = compressed_size * 100 / size,
    };
}

const Timer = struct {
    inner: std.time.Timer,

    fn start() Timer {
        return .{
            .inner = std.time.Timer.start() catch unreachable,
        };
    }

    fn gb_s(self: *Timer, size: usize) f64 {
        const s: f64 = @floatFromInt(size);
        const t: f64 = @floatFromInt(self.inner.lap());
        return s / t;
    }
};

fn alloc_thp(size: usize) ?[]align(1 << 12) u8 {
    if (size == 0) {
        return null;
    }
    const alloc_size = align_forward(size, 1 << 21);
    const page = mmap_wrapper(alloc_size, 0) orelse return null;
    posix.madvise(page.ptr, page.len, posix.MADV.HUGEPAGE) catch {
        posix.munmap(page);
        return null;
    };
    return page;
}

fn mmap_wrapper(size: usize, huge_page_flag: u32) ?[]align(1 << 12) u8 {
    if (size == 0) {
        return null;
    }
    const flags = linux.MAP{ .TYPE = .PRIVATE, .ANONYMOUS = true, .HUGETLB = huge_page_flag != 0, .POPULATE = true, .LOCKED = true };
    const flags_int: u32 = @bitCast(flags);
    const flags_f: linux.MAP = @bitCast(flags_int | huge_page_flag);
    const page = posix.mmap(null, size, posix.PROT.READ | posix.PROT.WRITE, flags_f, -1, 0) catch return null;
    return page;
}

pub fn align_forward(addr: usize, alignment: usize) usize {
    return (addr +% alignment -% 1) & ~(alignment -% 1);
}
