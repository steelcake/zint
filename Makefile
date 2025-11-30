fuzz:
	zig build fuzz --fuzz -Doptimize=ReleaseSafe -j64
fuzz_debug:
	zig build fuzz --fuzz -Doptimize=Debug -j64
