fuzz:
	zig build fuzz --fuzz -Doptimize=ReleaseSafe -j128
fuzz_debug:
	zig build fuzz --fuzz -Doptimize=Debug -j128
