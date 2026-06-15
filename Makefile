# prime-octal build
#
# Note on GENCODE: this machine pairs nvcc 12.0 with an RTX 5090 (sm_120,
# Blackwell), which nvcc 12.0 cannot target directly. We therefore embed
# compute_80 PTX only and let the (much newer) driver JIT it for the GPU.
# The first run after a build pays a few seconds of JIT; the driver caches it.
NVCC      ?= nvcc
GENCODE   ?= -gencode arch=compute_80,code=compute_80
NVCCFLAGS ?= -O3 -std=c++17
BIN       := bin

HDRS := src/octal_core.h src/sieve.cuh src/post.h src/primality.h

all: $(BIN)/prime_octal $(BIN)/test_prime_octal $(BIN)/lattice

$(BIN)/prime_octal: src/main.cu $(HDRS)
	@mkdir -p $(BIN)
	$(NVCC) $(NVCCFLAGS) $(GENCODE) -o $@ $<

$(BIN)/test_prime_octal: src/test_main.cu $(HDRS)
	@mkdir -p $(BIN)
	$(NVCC) $(NVCCFLAGS) $(GENCODE) -o $@ $<

$(BIN)/lattice: src/lattice.cu src/primality.h src/octal_core.h
	@mkdir -p $(BIN)
	$(NVCC) $(NVCCFLAGS) $(GENCODE) -o $@ $<

test: $(BIN)/test_prime_octal $(BIN)/lattice
	$(BIN)/test_prime_octal
	$(BIN)/lattice --selftest

run: $(BIN)/prime_octal
	$(BIN)/prime_octal --octal-digits 10 --out results

clean:
	rm -rf $(BIN)

.PHONY: all test run clean
