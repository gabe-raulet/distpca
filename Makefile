CC=clang
MPICC=mpicc
INCS=-I./inc -I/opt/homebrew/Cellar/openblas/0.3.24/include -I/opt/homebrew/include -I./inc
LIBS=-L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/openblas/0.3.24/lib
LINKS=-llapacke -lopenblas -lstdc++
PROGS=gen_svd dist_svd dist_pca
FILES=src/svd_algs.cpp src/svd_utils.cpp src/mmio.cpp src/mmio_dense.cpp src/kiss.cpp src/utils.cpp
CFLAGS=-Wall -std=c++17

D?=0

ifeq ($(D), 1)
CFLAGS+=-O0 -g -fsanitize=address -fno-omit-frame-pointer
else
CFLAGS+=-O2
endif

all: $(PROGS)

test: dist_svd gen_svd
	python gen_svd_cases.py
	python test_svd_cases.py

gen_svd: src/gen_svd.cpp $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

dist_svd: src/dist_svd.cpp $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

dist_pca: src/dist_pca.cpp $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

clear:
	rm -rf svd_matrix_cases
	git clean -i

clean:
	rm -rf $(PROGS) *.dSYM *.o *.out *.mtx *.diag

