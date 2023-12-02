MPICC=cc
INCS=-I./inc
LIBS=
LINKS=
PROGS=gen_svd dist_svd dist_pca
FILES=src/svd_algs.c src/svd_utils.c src/mmio.c src/mmio_dense.c src/kiss.c src/utils.c
CFLAGS=-Wall

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

gen_svd: src/gen_svd.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

dist_svd: src/dist_svd.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

dist_pca: src/dist_pca.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

clear:
	rm -rf svd_matrix_cases
	git clean -i

clean:
	rm -rf $(PROGS) *.dSYM *.o *.out *.mtx *.diag

