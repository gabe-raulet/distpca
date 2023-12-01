CC=clang
MPICC=mpicc
INCS=-I/opt/homebrew/Cellar/openblas/0.3.24/include -I/opt/homebrew/include -I./inc
LIBS=-L/opt/homebrew/Cellar/lapack/3.12.0/lib -L/opt/homebrew/Cellar/openblas/0.3.24/lib
LINKS=-llapacke -lopenblas
PROGS=gen_svd dist_svd proto_pca
FILES=svd_algs.c svd_utils.c mmio.c mmio_dense.c kiss.c utils.c
CFLAGS=-Wall

D?=0

ifeq ($(D), 1)
CFLAGS+=-O0 -g -fsanitize=address -fno-omit-frame-pointer
else
CFLAGS+=-O2
endif

all: $(PROGS)

gen_svd: gen_svd.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

proto_pca: proto_pca.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

dist_svd: dist_svd.c $(FILES)
	$(MPICC) $(CFLAGS) $(INCS) $(LIBS) $(LINKS) -o $@ $^

clear:
	git clean -i

clean:
	rm -rf $(PROGS) *.dSYM *.o *.out *.mtx *.diag

