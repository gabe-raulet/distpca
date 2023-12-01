MPICC=cc
INCS=-I./inc
LIBS=
LINKS=
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

test: gen_svd dist_svd
	python gen_svd_cases.py
	python test_svd_cases.py

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

