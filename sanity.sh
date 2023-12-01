#!/bin/bash

rm -f *.mtx *.diag

./gen_svd -m 256 -n 128 -u -1 -c 100 -d 2 -o sanity

mpirun -np 8 ./dist_svd A_sanity.mtx S_sanity.diag U_sanity.mtx Vt_sanity.mtx 10
