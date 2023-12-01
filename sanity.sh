#!/bin/bash

rm -f *.mtx *.diag

echo "./gen_svd -m 256 -n 128 -u -1 -c 100 -d 2 -o sanity"
./gen_svd -m 256 -n 128 -u -1 -c 100 -d 2 -o sanity

echo "srun -n 8 ./dist_svd sanity_A.mtx sanity_S.diag sanity_U.mtx sanity_Vt.mtx 10"
srun -n 8 ./dist_svd sanity_A.mtx sanity_S.diag sanity_U.mtx sanity_Vt.mtx 10
