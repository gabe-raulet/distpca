#!/bin/bash

rm -f *.mtx *.diag

echo "./gen_svd -m 1024 -n 1024 -u -1 -c 100 -d 2 -o sanity"
./gen_svd -m 1024 -n 1024 -u -1 -c 100 -d 2 -o sanity

echo "srun -n 32 -c 2 --cpu_bind=cores ./dist_svd sanity_A.mtx sanity_S.diag sanity_U.mtx sanity_Vt.mtx 10"
srun -n 32 -c 2 --cpu_bind=cores ./dist_svd sanity_A.mtx sanity_S.diag sanity_U.mtx sanity_Vt.mtx 10
