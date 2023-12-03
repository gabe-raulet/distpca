#!/bin/bash

P=$1
NPROCS=$2

python test_pca_fmi.py canonical_test/case_rows1024_cols512_cond100_dmax2_mode-1_A.mtx canonical_test/output $P $NPROCS
