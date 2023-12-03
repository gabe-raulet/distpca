import sys
import numpy as np
import subprocess as sp
import time
from pathlib import Path
from scipy.io import mmread
from sklearn.decomposition import PCA

def getdims(fname):
    with open(fname, "r") as f:
        next(f)
        next(f)
        m, n, = (int(v) for v in next(f).rstrip().split())
    return m, n

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write(f"Usage: {sys.argv[0]} <A.mtx> <label> <p> <nprocs>\n")
        sys.stderr.flush()
        sys.exit(1)

    fname = sys.argv[1]
    label = sys.argv[2]
    p = int(sys.argv[3])
    nprocs = int(sys.argv[4])

    n, d = getdims(fname)

    t1 = -time.perf_counter()

    tcpunch = sp.Popen(["fmi/extern/TCPunch/server/build/tcpunchd"], stdout=sp.PIPE, stderr=sp.PIPE)

    procs = [None]*nprocs
    for rank in range(nprocs):
        cmd = f"./build/fmi_pca {rank} {nprocs} {n} {d} {fname} {label}_princomps.mtx {label}_expvar.txt {p}"
        procs[rank] = sp.Popen(cmd.split(), stdout=sp.PIPE)

    for rank in range(nprocs):
        procs[rank].wait()
        if procs[rank].returncode != 0:
            sys.stderr.write("./build/fmi_pca exited with non-zero status\n")
            sys.stdout.flush()
            sys.exit(1)

    tcpunch.terminate()
    t1 += time.perf_counter()
    print(f"[fmi_pca::n={n},d={d},p={p},nprocs={nprocs}] completed in {t1:.5f} seconds\n")

    X = mmread(fname)
    components = mmread(f"{label}_princomps.mtx")
    expvar = np.array([np.double(line.rstrip()) for line in open(f"{label}_expvar.txt")])
    assert len(expvar) == p

    t2 = -time.perf_counter()
    pca = PCA(n_components=p, svd_solver="full")
    t = pca.fit_transform(X)
    t2 += time.perf_counter()

    print(f"[sklearn.decomposition.PCA::n={n},d={d},p={p}] completed in {t2:.5f} seconds\n")

    print(f"princomps_err = {np.linalg.norm(np.abs(pca.components_) - np.abs(components)):.18e}")
    print(f"expvar_err = {np.linalg.norm(pca.explained_variance_ - expvar):.18e}")

