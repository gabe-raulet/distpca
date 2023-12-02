import sys
import numpy as np
import subprocess as sp
from pathlib import Path
import glob
from scipy.io import mmread, mmwrite

p = 10

if __name__ == "__main__":

    target = Path("svd_matrix_cases")
    assert target.exists() and target.is_dir()

    for Afile in glob.glob(f"{target.name}/*_A.mtx"):
        label = Afile.split("A.mtx")[0]
        Sfile = f"{label}S.diag"
        Ufile = f"{label}U.mtx"
        Vtfile = f"{label}Vt.mtx"
        assert Path(Afile).is_file() and Path(Sfile).is_file() and Path(Ufile).is_file() and Path(Vtfile).is_file()

        m = int(Afile.split("_rows")[1].split("_")[0])
        n = int(Afile.split("_cols")[1].split("_")[0])

        tcpunch = sp.Popen(["fmi/extern/TCPunch/server/build/tcpunchd"], stdout=sp.PIPE, stderr=sp.PIPE)

        nprocs = 8
        procs = [None]*nprocs

        for rank in range(nprocs):
            # cmd = f"./build/fmi_svd {rank} {nprocs} {Afile} {Sfile} {Ufile} {Vtfile} {p}"
            cmd = f"./build/fmi_svd {rank} {nprocs} {m} {n} {Afile} {Sfile} {Ufile} {Vtfile} {p}"
            procs[rank] = sp.Popen(cmd.split(), stdout=sp.PIPE)

        for rank in range(nprocs):
            procs[rank].wait()
            if procs[rank].returncode != 0:
                sys.stderr.write("./build/fmi_svd exited with non-zero status\n")
                sys.stdout.flush()
                sys.exit(1)

        tcpunch.terminate()
        print("\n")
