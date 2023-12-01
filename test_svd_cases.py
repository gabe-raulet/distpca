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
        cmd = f"srun -n 8 -N 1 -c 2 --cpu_bind=cores ./dist_svd {Afile} {Sfile} {Ufile} {Vtfile} {p}"
        print(cmd)
        proc = sp.Popen(cmd.split(), stdout=sp.PIPE)
        proc.wait()
