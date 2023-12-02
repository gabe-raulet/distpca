import sys
import numpy as np
import subprocess as sp
from pathlib import Path
from scipy.io import mmread, mmwrite

def gen_params(rowcnts, colcnts, cvals, dvals, modes):
    params = []
    for m in rowcnts:
        for n in colcnts:
            for cond in cvals:
                for dmax in dvals:
                    for mode in modes:
                        params.append((m, n, cond, dmax, mode))
    return params

if __name__ == "__main__":

    target = Path("svd_matrix_cases")

    if target.exists():
        assert target.is_dir()
        for child in target.iterdir():
            child.unlink()
    else:
        target.mkdir()

    rowcnts = 2**np.arange(8,10)
    colcnts = 2**np.arange(8,10)
    cvals = [100]
    dvals = [2]
    modes = [-1, 2]
    params = gen_params(rowcnts, colcnts, cvals, dvals, modes)

    for m, n, cond, dmax, mode in params:
        label = str(target.joinpath(f"case_m{m}_n{n}_cond{cond}_dmax{dmax}_mode{mode}"))
        cmd = f"./gen_svd -m {m} -n {n} -u -1 -c {cond} -d {dmax} -u {mode} -o {label}"
        print(cmd)
        proc = sp.Popen(cmd.split(), stdout=sp.PIPE)
        proc.wait()
