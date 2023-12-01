import numpy as np
from scipy.io import mmread

def read_diag(fname):
    return np.array([np.double(line.rstrip()) for line in open(fname)])

A = mmread("A_out.mtx")
U = mmread("U_out.mtx")
Vt = mmread("Vt_out.mtx")
S = read_diag("S_out.diag")

Atest = U@np.diag(S)@Vt

print(f"A=U*S*Vt is {np.allclose(A, Atest)} with distance {np.linalg.norm(A-Atest):.10e}")
