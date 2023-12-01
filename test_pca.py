import numpy as np
from scipy.io import mmread
from sklearn.decomposition import PCA

def read_diag(fname):
    return np.array([np.double(line.rstrip()) for line in open(fname)])

components = mmread("princomps.mtx")
explained_variance = read_diag("expvar.txt")
p = len(explained_variance)

X = mmread("A_out.mtx")
pca = PCA(n_components=p, svd_solver="full")
t = pca.fit_transform(X)

print(f"{np.linalg.norm(np.abs(pca.components_) - np.abs(components)):.18e}")
print(f"{np.linalg.norm(pca.explained_variance_ - explained_variance):.18e}")
