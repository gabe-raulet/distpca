import numpy as np
from sklearn.decomposition import PCA
from scipy.io import mmread
from sklearn.preprocessing import StandardScaler

def mypca(X, p):
    m, n = X.shape
    if p <= 0: p = min(m,n)
    X = X.copy()
    mean = X.mean(axis=0)
    B = (X - mean)
    C = B.T@B / (m-1)
    D, V = np.linalg.eig(C) # V[:,i] is eigenvector for eigenvalue D[i], i.e. eigenvectors stored in columns.
    signs = np.sign(V[np.argmax(np.abs(V), axis=0), range(V.shape[0])])
    V = V*signs # make each column of V's largest component in absolute value be positive
    V = V.T # eigenvectors stored in rows now
    order = np.argsort(D)[::-1]
    D = D[order]
    V = V[order]
    components = V[:p,:]
    explained_variance = D[:p]
    explained_variance_ratio = explained_variance / np.sum(D)

    return components, explained_variance, explained_variance_ratio

def mypca2(X, p):
    m, n = X.shape
    if p <= 0: p = min(m,n)
    B = X - X.mean(axis=0)
    _, S, Vt = np.linalg.svd(B, full_matrices=False)
    D = (S*S)[:p] / (m-1) # explained variance
    Vt = Vt[:p,:] # principal components stored in rows
    Xproj = (B / X.std(axis=0))@Vt.T
    return Xproj, D, Vt


X = np.random.random((256,128))
m,n = X.shape
p = 10
components, explained_variance, explained_variance_ratio = mypca(X, p)


pca = PCA(n_components=p, svd_solver="full")
t = pca.fit_transform(X)

Xproj, D, Vt = mypca2(X, p)

C1 = components
C2 = pca.components_
E1 = explained_variance_ratio
E2 = pca.explained_variance_ratio_
F1 = explained_variance
F2 = pca.explained_variance_


