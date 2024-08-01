#!/usr/bin/env python3

from pymor.operators.numpy import NumpyMatrixOperator, NumpyToeplitzOperator
from pymor.algorithms.to_matrix import to_matrix

import numpy as np
import scipy.linalg as spla

np.set_printoptions(precision=2, suppress=True)

n = 4
alpha = 0.9
ai = alpha**np.repeat(np.arange(n).reshape(n,1)+1, n, axis=1)
P = np.minimum(ai, ai.T)
print(f'P =\n{P}\n')

Pinv = spla.inv(P)
print(f'Pinv =\n{Pinv}\n')

r = np.zeros(n)
r[:2] = [1, -1]
c = np.zeros_like(r)
c[0] = r[0]
Fop = NumpyToeplitzOperator(c=c, r=r)

F = to_matrix(Fop)
print(f'F =\n{F}\n')

R = np.eye(n)*(1+alpha)
R[0,0], R[-1,-1] = 1, 1
R += -np.sqrt(alpha)*np.diag(np.ones(n-1), k=1)
R += -np.sqrt(alpha)*np.diag(np.ones(n-1), k=-1)
R /= (1-alpha)
R /= (ai*ai.T)**0.5
print(f'R =\n{R}\n')
print(f'Rinv =\n{spla.inv(R)}\n')
chol = spla.cholesky(R, lower=False)
print(f'Cholesky =\n{chol}\n')

d = (ai[:,0]*(1-alpha))**-0.5
D = np.diag(d)
D[-1, -1] = ai[-1,0]**-0.5
print(f'D =\n{D}\n')
DF = D@F
print(f'DF =\n{DF}\n')
