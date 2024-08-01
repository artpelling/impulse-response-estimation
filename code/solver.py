#!/usr/bin/env python3
import numpy as np
import pyfar as pf
import scipy.linalg as spla
import scipy.optimize as spopt
from pymor.operators.numpy import NumpyToeplitzOperator, NumpyMatrixOperator


def solve_regularized_toeplitz_LS(u, y, **solver_kwargs):
    assert u.sampling_rate == y.sampling_rate
    sampling_rate = u.sampling_rate

    x0 = np.squeeze(pf.dsp.deconvolve(y, u).time)
    #x0 = np.zeros_like(x0)

    u, y = np.squeeze(u.time), np.squeeze(y.time)
    assert u.ndim == y.ndim == 1

    # construct Toeplitz operator
    c = np.zeros_like(y)
    c[:len(u)] = u
    r = np.zeros_like(y)
    r[0] = c[0]
    K = NumpyToeplitzOperator(c, r=r)
    f = pf.dsp.filter.butterworth(pf.signals.impulse(x0.shape[0]), 2, 2000, btype='highpass').time.T
    f = f.reshape(-1, 1, 1)

    c = np.zeros_like(f)
    c[0] = f[0]
    F = NumpyToeplitzOperator(c=c, r=f)  # upper triangular F

    import scipy.sparse as sps
    T60 = 0.9
    alpha = 1e-3**(1/T60/44100)

    alpha = 0.6
    print(alpha)
    n = len(y)
    ai = alpha**(np.arange(n)+1)
    d = (ai*(1-alpha))**-0.5
    d[-1] = ai[-1]**-0.5
    d *= 1
    m = 1e6
    d[d > m] = m
    print(d)
    d /= np.max(d)
    D = NumpyMatrixOperator(sps.dia_array((d, 0), shape=(n, n)))

    # define cost function
    def fun(x):
        xpm = K.source.from_numpy(x)
        Kx = K.apply(xpm).to_numpy().T
        f = D.apply(F.apply(xpm)).to_numpy().T
        return 0.5*(spla.norm(Kx-y)**2 + f.T@f)

    # define jacobian
    def jac(x):
        xpm = K.source.from_numpy(x)
        ypm = K.range.from_numpy(y)
        J = K.apply_adjoint(K.apply(xpm)-ypm) + F.apply_adjoint(D.apply_adjoint(D.apply(F.apply(xpm))))
        return np.squeeze(J.to_numpy().T)

    x0 = np.squeeze(K.source.random(1, distribution='normal').to_numpy().T)

    x = spopt.minimize(fun, x0, jac=jac, **solver_kwargs)['x']

    return pf.Signal(x, sampling_rate=sampling_rate)
