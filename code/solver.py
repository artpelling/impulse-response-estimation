#!/usr/bin/env python3
import numpy as np
import pyfar as pf
import scipy.linalg as spla
import scipy.optimize as spopt
from pymor.operators.numpy import NumpyToeplitzOperator


def solve_regularized_toeplitz_LS(u, y, **solver_kwargs):
    assert u.sampling_rate == y.sampling_rate
    sampling_rate = u.sampling_rate

    x0 = np.squeeze(pf.dsp.deconvolve(y, u).time)

    u, y = np.squeeze(u.time), np.squeeze(y.time)
    assert u.ndim == y.ndim == 1

    # construct Toeplitz operator
    c = np.zeros_like(y)
    c[:len(u)] = u
    r = np.zeros_like(y)
    r[0] = c[0]
    K = NumpyToeplitzOperator(c, r=r)

    # define cost function
    def fun(x):
        xpm = K.source.from_numpy(x)
        Kx = K.apply(xpm).to_numpy().T
        return spla.norm(Kx-y)**2/2

    # define jacobian
    def jac(x):
        xpm = K.source.from_numpy(x)
        ypm = K.range.from_numpy(y)
        return np.squeeze(K.apply_adjoint(K.apply(xpm)-ypm).to_numpy())

    #x0 = np.squeeze(K.source.random(1, distribution='normal').to_numpy().T)

    x = spopt.minimize(fun, x0, jac=jac, **solver_kwargs)['x']

    return pf.Signal(x, sampling_rate=sampling_rate)
