#!/usr/bin/env python3
import numpy as np
import pyfar as pf
import scipy.linalg as spla
import scipy.optimize as spopt

from pymor.algorithms.to_matrix import to_matrix

from kernel import construct_filter_operator, construct_toeplitz_operator


def solve_regularized_toeplitz_LS(u, y, n_taps, lam=1, alpha=0.9, beta=np.array([1]), init=0, method='qr', **solver_kwargs):
    assert u.sampling_rate == y.sampling_rate
    sampling_rate = u.sampling_rate

    # starting solution 0 or deconv
    x0 = np.squeeze(pf.dsp.deconvolve(y, u).time)[:n_taps]
    if init == 0:
        x0 = np.zeros_like(x0)
    elif init == 'rand':
        x0 = np.random.randn(*x0.shape)
    n = x0.shape[0]

    u, y = np.squeeze(u.time), np.squeeze(y.time)
    assert u.ndim == y.ndim == 1

    # construct Toeplitz operator
    K = construct_toeplitz_operator(u, len(y), n_taps)
    
    F = lam*construct_filter_operator(n_taps, alpha, beta, delay=0)

    if method == 'lstsq':
        K, F = to_matrix(K), to_matrix(F)
        x = spla.lstsq(K.T@K + F, K.T@y)[0]
    else:
        # define cost function
        def fun(x):
            xpm = K.source.from_numpy(x)
            Kx = K.apply(xpm).to_numpy().T
            f = F.apply(xpm).to_numpy().T
            return 0.5*(spla.norm(Kx-y)**2 + spla.norm(f)**2)

        # define jacobian
        def jac(x):
            xpm = K.source.from_numpy(x)
            ypm = K.range.from_numpy(y)
            J = K.apply_adjoint(K.apply(xpm)-ypm) + F.apply_adjoint(F.apply(xpm))
            return J.to_numpy().squeeze()

        x = spopt.minimize(fun, x0, jac=jac, **solver_kwargs)['x']

    return pf.Signal(x, sampling_rate=sampling_rate)
