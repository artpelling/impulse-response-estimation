#!/usr/bin/env python3

import numpy as np
import scipy.sparse as sps

from pymor.operators.numpy import NumpyMatrixOperator, NumpyToeplitzOperator


def construct_toeplitz_operator(u, m, n):
    c = np.zeros(m)
    c[:len(u)] = u
    r = np.zeros(n)
    r[0] = c[0]
    return NumpyToeplitzOperator(c, r=r)


def construct_filter_operator(n, alpha, beta, delay=0):
    gamma = 1e-2
    dd = alpha**(-np.arange(n-delay)/2)
    d = sps.diags_array(np.concatenate([np.ones(delay)/gamma, dd/np.max(dd)]))
    D = NumpyMatrixOperator(d)
    r = np.zeros(n)
    r[:len(beta)] = beta
    c = np.zeros_like(r)
    c[0] = r[0]
    F = NumpyToeplitzOperator(c, r=r)
    return D@F
