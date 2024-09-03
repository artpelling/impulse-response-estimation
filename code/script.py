#!/usr/bin/env python3

import numpy as np
import scipy.linalg as spla
import scipy.signal as ss
import matplotlib.pyplot as plt
import pyfar as pf
from pathlib import Path

from pymor.algorithms.to_matrix import to_matrix

from solver import solve_regularized_toeplitz_LS
from kernel import construct_toeplitz_operator


plt.close('all')

file = Path('/home/pelling/data/fabian/SOFA/FABIAN_HRIR_measured_HATO_0.sofa')
ir = pf.io.read_sofa(file)[0][0,0]

n_sweep = 512
u = pf.signals.exponential_sweep_time(n_sweep, (100, 22000))
y = pf.dsp.convolve(u, ir)
hexp = pf.dsp.deconvolve(y, u)

n_taps = 256
K = construct_toeplitz_operator(u.time[0], y.n_samples, n_taps)
k = to_matrix(K)
U, sv, Vh = spla.svd(k, full_matrices=False)
plot_style = {"lw": 1, "fillstyle": 'none', 'markersize': 5}
fig, ax = plt.subplots(figsize=(16/3, 3), dpi=300, constrained_layout=True)
ax.semilogy(sv, '.-', label='$\sigma_i(K)$', **plot_style)
uTy = y.time[0]@U
ax.semilogy(np.abs(uTy), '.-', label='$|u_i^Ty|$', **plot_style)
ax.semilogy(np.abs(uTy)/sv, '.-', label='$|u_i^Ty|/\sigma_i$', **plot_style)
ax.set_xlim((0, n_taps))
plt.legend()
plt.savefig('svs.png')


# ALGO
alpha = 0.825
frq = (100, 8000)
beta = ss.firls(101, np.array([0, frq[0], frq[0], frq[1], frq[1], 22050]), np.array([1, 1, 0, 0, 1, 1]), fs=44100)
beta = ss.minimum_phase(beta)
regularization = {'lam': 0, 'alpha': alpha, 'beta': beta}

# method can be either 'lstsq' or 'BFGS'
solver_kwargs = {'init': None, 'method': 'lstsq', 'tol': 1e-12, 'options': {'disp': True, 'gtol': 1e-8}}

hopt = solve_regularized_toeplitz_LS(u, y, n_taps, **solver_kwargs, **regularization)


# PLOTTING
dB = True
if n_taps > hexp.n_samples:
    ir = pf.dsp.pad_zeros(ir, hopt.n_samples-ir.n_samples)
    hexp = pf.dsp.pad_zeros(hexp, hopt.n_samples-hexp.n_samples)
else:
    ir = pf.dsp.pad_zeros(ir, hexp.n_samples-ir.n_samples)
    hopt = pf.dsp.pad_zeros(hopt, hexp.n_samples-hopt.n_samples)

fig, ax = plt.subplots(figsize=(16/3, 3), dpi=300, constrained_layout=True)
pf.plot.time((ir-hexp)/ir, dB=dB, unit='samples', label='pyfar', **plot_style)
pf.plot.time((ir-hopt)/ir, dB=dB, unit='samples', label=solver_kwargs['method'], **plot_style)
ax.set_xlim((0, n_taps))
plt.legend()
plt.savefig('errtime.png')

fig, ax = plt.subplots(figsize=(16/3, 3), dpi=300, constrained_layout=True)
pf.plot.freq(ir-hexp, dB=dB, label='pyfar', **plot_style)
pf.plot.freq(ir-hopt, dB=dB, label=solver_kwargs['method'], **plot_style)
plt.legend()
plt.savefig('errfreq.png')
