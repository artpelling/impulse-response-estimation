#!/usr/bin/env python3

from scipy.optimize import minimize
import scipy.linalg as spla
import matplotlib.pyplot as plt
import pyfar as pf
from pathlib import Path

from solver import solve_regularized_toeplitz_LS


plt.close('all')

files = {
    'basilica': 'KE_Basilica_of_Eberbach_Monastery.wav',  # T=4.6
    'gewandhaus': 'GH_Gewandhaus.wav',  # T=1.9
    'chamber': 'KE_Chamber_Music_Hall_of_Konzerthaus_Berlin.wav',  # T=1.2
    'murakuni': 'MK_Murakuni-Za.wav',   # T=0.9
}

file = Path('../resources/rirs_simulated') / files['murakuni']
ir = pf.io.read_audio(file)

file = Path('/home/pelling/data/fabian/SOFA/FABIAN_HRIR_measured_HATO_0.sofa')
ir = pf.io.read_sofa(file)[0][0,0]

n = 500
u = pf.signals.exponential_sweep_time(n-ir.n_samples+1, (10, 22000))
y = pf.dsp.convolve(u, ir)
hexp = pf.dsp.deconvolve(y, u)

ir = pf.dsp.pad_zeros(ir, n-ir.n_samples)
hexp = pf.dsp.pad_zeros(hexp, n-hexp.n_samples)

options = {'disp': True, 'gtol': 1e-8}
hopt = solve_regularized_toeplitz_LS(u, y, tol=0, method='BFGS', options=options)

pf.plot.time(ir-hexp)
pf.plot.time(ir-hopt)
#pf.plot.time(hopt)
#plt.xlim(np.array((-1,5))/44100)
#plt.ylim(np.array((-0.3,0.3))*1e-5)
#plt.plot(ir.times, np.exp(-np.linspace(0,10,n)))
plt.savefig('plot.png')
