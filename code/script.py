#!/usr/bin/env python3

from scipy.optimize import minimize
import scipy.linalg as spla
import matplotlib.pyplot as plt
import pyfar as pf
from pathlib import Path

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

n = 800
u = pf.signals.exponential_sweep_time(n-ir.n_samples+1, (10, 22000))
y = pf.dsp.convolve(u, ir)
hexp = pf.dsp.deconvolve(y, u)

ir = pf.dsp.pad_zeros(ir, n-ir.n_samples)
hexp = pf.dsp.pad_zeros(hexp, n-hexp.n_samples)


def fun(x):
    x = pf.Signal(x, sampling_rate=44100)
    yx = pf.dsp.convolve(u, x, mode='cut', method='fft')
    #e = np.max(np.abs((y-yx).time))
    s = (y-yx).time
    e = spla.norm(s*np.exp(-np.linspace(0, 10, len(s))))
    return e

method = 'BFGS'
method = None
method = 'L-BFGS-B'
x = minimize(fun, np.squeeze(hexp.time)*.99, options={'disp': True}, method=method, tol=1e-6)
hopt = pf.Signal(x['x'], sampling_rate=44100)
pf.plot.time(ir-hopt)
pf.plot.time(ir-hexp)
#plt.plot(ir.times, np.exp(-np.linspace(0,10,n)))
plt.savefig('plot.png')
