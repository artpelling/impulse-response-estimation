# %%
import pyfar as pf
import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
%matplotlib inline

# AK internal database (all RIRs simulated with RAVEN by David Ackermann)
path = '/Volumes/Forschungsprojekte/AK-Data/ROOMS'
# NOTE: RIR does not get shorter after room MK, because low freq reverberation
#       remains approx. constant.
rooms = [
    'KE_Basilica_of_Eberbach_Monastery/KE_RIR/KE_RIR_R1.wav',           # T=4.6
    'GH_Gewandhaus/GH_RIR/GH_RIR_R1.wav',                               # T=1.9
    'KH_Chamber_Music_Hall_of_Konzerthaus_Berlin/KH_RIR/KH_RIR_R1.wav', # T=1.2
    'MK_Murakuni-Za/MK_RIR/MK_RIR_R1.wav'                               # T=0.9
    #-  'HFT_Seminar_Room_HFT616_TU_Berlin/HFT_RIR/HFT_RIR_R1.wav',     # T=0.6
    #-  'ST_Studio/ST_RIR/ST_RIR_R1.wav',                               # T=0.3
    ]

# block length for analyzing rir length
n_block = 2048

# estimate RIR length: first block with a maximum absolute amplitude < 60 dB
# relative to the maximum value of the RIR
pf.plot.use()
for room_id, room in enumerate(rooms):

    # load RIR
    room_str = room.split('/')[0]
    rir = pf.io.read_audio(os.path.join(path, room))

    # reshape to blocks
    rir_max = np.argmax(np.abs(rir.time), -1)
    rir = pf.dsp.normalize(rir)
    rir = pf.dsp.pad_zeros(rir, n_block - rir.n_samples % n_block)
    blocks = np.reshape(rir.time, (-1, n_block))
    blocks_max = pf.dsp.power(pf.Signal(blocks, rir.sampling_rate))

    # estimate length
    rir_length = np.where(blocks_max < 10**(-90/10))[0] * n_block
    rir_length = rir_length[rir_length > rir_max + .1 * rir.sampling_rate][0]

    # truncate RIR
    rir_trunc = pf.dsp.time_window(
        rir, [rir_length, rir_length + n_block], shape='right', crop='end')

    # save T30 values
    room_str_id = room_str.split('_')[0]
    t30_file = os.path.join(
        path, os.path.dirname(room), '..', f'{room_str_id}_RAParameter',
        f'{room_str_id}_R1_T30.csv')
    shutil.copy(t30_file, f'{room_str}_T30.csv')

    # get T30_mid (average at 0.5 and 1 kHz)
    t30_file = f'{room_str}_T30.csv'
    t30_data = np.genfromtxt(t30_file, delimiter=',')[:, :2]
    t30_mid = np.mean(t30_data[2:3, 1])

    # save figure
    plt.figure()
    ax = pf.plot.time(rir, dB=True, color=[.5, .5, .5])
    ax = pf.plot.time(rir_trunc, dB=True, color='k')
    ax.set_ylim(-150, 5)
    ax.axvline(rir_length/rir.sampling_rate)
    ax.set_title((f'{room_str}\nT_m={t30_mid:.2f} s, '
                  f'duration {rir_length/rir.sampling_rate:.2f} s'))
    plt.savefig(f'{room_str}.png', bbox_inches='tight')

    # save wav
    pf.io.write_audio(rir_trunc, f'{room_str}.wav', 'FLOAT')
