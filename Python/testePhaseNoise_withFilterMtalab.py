import math

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import decimal
import locale
import datetime
import time
import os
import glob
import pandas as pd
import csv
def plot_DCO_signal():
    '''
    plotar sinal de saída do DCO e comparar ao sinal de inicio

    INPUT arguments
    none
    '''
    global len_simulation , dco_init_time , dco_freq_init , f_CKV , T0 , OVERSAMPLE , fs
    fs = OVERSAMPLE * f_CKV
    t = np.arange(0 , 10 * T0 , 1 / fs)
    # x_noise = np.sin(2 * np.pi * 1/(tc + jitter + wander) * t)
    x = np.sin(2 * np.pi * f_CKV * t)
    plt.figure()
    plt.plot(dco_init_time[:len_simulation] / 1e-9 , dco_freq_init[:len_simulation] , label="DCO Inicial")
    plt.plot(t[:len_simulation] / 1e-9 , x[:len_simulation] , label="DCO out")
    plt.grid(visible=True)
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (V)')
    plt.show()

def plot_phaseNoiseMask():
    # Geração das frequências
    frequencies_1k_to_3_5M = np.linspace(1000, 3.5e6, 300)  # De 1k a 3,5M
    frequencies_3_5M_to_10M= np.linspace(3.5e6, 10e6, num=200)  # De 3,5M a 10M
    frequencies_above_10M = np.linspace(10e6, 2.2e9, num=300)  # Acima de 10M

    # Geração dos níveis de fase noise correspondentes
    phase_noise_1k_to_1M = np.full_like(frequencies_1k_to_3_5M, -89)
    phase_noise_1M_to_10M = np.full_like(frequencies_3_5M_to_10M, -124)
    phase_noise_above_10M = np.full_like(frequencies_above_10M, -132)

    # Concatenação dos dados
    freq = np.concatenate((frequencies_1k_to_3_5M, frequencies_3_5M_to_10M, frequencies_above_10M))
    phase_noise = np.concatenate((phase_noise_1k_to_1M, phase_noise_1M_to_10M, phase_noise_above_10M))

    return phase_noise, freq
    
def fun_calc_psd(x , fs=1 , rbw=100e3 , fstep=None):
    '''
    Calculate power spectral density

    INPUT arguments
    x     :  data vector
    fs    :  sample rate [Hz]
    rbw   :  resolution bandwidth [Hz]
    fstep :  FFT frequency bin separation [Hz]
    OUTPUT
    XdB	: spectrum of x [dB]
    f	: frequency vector [Hz]
    '''

    if fstep is None:
        fstep = rbw / 1.62
    len_x = len(x)
    nwin = round(fs * 1.62 / rbw)
    nfft = round(fs / fstep)
    if nwin > len_x:
        nwin = len_x
        rbw = fs * 1.62 / nwin
    num_segments = 4
    nwin = math.floor(len(x) / num_segments)
    fftstr = (f'len(x)={len_x:.2f}, rbw={rbw / 1e3:.2f}kHz, fstep={fstep / 1e3:.2f}kHz, nfft={nfft:d}, nwin={nwin:d}')
    print(f'Calculating the PSD: {fftstr} ...')
    f , X = signal.welch(x , fs=fs , window=signal.windows.blackman(nwin) , nperseg=nwin , nfft=nfft ,
                         scaling='density')
    # f, X = signal.welch(x, nperseg=nwin, scaling='density')
    X *= (np.sinc(f / fs)) ** 2  # correct for ZOH
    XdB = 10 * np.log10(X)
    XdB_sig = np.max(XdB)
    print(f'Signal PSD peak = {XdB_sig:.2f} dB, 10log(rbw) = {10 * np.log10(rbw):.1f}')
    return XdB , f

filephaseFilteresd = 'C:/Users/ander/OneDrive - Associacao Antonio Vieira/UNISINOS/TCC/Python/phase_filtered.csv'
# x =[]
# with open(filephaseFilteresd, 'r') as csvfile:
#     a =list(csv.reader(filephaseFilteresd, quoting = csv.QUOTE_NONNUMERIC, delimiter = ' '))
#     for row in a:
#         x.append(row[:])
x = pd.read_csv(filephaseFilteresd, header=None)
    # x.append(a)
# x = np.array(a)
x = x.to_numpy()
print("cálculo da PSD")
Xdb_o , f = fun_calc_psd(x[0] , 2e9 , 2e3 , 700)
mask_phase_noise, freq  = plot_phaseNoiseMask() # obter a mascara de phase noise
plt.figure()
plt.semilogx(f , Xdb_o , label="Phase Noise")
plt.semilogx(freq, mask_phase_noise, label='Phase Noise MASK')
plt.grid(visible=True)
plt.legend()
plt.xlabel('Freq (Hz)')
plt.ylabel('Phase Noise [dBc/Hz]')
plt.show()