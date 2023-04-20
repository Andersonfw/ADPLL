from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from deltasigma import *

OSR = 32
order = 5
H = synthesizeNTF(order, OSR, 1)
def fun_calc_psd(x, fs=1, rbw=100e3, fstep=None):
    # Calculate power spectral density
    #
    # INPUT arguments
    # x     :  data vector
    # fs    :  sample rate [Hz]
    # rbw   :  resolution bandwidth [Hz]
    # fstep :  FFT frequency bin separation [Hz]
    # OUTPUT
    # XdB	: spectrum of x [dB]
    # f	: frequency vector [Hz]

    if fstep is None:
        fstep = rbw / 1.62
    len_x = len(x)
    nwin = round(fs * 1.62 / rbw)
    nfft = round(fs / fstep)
    if nwin > len_x:
        nwin = len_x
        rbw = fs * 1.62 / nwin
    fftstr = f'len(x)={len_x:.2f}, rbw={rbw / 1e3:.2f}kHz, fstep={fstep / 1e3:.2f}kHz, nfft={nfft:d}, nwin={nwin:d}'
    print(f'Calculating the PSD: {fftstr} ...')
    f, X = signal.welch(x, fs=fs, window=signal.windows.blackman(nwin), nperseg=nwin, nfft=nfft, scaling='density')
    X *= (np.sinc(f / fs)) ** 2  # correct for ZOH
    XdB = 10 * np.log10(X)
    XdB_sig = np.max(XdB)
    print(f'Signal PSD peak = {XdB_sig:.2f} dB, 10log(rbw) = {10 * np.log10(rbw):.1f}')
    return XdB, f

# Parâmetros do modulador delta-sigma
fs = 100e3 # Frequência de amostragem (Hz)
f0 = 1e3 # Frequência do sinal analógico de entrada (Hz)
A = 1.0 # Amplitude do sinal analógico de entrada (V)
order = 1 # Ordem do modulador delta-sigma

# Geração do sinal analógico de entrada
t = np.arange(0, 400*1/f0, 1/fs) # Vetor de tempo
x = A * np.sin(2 * np.pi * f0 * t) # Sinal analógico de entrada (seno)
# Set the length of the white noise signal
length = len(x)
# Set the mean and standard deviation of the white noise distribution
mean = 0.2
std_dev = 0.5

# Generate the white noise signal
white_noise = np.random.normal(mean, std_dev, length)
x = x * white_noise
# Quantização delta
quantizer = np.zeros_like(x) # Inicializa o vetor de saída do quantizador
# for i in range(len(x)):
#     if x[i] >= 0:
#         quantizer[i] = A
#     else:
#         quantizer[i] = -A
N = 3

res_quantization = A * 2 / 2 ** N

# x_n = np.round(x_t / res_quantization) * res_quantization
# divide o valor de cada amostra de x_n pela resolução de N bits, e pega o valor arredondado. Será gerado valores de 0 até 8 para N = 3
x_q = np.round(x / res_quantization)
# Do valor arredondado multiplica pela resolução para normalizar a amplitude V novamente.
x_q = x_q * res_quantization

# Modulação delta-sigma
delta_sigma_output = np.zeros_like(x) # Inicializa o vetor de saída do modulador delta-sigma
prev = 0 # Valor anterior do modulador delta-sigma
sum = 0
out = 0
ADC = 1.1
for i in range(len(t)):
    sum = x_q[i] - prev # Cálculo do delta
    out += sum # Atualiza o valor anterior do modulador delta-sigma
    if out >= 0:
        quantizer[i] = A
        prev = ADC
    else:
        quantizer[i] = 0
        prev = -ADC


XdB = np.fft.fft(quantizer) / len(quantizer)
fft_q = np.fft.fft(x_q) / len(x_q)

N_BINS = len(quantizer)
F_RES = fs / N_BINS
f = np.arange(0, N_BINS * F_RES, F_RES)

XdB_q, fq = fun_calc_psd(quantizer, fs, 100)
XdB_o, fo = fun_calc_psd(x, fs, 100)

plt.figure()
plt.plot(fo, XdB_o, label="sinal quantizado")
# plt.semilogy(fo/1e6, XdB_o,label="full resolution")
plt.plot(fq, XdB_q, label="Saida do modulador DELTA")
# plt.semilogy(fq/1e6, XdB_q,label=""full resolution")
# plt.plot(fq / 1e6, np.full_like(fq, np.mean(XdB_q)), 'black', label="mean")
plt.grid(visible=True)
plt.legend()
plt.xlabel('Feq (MHz)')
plt.ylabel('Amplitude')

# Plot do sinal analógico de entrada, sinal quantizado e saída do modulador delta-sigma
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(t, x)
plt.plot(t, x_q)
plt.title('Sinal Analógico de Entrada')
plt.xlabel('Tempo (s)')
plt.ylabel('Amplitude (V)')
plt.subplot(3, 1, 2)
plt.step(t, quantizer, where='post')
plt.title('Sinal Quantizado (Delta)')
plt.xlabel('Tempo (s)')
plt.ylabel('Amplitude (V)')
plt.subplot(3, 1, 3)
plt.step(f, abs(XdB),"r",label="fft DELTA")
plt.step(f, abs(fft_q), "black", label="fft quantizado")
plt.title('Saída do Modulador Delta-Sigma')
plt.xlabel('Tempo (s)')
plt.ylabel('Amplitude (V)')
plt.grid(visible=True)
plt.legend()
plt.tight_layout()
plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import deltasigma
# OSR = 32
# order = 5
# H = deltasigma.synthesizeNTF(order, OSR, 1)
#
#
# plt.figure(figsize=(20, 4))
# N = 8192
# fB = int(np.ceil(N/(2.*OSR)))
# ftest = np.floor(2./3.*fB)
# u = 0.5*np.sin(2*np.pi*ftest/N*np.arange(N))
# v, xn, xmax, y = deltasigma.simulateDSM(u, H)
# t = np.arange(301)
# plt.step(t, u[t],'r')
# # plt.hold(True)
# plt.step(t, v[t], 'g')
# plt.axis([0, 300, -1.2, 1.2])
# plt.xlabel('Sample Number')
# plt.ylabel('u, v')
# plt.title('Modulator Input & Output');