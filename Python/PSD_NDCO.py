#  Copyright (c) 2023. Ânderson Felipe Weschenfelder

"""
Created on abril 01 21:23:40 2023

@author: Ânderson Felipe Weschenfelder
"""
from _ast import expr

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


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


N = 3
fs = 450e6
ts = 1 / fs
f1 = 1.33e6
V = 1
t = np.arange(0, 400 * 1 / f1, ts)
# resolução de cada nível. Para N = 3bits são 8 níveis, e se V=1V, cada nível vale 1/8V.,
res_quantization = V * 2 / 2 ** N

x_t = V * np.sin(2 * np.pi * f1 * t)
# x_n = np.round(x_t / res_quantization) * res_quantization
# divide o valor de cada amostra de x_n pela resolução de N bits, e pega o valor arredondado. Será gerado valores de 0 até 8 para N = 3
x_n = np.round(x_t / res_quantization)
# Do valor arredondado multiplica pela resolução para normalizar a amplitude V novamente.
x_n = x_n * res_quantization

# Calculando o erro
err = x_t - x_n

# Quanto mais níveis melhor a resolução e menor o erro de quantização

# Plotando o sinal original, quantizado e o erro de quantização
fig, axs = plt.subplots(3, 1, figsize=(10, 8))

axs[0].plot(t[:1000], x_t[:1000])
axs[0].set_title('Sinal Original')

axs[1].step(t[:1000], x_n[:1000], where='post')
axs[1].set_title('Sinal Quantizado (3 bits)')

axs[2].plot(t[:1000], err[:1000])
axs[2].set_title('Erro de Quantização')

plt.tight_layout()  # ajustar automaticamente as posições das subplots de um gráfico para que não haja sobreposição de legendas, rótulos de eixos ou títulos.
# plt.show()

# XdB_o, fo = fun_calc_psd((x_t), fs, 100e3)
# XdB_q, fq = fun_calc_psd((x_n), fs, 100e3)

XdB_o, fo = fun_calc_psd((x_t), fs, 100e3)
XdB_q, fq = fun_calc_psd((x_n), fs, 100e3)

# plt.figure()
# f , x = signal.welch(x_t, fs, nperseg=1024)
# plt.semilogy(f, x)

plt.figure()
plt.plot(fo / 1e6, XdB_o, label="Full resolution")
# plt.semilogy(fo/1e6, XdB_o,label="full resolution")
plt.plot(fq / 1e6, XdB_q, label="Quantized")
# plt.semilogy(fq/1e6, XdB_q,label=""full resolution")
plt.plot(fq / 1e6, np.full_like(fq, np.mean(XdB_q)), 'black', label="mean")
plt.grid(visible=True)
plt.legend()
plt.xlabel('Feq (MHz)')
plt.ylabel('Amplitude')
# plt.show()

Xq_value_raw = 10 ** (XdB_q / 10)   # PSD do sinal quantizado em valores de pico não em DB
Xt_value_raw = 10 ** (XdB_o / 10)   # PSD do sinal original em valores de pico não em DB
# X_value_q = XdB_q
# X_value_t = XdB_o

# Identificar a banda de frequência do sinal
f_signal_min = f1/2  # margem esquerda da banda
f_signal_max = f1 + f1/2  # margem direita da banda

# Identificar a banda de frequência do ruído
f_noise_min = 0   # margem esquerda da banda
f_noise_max = fs/2   # margem direita da banda

# Calcular a PSD média do sinal na banda de interesse
Pxx_signal = np.mean(Xq_value_raw[(fq >= f_signal_min) & (fq <= f_signal_max)])

# Calcular a potência do sinal na banda de interesse
# Potência téorica e dada por V^2 /2
Ps_calc = V**2 / 2
print("potência do sinal calculado", Ps_calc)
Ps = np.trapz(Xq_value_raw[(fq >= f_signal_min) & (fq <= f_signal_max)], fq[(fq >= f_signal_min) & (fq <= f_signal_max)])
print("potência do sinal medida", Ps)

# Calcular a potência do ruído fora da banda de interesse
Pn = np.trapz(Xq_value_raw[(fq >= f_noise_min) & (fq <= f_noise_max)], fq[(fq >= f_noise_min) & (fq <= f_noise_max)])
# print("pn total:", Pn)
Pn = (Pn - Ps)
print("potência do ruido", Pn)

# Calcular a SNR em dB
SNR = 10 * np.log10(Ps / Pn)
# SNR em dBc/Hz (decibéis em relação à potência do sinal no canal por Hertz)
SNR_dbc = - SNR - 10* np.log10(fs)
print('SNR: {:.2f} DB'.format(SNR))
print('SNR: {:.2f} dBc/Hz'.format(SNR_dbc))

# SNR em dB teórico
SNR_Calc = 6.02 * N + 1.76

# SNR em dBc teórico
SNR_Calc_dbc = - SNR_Calc - 10 * np.log10(fs)
print('SNR Teorico: {:.2f} DB'.format(SNR_Calc))
print('SNR Teorico: {:.2f} dBc/Hz'.format(SNR_Calc_dbc))
