#  Copyright (c) 2023. Vieira Filho Tecnológia

"""
Created on março 17 14:17:40 2023
@author: Ânderson Felipe Weschenfelder
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

osrate = 8  # oversample ratio
freq = 100e6  # Signal Frequency
Fs = osrate * freq  # Frequency sampling rate
# sr = 1000e9
ts = 1 / Fs  # sampling interval
repeat = 10000  # number of times repeat the signal.
points_number = Fs * repeat / freq  # numero de amostras que será obtido do sinal
# em Fs(Hz) * "n" repetições tem a largura de banda, e divindo pela freq. do sinal tem-se o número de amostras
# Com Fs = 8GHz, em 10 repetições tem-se 80GHz, e sendo o sinal de 1GHz gerando 80 pontos

timesimulation = 1 / freq * repeat  # tempo total de simualção
Freq_resolution = 1 / timesimulation  # freq. de resolução da DFT. Para 10 repet. de 1GHz, a resolução será 100MHz

t = np.arange(0, timesimulation, ts)  # array of sampling com intervalo de tempo conforme a resolução definida
x = np.sin(2 * np.pi * freq * t)  # signal sampled
# x = x + np.sin(2 * np.pi * (1.1e9) * t)
# x = np.cos(2* np.pi*freq*t + 0.5*np.sin(2*np.pi*10e6*t))
f_triangular = 50e3
x1 = signal.sawtooth(2*np.pi*f_triangular*t, width=0.5)
x = x * x1

# plot signal sampled
plt.figure(figsize = (8, 6))
plt.plot(t, x, 'bo-')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.legend(["SR={} Hz".format(Fs)])
plt.show()

# calc fft of signal x
X = np.fft.fft(x) / len(x)  # calcula a fft do sinal e normaliza pelo número de pontos do sinal
N = len(X)  # tamanho da fft
k = np.arange(N)  # vetor em k do tamanho do sinal
T = N / Fs  # resolução de frequência Fs / N = Res.; cada N ponto representa N*Res (Hz)
# frq = np.fft.fftfreq(1000, ts)
frq = k / T  # k representa um vetor que multiplica pelo valor de resolução da FFT

# se for SSB
frq = frq[range(int(N / 2))]  # FFT é uma função par, logo só é necessário metade dos pontos
X = X[range(int(N / 2))]  # FFT é uma função par, logo só é necessário metade dos pontos
# X2 = np.fft.fftshift(X)

plt.figure(figsize=(12, 6))

# show the response in frequency of signal sampled
plt.subplot(141)
# plt.stem(freq, np.abs(X[np.arange(0,len(n),1)]), linefmt='b', markerfmt=' ',basefmt="-b")
plt.stem(frq, 2 * abs(X), linefmt='b', markerfmt=' ', basefmt="-b")  # multiplica por 2 o modúlo de X, pois sendo
# função par divide o sinal pela metade
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT SSB Amplitude |X(freq)|')
# plt.ylabel('FFT DSB Amplitude |X(freq)|')
# plt.xlim(0, 10)

# Power Spectre Density PSD
fft_x = np.fft.fft(x)  # Calculando a FFT do sinal
P = np.abs(fft_x) ** 2  # Calculando o espectro de potência
# Normalizando o espectro de potência
n = len(x)  # número de amostras
bw = Fs / 2  # largura de banda
P_norm = P / n / bw
f = np.linspace(0, bw, n // 2)  # cria vetor de frequências
#   Plot do gráfico
plt.subplot(142)
plt.semilogy(f, P_norm[:n // 2])
plt.xlabel('Freq (Hz)')
plt.ylabel('PSD (dB/Hz)')
plt.legend(['PSD with calc'])

# f, psd = signal.periodogram(x, fs=Fs)
# plt.semilogy(f, psd)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('PSD (dB/Hz)')
# plt.legend(['PSD with LIB'])


plt.subplot(143)

# Energy  Spectre Density ESD
esdx = P_norm[:n // 2] * ts  # time interval
plt.subplot(143)
plt.plot(f, np.log10(esdx))
plt.xlabel('Freq (Hz)')
plt.ylabel('ESD (dB s/Hz)')
plt.legend(['ESD with calc'])

# apply de inverse fft to reconstuct the base signal
plt.subplot(144)
y = np.fft.ifft(np.fft.fft(x))
plt.plot(t, np.abs(y), 'r')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.tight_layout()
plt.show()
