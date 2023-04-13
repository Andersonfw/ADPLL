#  Copyright (c) 2023. Vieira Filho Tecnológia

"""
Created on março 20 09:04:46 2023
@author: Ânderson F.W
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import ndarray
from scipy import signal

# Definir as constantes do TDC

#Para demostração
delay = 0.1 # tempo de delay
f_ref = 1 / 1.6 # frequência de referência
f_clk = 2 * f_ref # frequência do clock DCO
n_bits = int(1 / f_clk / delay)  # número de bits de saída / estágios do TDC
time_simulation = 1.5 * 1 / f_ref  # tempo de simulação do TDC

# delay = 1e-12  # tempo de delay
# f_ref = 40e6  # frequência de referência
# f_clk = 2.4e9  # frequência do clock DCO
# n_bits = int(1 / f_clk / delay)  # número de bits de saída / estágios do TDC
# time_simulation = 1.5 * 1 / f_ref  # tempo de simulação do TDC

print("número de bits:", n_bits)

# Simular o TDC para um sinal de entrada quadrado
t = np.linspace(0, time_simulation, 1000)  # vetor de tempo
v_clk = np.where(np.sin(2 * np.pi * f_clk * t) > 0, 1, -1)  # sinal de entrada CLK DCO
v_ref = signal.square(2 * np.pi * f_ref * t - 0 * np.pi / 2)  # sinal de referẽncia FREF deslocado em 90 graus

v_clk_delayed: ndarray = np.zeros((n_bits, len(t)))
# cria uma matriz de n_bits vetores, cada um com um deslocamento definidido pelo valor do delay;
for i in range(n_bits):
    v_clk_delayed[i, :] = signal.square((2 * np.pi * f_clk * t) + (i * np.pi / (n_bits / 2)))

tdc_bits = np.zeros(n_bits)  # array de saída do valor binario;
for i in range(len(t)):
    if v_ref[i - 1] < 0 < v_ref[i]:  # if v_ref[i - 1] < 0 and v_ref[i] > 0:
        for n in range(n_bits):
            tdc_bits[n] = np.where(v_clk_delayed[n_bits - 1 - n, i] > 0, 1, 0)

print("valor de sáida do TDC em bits: ", tdc_bits)

# Plotar os resultados
norm = 1e9  # valor de normalização do tempo
time = "ns"  # string de texto do eixo de tempo

# plot dos resultados binários
if n_bits < 10:
    fig, ax = plt.subplots(n_bits + 2, 1)  # plotar sinais
    plt.title("TDC binário")
    ax[0].plot(t * norm, v_clk, 'r', label="CLK")
    ax[0].legend()
    ax[1].plot(t * norm, v_ref, 'g', label="REF")
    ax[1].legend()
    for i in range(n_bits):
        ax[i + 2].plot(t * norm, v_clk_delayed[n_bits - i - 1, :], label='D={} '.format(i))
        ax[i + 2].legend()
    plt.xlabel('Tempo ({})'.format(time))
    plt.ylabel('Saída do TDC (bits)')

# Reconstrução do sinal a partir do valor do TDC
print('Frequência do sinal TDC:{}'.format(1 / (delay * n_bits) / norm), 'GHz')
amp_tdc = np.zeros(2 * n_bits)  # vetor de amplitudes do sinal com o dobro do tamanho para mostrar dois peridos
amp_tdc[0:n_bits] = tdc_bits  # copia o valor binarios para a primeira metade
amp_tdc[n_bits:] = tdc_bits  # copia o valor binario para a segunda metade
t_r = np.arange(0, 2 * n_bits, 1)  # gera o vetor de tempo
plt.figure()
plt.title("Sinal reconstruido pelo TDC")
plt.xlabel('Tempo ({})'.format(time))
plt.ylabel('Amplitude')
plt.plot(t_r * delay * norm, amp_tdc)  # plota o gráfico normalizando o tempo pelo delay definido
plt.show()
