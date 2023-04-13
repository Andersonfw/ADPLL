#  Copyright (c) 2023. Ânderson Felipe Weschenfelder

"""
Created on março 30 19:41:30 2023

@author: Ânderson Felipe Weschenfelder
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

ts = 100e-12  # time step of time vector
t_max = 250e-6  # max value of time vector
fc = 10e6  # carrier frequency
t0 = 1 / fc  # period of carrier
kvco = 1e6  # kvco MHz/V
t = np.arange(0, t_max, ts)  # time vector
fm = 100e3  # frequency modulation for vtune
vtune = 5  # np.cos(2 * np.pi * fm * t)  # tune control of DCO
hysteresis = 0  # hysteresis of rise edge detect

f0 = np.sin(2 * np.pi * fc * t)  # signal of carrier frequency
v_t = np.sin(2 * np.pi * fc * t + 2 * np.pi * kvco * vtune * t)  # output signal DCO
noise = np.random.normal(0, 0.1, len(v_t)) + np.cos(2 * np.pi * 3e6)  # noise signal
# v_t = v_t + noise

len_simulation = 10000  # size of vector to simulate faster
t_k = np.zeros(len_simulation)  # vector to save the time value when occur the rise edge of DCO clock
ckr = np.zeros(len_simulation)  # output sigan of rise and falling edge of DCO
delta_t = np.zeros(int(len_simulation * ts / t0))  # delta time between the ideal time instance and the DCO time
TDEV = np.zeros(int(len_simulation * ts / t0) + 1)  # time deviation accumulated
k = 0  # index for time vector
count = 0  # int(t0 / ts / 4)  # counter to identify when new period of carrier frequency start. O to sine, t0/ts/4 to cos
i = 0  # index of TDEV
start = 0   # detecção do primeiro cruzamento por zero, para correto calculo de delta_t, considerando ao menor um valor de T0/2
for n in range(len_simulation):
    if v_t[n] == 0 or (v_t[n - 1] < 0 < v_t[n]) or (v_t[n - 1] > 0 > v_t[n]):
        start = n * ts
        break
print('start',start/ts)
for n in range(len_simulation):
    if v_t[n - 1] < 0 and v_t[n] > hysteresis or (
            v_t[0] <= 0 and v_t[n] > hysteresis and n < 2):  # detectar se iniciar em zero e for borda de subida
        t_k[k] = n * ts - start
        k = k + 1
        ckr[n] = 1
    elif v_t[n - 1] > 0 and v_t[n] < - hysteresis:
        ckr[n] = 0
    else:
        ckr[n] = ckr[n - 1]
    count += 1
    # if count * ts >= t0:
    if f0[n - 1] < 0 < f0[n]:
        count = 0
        i += 1
        # delta_t[i - 1] = n * ts - t_k[k - 1]
        delta_t[i - 1] = i * t0 - t_k[k - 1]
        TDEV[i] = TDEV[i - 1] + delta_t[i - 1]

print("k", k)
print("t_k", t_k)
print('delta_t', delta_t)
print('TDEV', TDEV)

# plot v_t
plt.figure()
plt.plot(t[:len_simulation] / 1e-9, v_t[:len_simulation], label="V(t)")
plt.plot(t[:len_simulation] / 1e-9, f0[:len_simulation], label="f0")
plt.plot(t[:len_simulation] / 1e-9, ckr[:len_simulation], label="ckr")
plt.grid(visible=True)
plt.xticks(range(0, int(t[len_simulation] / 1e-9), 50))  # Definindo valores e step do eixo x
# plt.yticks(range(0, 11, 2)) # Definindo valores e step do eixo y
plt.legend()
# plt.plot(t[:10000]/1e-6, t_k[:10000])
plt.xlabel('Time (ns)')
plt.ylabel('Amplitude')

# plt.figure()
# plt.stem(np.arange(0, len(TDEV), 1), TDEV, linefmt='b', markerfmt='+', basefmt="-b", label="TDEV")
# plt.stem(np.arange(0, len(delta_t), 1), delta_t, linefmt='r', markerfmt='*', basefmt="-b",label="delta")
# plt.legend()
# plt.xlabel('k')
# plt.ylabel('count')
plt.show()

'''
DÚVIDAS:
Sempre utilizar a borda do sinal mais próxima ao periodo desejado? 
E quando a sáida conter duas bordas de subida dentro de um periodo T0, se pegar somente a ultima borda
terá um delta pequeno, em relação ao correto.

'''