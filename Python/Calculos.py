
import numpy as np

L =  1e-9
C0 = 4.809222646671402e-12
PVT_LSB = 1.1923378339857862e-14

fc = 2045e6
nb = 5
fr = 31e3
fmin = fc - fr / 2
fmax = fc + fr / 2
cmax = 1 / (L * (2 * np.pi * fmin) ** 2)
cmin = 1 / (L * (2 * np.pi * fmax) ** 2)
lsb = (cmax - cmin) / 2 ** nb
freq_lsb = fr / 2 ** nb

F0 = 1 / (2 * np.pi * np.sqrt(L * (C0)))

print(F0)

F = 1 / (2 * np.pi * np.sqrt(L * (C0 + 1 * lsb)))

print("F1",F, "  Dif", (F0 - F))