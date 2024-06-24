import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
from scipy import signal
import pandas as pd
import decimal

noise_floor = -150  # -150  # noise floor [dBc)
L_j = 10 ** (noise_floor / 10)  # noise level
f_desired = 2.4e9 * 2#* DIVISION_OUTPUT # F0  # desired frequency
t_required = 1 / f_desired  # period of frequency
Thermal_noise = -105   # 6dB acima do desjado para dobro de freq.  # Up converted Thermal noise with deltaf frequency offset [dBc]
L_w = 10 ** (Thermal_noise / 10)  # noise level
deltaf = 0.5e6  # offset frequency

j_noise = (t_required / (2 * np.pi)) * np.sqrt(L_j * f_desired)  # Jitter noise standard deviation
W_noise = deltaf / f_desired * np.sqrt(t_required) * np.sqrt(L_w)  # Wander noise standard deviation (including the 1/f noise)
#W_noise = 1/np.sqrt(2) * W_noise
# Converte o número em um decimal
j_decimal = decimal.Decimal(j_noise)
w_decimal = decimal.Decimal(W_noise)
# Arredonda o número com uma precisão de 15 dígitos
j_decimal = j_decimal.quantize(decimal.Decimal('1e-15') , rounding=decimal.ROUND_HALF_EVEN)
w_decimal = w_decimal.quantize(decimal.Decimal('1e-15') , rounding=decimal.ROUND_HALF_EVEN)
# Converte o número arredondado de volta para notação científica
Jt_noise = float('{:e}'.format(j_decimal))
Wt_noise = float('{:e}'.format(w_decimal))
print("jitter noise" , Jt_noise)
print("Wander noise" , Wt_noise)