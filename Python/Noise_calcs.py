"""
Created on setembro 12 17:46:23 2023

@author: Ânderson Felipe Weschenfelder
"""
import numpy as np


DELTA_F = 1e6
TDC_RES = 12e-12
FRES_DCO = 31.25e3
FREF = 26e6
TREF = 1/FREF
FDCO = 4.8e9/2
TDCO = int(1/FDCO / TDC_RES) * TDC_RES



noise_DCO_time = (1/12) * ((FRES_DCO/DELTA_F)**2) * TREF * (np.sinc(DELTA_F/FREF)**2)
noise_DCO_DB = 10*np.log10(noise_DCO_time)

noise_TDC_time = ( ( (2*np.pi)**2) / 12) * (TDC_RES/TDCO)**2 * TREF
noise_TDC_DB = 10*np.log10(noise_TDC_time)


print("TDC   ---- time:", noise_TDC_time, "  ", noise_TDC_DB, "DB")
print("DCO   ---- time:", noise_DCO_time, "  ", noise_DCO_DB, "DB")

import decimal

F_DESIRED = 2.4e9
DIVISION_OUTPUT = 2


noise_floor = -150  # -150  # noise floor [dBc)
L_j = 10 ** (noise_floor / 10)  # noise level
f_desired = F_DESIRED * DIVISION_OUTPUT # F0  # desired frequency
t_required = 1 / f_desired  # period of frequency
Thermal_noise = -105  # 6dB acima do desjado para dobro de freq.  # Up converted Thermal noise with deltaf frequency offset [dBc]
L_w = 10 ** (Thermal_noise / 10)  # noise level
deltaf = 0.5e6  # offset frequency

j_noise = (t_required / (2 * np.pi)) * np.sqrt(L_j * f_desired)  # Jitter noise standard deviation
W_noise = deltaf / f_desired * np.sqrt(t_required) * np.sqrt(L_w)  # Wander noise standard deviation (including the 1/f noise)
#W_noise = 1 / np.sqrt(DIVISION_OUTPUT) * W_noise
# Converte o número em um decimal
#j_noise *= np.sqrt(DIVISION_OUTPUT)
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
# Wt_noise = 12e-15  # Wander noise time
# Jt_noise = 111e-15  # jitter noise time