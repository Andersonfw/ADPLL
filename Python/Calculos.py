
import numpy as np
import decimal



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


'''
    Noises
'''

noise_floor = -150      # noise floor [dBc)
L_j = 10 ** (noise_floor/10)    # noise level
# f0 = 2.4e9  # desired frequency
f0 = 2045e6
t0 = 1/f0   # period of frequency

jitter = (t0 / (2*np.pi)) * np.sqrt(L_j * f0)   # Jitter noise standard deviation

print("jitter noise",jitter)

Thermal_noise = -130        # Up converted Thermal noise with deltaf frequency offset [dBc]
L_w = 10 ** (Thermal_noise/10)  # noise level
deltaf = 3.5e6  # offset frequency
wander = deltaf/f0 * np.sqrt(t0) * np.sqrt(L_w)   # Wander noise standard deviation
print("Wander noise",wander)


numero_decimal = decimal.Decimal(jitter)
# Arredonda o número para 12e-15 com uma precisão de 15 dígitos
numero_arredondado = numero_decimal.quantize(decimal.Decimal('1e-15'), rounding=decimal.ROUND_HALF_EVEN)

# Converte o número arredondado de volta para notação científica
numero_notacao_cientifica = '{:e}'.format(numero_arredondado)

# Imprime o número arredondado em notação científica
print(numero_notacao_cientifica)


L = 1e-9  # Indutor utilizado
C0 = 4.809222646671402e-12  # valor de capacitância inicial
pvt_lsb = 1.1923378339857862e-14  # valor do LSB em PVT mode
acq_lsb = 2.316699747480977e-15  # valor do LSB em acquisition mode
trk_i_lsb = 1.8511454815068396e-16  # valor do LSB em Trekking integer mode
pvt = 256 - 78 #
acq = 256 - 130 #130
trk_i = 64 - 39 #44
f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb)))
f2 = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + (trk_i-3) * trk_i_lsb)))
print("freq: ",f)
print("freq: ",f2)
print("freq diff: ",(1/f-1/f2))