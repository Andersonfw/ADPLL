"""
Created on julho 10 22:00:20 2023

@author: Ânderson Felipe Weschenfelder
"""
import matplotlib.pyplot as plt
from sympy import symbols, solve, Rational, simplify
import sympy as sp
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
# Defina as variáveis simbólicas
# r1, r2, R, Cgs2, CL, gm2, Rs, cgs1 , s = symbols('r1 r2 R Cgs2 CL gm2 Rs cgs1 s')
s, Rs, cgs1= symbols('s Rs cgs1')

r1 = 20e3
r2 = r1
Cgs2 = 100e-15
R = 100e3
CL = 10e-12
gm2 = 1e-3
Rs = 100e3
cgs1 = Cgs2

L =(Cgs2 / gm2)*(R - 1/gm2)*s

# Defina a função de transferência
# H = (R * C * s**2 + 5*R*C*s + 1)
num = r1*r2* Cgs2*(gm2*R -1)*s
den = r1*r2*(CL * Cgs2 * s**2*(R*gm2 -1) + gm2) + Cgs2*(gm2*R -1)*s*(r1+r2)
# H = simplify(num / den)
H = simplify( (1 /(1/L + 1/r1 + 1/r2 + CL*s))*(1/(Rs*cgs1*s +1)))
print("H",H)
# H = simplify(H * 1/(Rs*cgs1*s +1))
# print(H)

H2 = s / (s**2*CL + s/r1 + s/r2 + gm2**2/(Cgs2*(R*gm2 -1)))
print ("H2",simplify(H2))
# H = 1 / (Rational(1, R * C) * s + 1)
# Calcule os polos
# polos = solve(H_simplified., s)

# Calcule os zeros
# zeros = solve(H.as_numer_denom()[0], s)
expressao = sp.sympify(H)
# Obter os numeradores e denominadores como objetos simbólicos
num_Hol, den_Hol = expressao.as_numer_denom()
# print("numerador em s:", num_Hol)
# print("denominador em s:", den_Hol)
# Extrair coeficientes do numerador e denominador de Hs
num_Hol = sp.Poly(num_Hol, s).all_coeffs()
den_Hol = sp.Poly(den_Hol, s).all_coeffs()
num = np.array(num_Hol, dtype=float)
den = np.array(den_Hol, dtype=float)
sys = signal.TransferFunction(num, den)

# Frequências para avaliar a resposta em frequência
frequencies = np.logspace(1, 15, num=1000)
w, mag, phase = signal.bode(sys, frequencies)

# Plotar o diagrama de magnitude
plt.figure()
plt.semilogx(w, mag)
plt.xlabel('Frequência (rad/s)')
plt.ylabel('Magnitude (dB)')
plt.title('Diagrama de Bode - Magnitude')
plt.grid(True)
# Imprima os polos e zeros
print("Polos:", num)
print("Zeros:", den)
print("\r\n")
polos = solve(H.as_numer_denom()[1], s)
zeros = solve(H.as_numer_denom()[0], s)

# Imprima os polos e zeros
print("Polos:", abs(polos[0]))
print("Polos:", polos)
print("Zeros:", zeros)

H3 = s / (s**2 +  gm2/(CL*Cgs2*(R)))
polos = solve(H3.as_numer_denom()[1], s)
zeros = solve(H3.as_numer_denom()[0], s)
print("Polos:", abs(polos[0]))
print("Zeros:", zeros)


hR = []
# r1 = 5e3
# r2 = r1
R = 10e3
for i in range(10):
    eq = sp.sympify(1 /(1/L + 1/r1 + 1/r2 + CL*s))
    L = (Cgs2 / gm2) * (R - 1 / gm2) * s
    hR.append(abs(solve(eq.as_numer_denom()[1], s)[0]))
    # r1 = r1 + 15e3
    # r2 = r1
    R = R + 20e3

hcl = []
R = 100e3
CL = 10e-12
for i in range(10):
    eq = sp.sympify(1 /(1/L + 1/r1 + 1/r2 + CL*s))
    L = (Cgs2 / gm2) * (R - 1 / gm2) * s
    hcl.append(abs(solve(eq.as_numer_denom()[1], s)[0]))
    # r1 = r1 + 15e3
    # r2 = r1
    CL = CL + 20e-12
# print(hcgs)
# plt.Figure
# plt.plot(range(0,10,1), hR, 'b')
# plt.plot(range(0,10,1), hcl,'g')
plt.show()