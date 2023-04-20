import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
warnings.filterwarnings('ignore')

# Símbolo 's' para a variável de Laplace
s = sp.symbols('s')
N = 1   # relação de f/f_r
fr = 26e6   # Frequeência de referência
a = 2**-7    # alpha value
p = 2**-15    # rho value
# fr = 40e6   # Frequeência de referência
# a = 2**-6    # alpha value
# p = 2**-14    # rho value

# Função de tranferência do filtro IIR
IRR_order = 4
l = [2**-3, 2**-3, 2**-3, 2**-4]
IRR = (1 + s/fr) / (1 + (s / (l[0]*fr)))
for i in range(1, IRR_order - 1):
    IRR *= (1 + s/fr) / (1 + (s / (l[i]*fr)))
IRR *= (1 + s/fr) / (1 + (s / (l[3]*fr)))
print(IRR)

# Função de tranferência de loop aberto
Hol = (a + p * fr/s) * (fr/s)
# Hol = (p * fr**2)/s * (1 + s /(p * fr / a))/s
print(Hol)

Hol = Hol * IRR
print(Hol)
# Converter a string em uma expressão simbólica
expressao = sp.sympify(Hol)
# Obter os numeradores e denominadores como objetos simbólicos
num_Hol, den_Hol = expressao.as_numer_denom()
print("numerador em s:",num_Hol)
print("denominador em s:",den_Hol)
# Extrair coeficientes do numerador e denominador de Hs
num_Hol = sp.Poly(num_Hol, s).all_coeffs()
den_Hol = sp.Poly(den_Hol, s).all_coeffs()
print("numerador:",num_Hol)
print("denominador:",den_Hol)
# Tranforma em array para aplicar a transformada de laplace
num_array_hol = np.array(num_Hol, dtype=float)
den_array_hol = np.array(den_Hol, dtype=float)

# # Aplica a tranformada de laplace
Hol = control.tf(num_array_hol,den_array_hol) # Trarnsformada de loop aberto
print("função de transferencia em Laplace: ", Hol)
print(Hol.poles())
print(Hol.zeros())
w = np.logspace(4,9,10000)    #  List of frequencies in rad/sec to be used for frequency response ( 10^-1 até 10^3)
# Plotar bode da TF
mag,phase,f = control.bode(Hol,w,Hz=True,dB=True,deg=True, margins=True)
plt.tight_layout()
plt.show()


