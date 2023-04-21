import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings

warnings.filterwarnings('ignore')

'''
            Extrai os parêmetros do numerador e denominador de uma função de transferência em S
            Parameters
            ----------
            H : Função de transferência no dominio de frequência S
'''


def extract_parameters(H):
    # Converter a string em uma expressão simbólica
    expressao = sp.sympify(H)
    # Obter os numeradores e denominadores como objetos simbólicos
    num_Hol, den_Hol = expressao.as_numer_denom()
    # print("numerador em s:", num_Hol)
    # print("denominador em s:", den_Hol)
    # Extrair coeficientes do numerador e denominador de Hs
    num_Hol = sp.Poly(num_Hol, s).all_coeffs()
    den_Hol = sp.Poly(den_Hol, s).all_coeffs()
    # print("numerador:", num_Hol)
    # print("denominador:", den_Hol)
    # Tranforma em array para aplicar a transformada de laplace
    num = np.array(num_Hol, dtype=float)
    den = np.array(den_Hol, dtype=float)
    return num, den


'''
            Calcula o filtro IRR de acordo com os coeficientes recebidos
            Parameters
            ----------
            lp : array de coeficientes
'''


def IRR_filter(lp):
    IRR_order = len(lp)
    IRR = 1
    for i in range(IRR_order - 1):
        IRR *= (1 + s / fr) / (1 + (s / (lp[i] * fr)))
    IRR *= (1 + s / fr) / (1 + (s / (lp[IRR_order - 1] * fr)))

    return IRR


'''
            Plota o gráfico de Bode na função tranferência no dominío S 
            Parameters
            ----------
            H : (opcional) Função de transferência no dominío de frequência S
            numerador : (opcional) array dos parâmetros do numerador da Função de transferência sem "S"
            denominador : (opcional) array dos parâmetros do denominador da Função de transferência sem "S"
            irr : (opcional) array de coeficientes do filtro IRR
            w : (opcional) Kista de frequências para ser usada pela resposta em frequência
            margem : (opcional) Plotar a margem e ganho de fase
            prints : (opcional) Printar informações
'''


def h_to_bode(H=None, numerador=None, denominador=None, irr=None, freq=None, margem=False, prints=False):
    if irr:
        IRR = IRR_filter(irr)
    else:
        IRR = 1
    if H:
        num_irr, den_irr = extract_parameters(IRR)
        IRR = control.tf(num_irr, den_irr)
        num_H, den_H = extract_parameters(H)
        tf = control.tf(num_H, den_H)
    else:
        num_irr, den_irr = extract_parameters(IRR)
        IRR = control.tf(num_irr, den_irr)
        tf = control.tf(numerador, denominador)
    if prints:
        print("Função de transfêrencia sem filtro IRR: ", tf)
    tf = tf * IRR
    # tf = (tf /(1 + tf))
    # tf = (1 / (1 + tf))
    # tf = control.feedback(tf, 1)
    if prints:
        print("função de transferencia em Laplace com filtro IRR: ", tf)
    if freq is not None:
        mag, phase, f = control.bode(tf, omega=freq, Hz=True, dB=True, deg=True, margins=margem)
    else:
        mag, phase, f = control.bode(tf, Hz=True, dB=True, deg=True, margins=margem)
    return mag, phase, f, tf


'''
        DEFINIÇÕES
'''
# Símbolo 's' para a variável de Laplace
s = sp.symbols('s')
N = 69.23  # relação de f/f_r

# Livro Bogdan PG 137/152
# fr = 26e6  # Frequeência de referência
# a = 2 ** -7  # alpha value
# p = 2 ** -15  # rho value
# #   Coeficiêntes do filtro IIR
# l = [2**-3, 2**-3, 2**-3, 2**-4]

# Aassignment 4
fr = 40e6  # Frequeência de referência
a = 2 ** -6  # alpha value
p = 2 ** -14  # rho value
# Coeficiêntes do filtro IIR
l = [2**-2, 2**-3, 2**-2, 2**-3]

# Função de tranferência de loop aberto
Hol = (a + p * fr / s) * (fr / s)


if __name__ == "__main__":
    w = np.logspace(3, 10, 10000)  # List of frequencies in rad/sec to be used for frequency response ( 10^-1 até 10^3)
    # # Plotar bode da TF
    # mag, phase, f, Hol = h_to_bode(H=Hol, freq=w, irr=l, margem=True, prints=True)
    mag, phase, f, Hol = h_to_bode(H=Hol, freq=w, irr=None, margem=True, prints=True)
    ax1, ax2 = plt.gcf().axes  # get subplot axes
    ax1.set_title('Magnitude Open Loop Hol(s)')
    ax2.set_title('Phase Open Loop Hol(s)')

    plt.figure()
    # Função de tranferência de loop fechado para referência
    Hol_ref = N * (Hol / (1 + Hol))
    control.bode(Hol_ref, omega=w, Hz=True, dB=True, deg=True, margins=True)
    # mag, phase, f, Hol = h_to_bode(H=Hol_ref, freq=w, prints=True)
    ax1, ax2 = plt.gcf().axes  # get subplot axes
    ax1.set_title('Magnitude Closed Loop Hcl reference(s)')
    ax2.set_title('Phase Closed Loop Hcl reference (s)')

    plt.figure()
    # Função de tranferência de loop fechado para o TDC
    Hol_TDC = (Hol / (1 + Hol))
    control.bode(Hol_TDC, omega=w, Hz=True, dB=True, deg=True, margins=True)
    ax1, ax2 = plt.gcf().axes  # get subplot axes
    ax1.set_title('Magnitude Closed Loop Hcl TDC(s)')
    ax2.set_title('Phase Closed Loop Hcl TDC (s)')

    plt.figure()
    # Função de tranferência de loop fechado para o DCO
    Hol_DCO = ( 1 / (1 + Hol))
    control.bode(Hol_DCO, omega=w, Hz=True, dB=True, deg=True, margins=True)
    ax1, ax2 = plt.gcf().axes  # get subplot axes
    ax1.set_title('Magnitude Closed Loop Hcl DCO(s)')
    ax2.set_title('Phase Closed Loop Hcl DCO (s)')

    plt.tight_layout()
    plt.show()
