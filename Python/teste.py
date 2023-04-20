import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
warnings.filterwarnings('ignore')

# Símbolo 's' para a variável de Laplace
s = sp.symbols('s')
K = 1
s_num = 100 * K
# s_den = s * (s + 36) * (s + 100)
s_den = (1 + 100*2*np.pi)

# Função de transferecia em s
Hs = s_num / s_den
print(Hs)

# Extrair coeficientes do numerador e denominador de Hs
num = sp.Poly(s_num, s).all_coeffs()
den = sp.Poly(s_den, s).all_coeffs()
print("numerador:",num)
print("denominador:",den)

# Tranforma em array para aplicar a transformada de laplace
num_array = np.array(num, dtype=float)
den_array = np.array(den, dtype=float)

# Aplica a tranformada de laplace
G = control.tf(num_array,den_array) # Trarnsformada de loop aberto
# G = control.feedback(G, 1)          # Transformada de loop fechado
print(G)

'''
Create a linear system that approximates a delay
'''
# (num,den) = control.pade(0.1,2)
# Gp = control.tf(num,den)*G
# print(Gp)
Gp = G

w = np.logspace(-1,3,100)    #  List of frequencies in rad/sec to be used for frequency response ( 10^-1 até 10^3)
# Plotar bode da TF
mag,phase,f = control.bode(Gp,w,Hz=True,dB=True,deg=True, margins=False)
plt.tight_layout()
plt.show()


# ax1, ax2 = plt.gcf().axes  # get subplot axes
# plt.sca(ax1)  # magnitude plot
# Encontre a frequência correspondente a -3dB
max_value = max(mag)
print("valor maximo = ", 20*np.log10(max_value))
max_index = np.argmax(max(mag))
aux_vector = np.argwhere(mag <= (max_value * np.sqrt(2)/2))
idx_3db = aux_vector[0]
print("index = ",idx_3db, "Valor no ponto de -3dB =", 20*np.log10(mag[idx_3db]))#, "dB. frequência =",f[idx_3db],"Hz")
# indice_x = np.argmin(np.abs(f - (max_value * np.sqrt(2)/2)))
# print(indice_x)
# f_index = int(idx_3db/np.sqrt(2))
# print(f_index)
# # Marque o ponto de -3dB no gráfico
# ax1.axvline(x=f[idx_3db], ymin=-100, ymax=100, color='red', linestyle='--')
# ax1.scatter(f[idx_3db], 20*np.log10(mag[idx_3db]), color='red', marker='o')
# # Adicionar o marcador com valor e posição específicos
# posicao_marcador = (f[idx_3db], 20*np.log10(mag[idx_3db]))
# ax1.annotate(f'Frequência {f[idx_3db]}Hz',  # Valor do marcador
#              xy=posicao_marcador,  # Posição do marcador
#              # xycoords='data',  # Coordenadas em relação aos dados do gráfico
#              xytext=(f[idx_3db - 10], mag[idx_3db - 10]),  # Posição do texto do valor do marcador
#              textcoords='offset points',  # Coordenadas do texto em relação ao marcador
#              arrowprops=dict(arrowstyle="->"))
# plt.show()