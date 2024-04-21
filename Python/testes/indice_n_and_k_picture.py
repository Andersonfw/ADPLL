import numpy as np
import matplotlib.pyplot as plt
import scienceplots as sp

FCW = 4.5
# Gerar os dados para a curva A
x_a = np.arange(0, 3, 1/FCW) #np.arange(0, 21)  # Variando de 0 a 20
y_a = np.ones_like(x_a)  # Valor de y é inicialmente definido como 1 para a curva A

# Definir y_a como 0 quando x_a é um número fracionário
""" for i in range(len(x_a)):
    if not x_a[i].is_integer():
        y_a[i] = 0 """

# Gerar os dados para a curva B
x_b = np.arange(1, 3, 1) #np.arange(5.5, 21, 5.5)  # Variando de 0 a 20 em incrementos de 5.5
y_b = 1.5*np.ones_like(x_b)  # Valor de y é sempre 1 para a curva B

# Plotar o gráfico
plt.style.use(['science','ieee'])
# plt.figure()
plt.figure(figsize=(6.5,3.5), dpi=600)
plt.stem(x_a, y_a, label='CKV índice n',linefmt='black', markerfmt='.', basefmt="black")  # Linha sólida para curva A
plt.stem(x_b, y_b, label='$f_{ref}$ índice k',linefmt='red', markerfmt='.', basefmt="red")
#plt.plot(x_b, y_b, label='Curva B', color='red',  marker='|')  # Linha tracejada para curva B

# Adicionar flechas e rótulos entre dois pontos específicos
plt.annotate('$T_{Out}$', xy=(0.42, 0.85), xytext=(0.7, 0.83),
             arrowprops=dict(facecolor='black', arrowstyle='<|-|>', lw=0.5), fontsize=12)

plt.annotate('$T_{ref}$', xy=(0.98, 1.35), xytext=(2.026, 1.33),
             arrowprops=dict(facecolor='black', arrowstyle='<|-|>', lw=0.5), fontsize=12)

# plt.annotate('$T_{Out}$', xy=(2/FCW-0.0, 1), xytext=(3/FCW+0.01, 1-0.03),
#              arrowprops=dict(facecolor='blue', arrowstyle='<|-|>', lw=2.5), fontsize=15)

# plt.annotate('$T_{ref}$', xy=(1, 2), xytext=(2.013, 1.975),
#              arrowprops=dict(facecolor='blue', arrowstyle='<|-|>', lw=2.5), fontsize=15)

# Configurar os limites dos eixos
plt.xlim(0, 3)
plt.yticks([0, 0.5, 1, 1.5 ])
# plt.ylim(0, 1.5)

# Adicionar rótulos e título
#plt.xlabel('Eixo X')
#plt.ylabel('Eixo Y')
tick_positions = np.arange(0, 3, 1/FCW)
plt.tick_params(axis='both', which='major')#, labelsize=12)  # Defina o tamanho da fonte para os números nos eixos x e y
tick_labels = [f'CKV({i})' for i in range(len(tick_positions))]
# tick_labels = [f'{i}' for i in range(len(tick_positions))]
# tick_labels = [0,4,8,16,20]
# plt.xticks(np.arange(0, 3, 1/FCW), tick_labels)
plt.xticks(np.arange(0, 3, 1/FCW), tick_labels,  rotation=45, fontsize=10)
plt.yticks(fontsize=12)
#plt.title('Curvas A e B')

# Adicionar legenda
plt.legend(loc='upper left',fontsize=12)
# plt.xlabel('CKV')#,fontsize=12)
# plt.ylabel(fontsize=12)
# plt.ylabel(label1 + ' [dBc/Hz]')
# Mostrar o gráfico
# plt.grid(True)
plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\Imagens\ckv_tref_explain.eps', bbox_inches='tight',format='eps')
plt.show()
