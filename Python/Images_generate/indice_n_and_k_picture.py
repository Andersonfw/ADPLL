import numpy as np
import matplotlib.pyplot as plt

FCW = 4.5
# Gerar os dados para a curva A
x_a = np.arange(0, 3, 1/FCW) #np.arange(0, 21)  # Variando de 0 a 20
y_a = np.ones_like(x_a)  # Valor de y é inicialmente definido como 1 para a curva A

# Definir y_a como 0 quando x_a é um número fracionário
""" for i in range(len(x_a)):
    if not x_a[i].is_integer():
        y_a[i] = 0 """

# Gerar os dados para a curva B
x_b = np.arange(1, 4, 1) #np.arange(5.5, 21, 5.5)  # Variando de 0 a 20 em incrementos de 5.5
y_b = 2*np.ones_like(x_b)  # Valor de y é sempre 1 para a curva B

# Plotar o gráfico
plt.figure()
plt.stem(x_a, y_a, label='CKV índice n',linefmt='b', markerfmt='o', basefmt="-b")  # Linha sólida para curva A
plt.stem(x_b, y_b, label='$f_{ref}$ índice k',linefmt='r', markerfmt='o', basefmt="-r")
#plt.plot(x_b, y_b, label='Curva B', color='red',  marker='|')  # Linha tracejada para curva B

# Adicionar flechas e rótulos entre dois pontos específicos
plt.annotate('$T_{Out}$', xy=(2/FCW-0.0, 1), xytext=(3/FCW+0.01, 1-0.03),
             arrowprops=dict(facecolor='blue', arrowstyle='<|-|>', lw=2.5), fontsize=15)

plt.annotate('$T_{ref}$', xy=(1, 2), xytext=(2.013, 1.975),
             arrowprops=dict(facecolor='blue', arrowstyle='<|-|>', lw=2.5), fontsize=15)

# Configurar os limites dos eixos
plt.xlim(0, 3)
plt.ylim(0, 2.5)

# Adicionar rótulos e título
#plt.xlabel('Eixo X')
#plt.ylabel('Eixo Y')
tick_positions = np.arange(0, 3, 1/FCW)
tick_labels = [f'CKV({i})' for i in range(len(tick_positions))]
plt.xticks(np.arange(0, 3, 1/FCW), tick_labels,  rotation=45)
#plt.title('Curvas A e B')

# Adicionar legenda
plt.legend()

# Mostrar o gráfico
plt.grid(True)
plt.show()
