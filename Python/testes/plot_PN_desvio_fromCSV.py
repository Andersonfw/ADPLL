import pandas as pd
import matplotlib.pyplot as plt
import os

ENGLISH = True
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
PN_without_des = os.path.join(dirr, '..', 'PN_sem_desvio_1.csv')
PN_with_des = os.path.join(dirr, '..', 'PN_com_desvio_1.csv')

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN_without_des, sep=';')
df1 = pd.read_csv(PN_with_des, sep=';')
# Extrair os dados das colunas
print(df['y'])
Xdb_o = df['y']
f = df['x']

Xdb_1 = df1['y']
f_1 = df1['x']

ponto_especifico_f = 1e6  # Substitua pelo valor específico de frequência desejado
ponto_especifico_Xdb_o = Xdb_1[f == ponto_especifico_f].values[0]

# Plotar um gráfico de barras simples
plt.figure()
if ENGLISH:
    label1 = "IDEAL"
    label2 = "5% INACCURACY"
else:
    label1 = "Ruído de fase"
    label2 = "Ruído de fase + desvio"
plt.semilogx(f , Xdb_o , label=label1)
plt.semilogx(f_1 , Xdb_1 , label=label2)
# plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @1 MHz')
plt.grid(visible=True)
plt.legend()
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Freq. (Hz)')
plt.ylabel(label1 +' [dBc/Hz]')
plt.show()