import pandas as pd
import matplotlib.pyplot as plt
import os

ENGLISH = False
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
PN_with_IRR = os.path.join(dirr, '..', 'OTW_ACQ.csv')
PN_witout_IRR = os.path.join(dirr, '..', 'OTW.csv')

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN_witout_IRR, sep=';')
df1 = pd.read_csv(PN_with_IRR, sep=';')
# Extrair os dados das colunas
print(df['y'])
Xdb_o = df['x']
f = df['y']

Xdb_1 = df1['x']
f_1 = df1['y']

ponto_especifico_f = 1e6  # Substitua pelo valor específico de frequência desejado
# ponto_especifico_Xdb_o = Xdb_1[f == ponto_especifico_f].values[0]

# Plotar um gráfico de barras simples
plt.figure()
if ENGLISH:
    label1 = "Reference Clock Cycles"
    label2 = 'Frequency out of DCO (Hz)'
else:
    label1 = "$K_{P_{TRK}} = 2^{-5}$"
    label2 = '$K_{P_{TRK}} = 2^{-3}$'
plt.plot(f , Xdb_o , label=label1)
plt.plot(f_1 , Xdb_1 , label=label2)
# plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @1 MHz')
plt.grid(visible=True)
plt.legend(loc='lower right')
# plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Ciclos de clock de referência')
plt.ylabel("OTW")
plt.show()
