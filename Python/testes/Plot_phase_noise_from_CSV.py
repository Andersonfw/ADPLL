import pandas as pd
import matplotlib.pyplot as plt
import os

ENGLISH = False
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
PN_with_IRR = os.path.join(dirr, '..', 'phasenoise_PSD.csv')
PN_witout_IRR = os.path.join(dirr, '..', 'phasenoise_PSD_without_IRR.csv')

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN_witout_IRR, sep=';')
df1 = pd.read_csv(PN_with_IRR, sep=';')
# Extrair os dados das colunas
print(df['PN'])
Xdb_o = df['PN']
f = df['freq']

Xdb_1 = df1['PN']
f_1 = df1['freq']

ponto_especifico_f = 500e3  # Substitua pelo valor específico de frequência desejado
ponto_especifico_Xdb_o = Xdb_1[f == ponto_especifico_f].values[0]

# Plotar um gráfico de barras simples
plt.figure()
if ENGLISH:
    label1 = "Phase Noise"
    label2 = "Phase Noise + IRR Filter"
else:
    label1 = "Ruído de fase"
    label2 = "Ruído de fase + Filtro IRR"
plt.semilogx(f , Xdb_o , label=label1)
plt.semilogx(f_1 , Xdb_1 , label=label2)
plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @500 kHz')
plt.grid(visible=True)
plt.legend()
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Freq. (Hz)')
plt.ylabel(label1 +' [dBc/Hz]')
plt.show()
