import pandas as pd
import matplotlib.pyplot as plt
import os
import scienceplots


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
# plt.figure()
plt.style.use(['science','ieee'])
plt.figure(figsize=(5.8,4), dpi=600)
plt.rcParams['legend.frameon'] = True  # Mostrar a moldura da legenda
plt.rcParams['legend.edgecolor'] = 'lightgray'  # Cor da borda da legenda
plt.rcParams['legend.facecolor'] = 'lightgray'  # Cor do fundo da legenda
plt.rcParams['legend.shadow'] = False  # Adicionar sombra à legenda

if ENGLISH:
    label1 = "Reference Clock Cycles"
    label2 = 'Frequency out of DCO (Hz)'
else:
    label1 = "$K_{P_{TRK}} = 2^{-5}$"
    label2 = '$K_{P_{TRK}} = 2^{-3}$'
plt.plot(f , Xdb_o , label=label1)
plt.plot(f_1 , Xdb_1 , label=label2)

plt.annotate('Modo\nPVT', xy=(6, 40), xytext=(36, 20),
             arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=2), fontsize=12)

plt.annotate('Modo\nACQ', xy=(38, 67), xytext=(68, 87),
             arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=2), fontsize=12)

plt.annotate('Modo\nTRQ', xy=(164, 47), xytext=(194, 27),
             arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=2), fontsize=12)
# plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @1 MHz')
plt.grid(visible=True)
plt.legend(facecolor='white', framealpha=1,loc='lower right',fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
# plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Ciclos de clock de referência',fontsize=12)
plt.ylabel("OTW",fontsize=12)
plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\Imagens\OTW_plot_style_both_.eps', bbox_inches='tight',format='eps')
plt.show()
