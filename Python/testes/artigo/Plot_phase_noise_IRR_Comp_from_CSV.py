import pandas as pd
import matplotlib.pyplot as plt
import os
import scienceplots as sp

ENGLISH = True
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
PN_with_IRR = os.path.join(dirr, '', 'phasenoise_PSD_plus_IRR.csv')
PN_witout_IRR = os.path.join(dirr, '', 'phasenoise_PSD_without_IRR.csv')

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN_witout_IRR, sep=';')
df1 = pd.read_csv(PN_with_IRR, sep=';')
# Extrair os dados das colunas
print(df['y'])
Xdb_o = df['y']
f = df['x']

Xdb_1 = df1['y']
f_1 = df1['x']

ponto_especifico_f = 1e6  # Substitua pelo valor específico de frequência desejado
ponto_especifico_Xdb_o = Xdb_1[f == ponto_especifico_f].values[0]

# Plotar um gráfico de barras simples
plt.style.use(['science','ieee'])
# plt.figure()
plt.figure(figsize=(3.54,2.5), dpi=600)
if ENGLISH:
    label1 = "Phase Noise"
    label2 = "Phase Noise + IIR Filter"
else:
    label1 = "Ruído de fase"
    label2 = "Ruído de fase + Filtro IIR"
plt.semilogx(f , Xdb_o , label=label1)
plt.semilogx(f_1 , Xdb_1 , label=label2)
# plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @1 MHz')
plt.grid(visible=True)
plt.legend()#fontsize=12)
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Freq. (Hz)')#, fontsize=12)
plt.ylabel(label1 +' [dBc/Hz]')#, fontsize=12)

plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\artigo\PN_IRR_COMPARE.png', bbox_inches='tight')
plt.show()