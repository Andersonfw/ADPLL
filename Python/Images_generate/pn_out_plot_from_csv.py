import pandas as pd
import matplotlib.pyplot as plt
import os
import scienceplots as sp
import numpy as np

ENGLISH = False
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
# PN = os.path.join(dirr, '', '..\pn_2_4_final.csv')
PN = os.path.join(dirr, '', '..\PN_noise_model_valid.csv')


# savefig_path = os.path.join("C:\Users\ander\OneDrive\Área de Trabalho\artigo\", 'pn_plus_desvio')

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN, sep=';')
# Extrair os dados das colunas
print(df['y'])
Xdb_o = df['y']
f = df['x']

marker = 525e3  # Substitua pelo valor específico de frequência desejado

# marker = 1e6  # Substitua pelo valor específico de frequência desejado
indice = np.where(f == marker)[0][0]
marker_dB = Xdb_o[indice]
# Plotar um gráfico de barras simples
plt.figure(figsize=(6,4), dpi=600)
plt.style.use(['science','ieee'])
plt.rcParams['legend.frameon'] = True  # Mostrar a moldura da legenda
plt.rcParams['legend.edgecolor'] = 'lightgray'  # Cor da borda da legenda
plt.rcParams['legend.facecolor'] = 'lightgray'  # Cor do fundo da legenda
plt.rcParams['legend.shadow'] = False  # Adicionar sombra à legenda
# plt.figure()

# plt.figure(figsize=(1,1), dpi=600)
# sp.setup_figure(size=(1, 1), dpi=600)
if ENGLISH:
    label1 = "Phase Noise"
else:
    label1 = "Ruído de fase"
plt.semilogx(f , Xdb_o , label=label1)
# plt.scatter(marker, marker_dB, color='black', marker='o', label=f'{marker_dB:.2f} dBc/Hz  @1 MHz')
plt.scatter(marker, marker_dB, color='black', marker='o', label=f'{marker_dB:.2f} dBc/Hz  @500 kHz')
plt.grid(visible=True)
# plt.legend(fontsize=12)
plt.legend(facecolor='white', framealpha=1, fontsize=12)
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70],fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel('Freq. (Hz)', fontsize=12)
if ENGLISH:
    plt.ylabel('Phase noise' +' [dBc/Hz]', fontsize=12)
else:
    plt.ylabel('Ruído de fase' +' [dBc/Hz]', fontsize=12)
if ENGLISH:
    plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\artigo\pn_plus_desvio.png', bbox_inches='tight')
else:
    plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\Imagens\PN_ruidos_sem_TDC.eps', bbox_inches='tight',format='eps')
plt.show()
