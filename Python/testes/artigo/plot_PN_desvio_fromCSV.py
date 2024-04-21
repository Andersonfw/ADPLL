import pandas as pd
import matplotlib.pyplot as plt
import os
import scienceplots as sp

ENGLISH = True
# Especificar o caminho do arquivo CSV
dirr = os.path.dirname(__file__)
PN_without_des = os.path.join(dirr, '', 'PN_sem_desvio_1.csv')
PN_with_des = os.path.join(dirr, '', 'PN_com_desvio_1.csv')
# savefig_path = os.path.join("C:\Users\ander\OneDrive\Área de Trabalho\artigo\", 'pn_plus_desvio')

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

plt.style.use(['science','ieee'])
plt.rcParams['legend.frameon'] = True  # Mostrar a moldura da legenda
plt.rcParams['legend.edgecolor'] = 'lightgray'  # Cor da borda da legenda
plt.rcParams['legend.facecolor'] = 'lightgray'  # Cor do fundo da legenda
plt.rcParams['legend.shadow'] = False  # Adicionar sombra à legenda
# plt.figure()

# plt.figure(figsize=(1,1), dpi=600)
# sp.setup_figure(size=(1, 1), dpi=600)
if ENGLISH:
    plt.figure(figsize=(3.54,2.5), dpi=600)
    label1 = "Ideal"
    label2 = "5\% Variation"
    size_font = 8
else:
    plt.figure(figsize=(6,4), dpi=600)
    label1 = "Ideal"
    label2 = "5\% Variação"
    size_font = 12
plt.semilogx(f , Xdb_o , label=label1)
plt.semilogx(f_1 , Xdb_1 , label=label2)
# plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz @1 MHz')
plt.grid(visible=True)
# plt.legend(fontsize=12)
plt.legend(facecolor='white', framealpha=1, fontsize=size_font)
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70],fontsize=size_font)
plt.xticks(fontsize=size_font)
plt.xlabel('Freq. (Hz)', fontsize=size_font)
if ENGLISH:
    plt.ylabel('Phase noise' +' [dBc/Hz]', fontsize=size_font)
else:
    plt.ylabel('Ruído de fase' +' [dBc/Hz]', fontsize=size_font)
if ENGLISH:
    plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\artigo\pn_plus_desvio.eps', bbox_inches='tight',format='eps')
else:
    plt.savefig(r'C:\Users\ander\OneDrive\Área de Trabalho\Imagens\pn_desvio_compare1.eps', bbox_inches='tight',format='eps')
plt.show()
