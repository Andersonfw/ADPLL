import pandas as pd
import matplotlib.pyplot as plt

# Especificar o caminho do arquivo CSV
PN_witout_IRR = 'Phaseresultssimulations.csv'
PN_with_IRR = 'PN_with_IRR.csv'

# Carregar o arquivo CSV em um DataFrame do pandas
df = pd.read_csv(PN_witout_IRR, sep=';')
df1 = pd.read_csv(PN_with_IRR, sep=';')
# Extrair os dados das colunas
print(df['PN'])
Xdb_o = df['PN']
f = df['freq']

Xdb_1 = df1['PN']
f_1 = df1['freq']

ponto_especifico_f = 1e6  # Substitua pelo valor específico de frequência desejado
ponto_especifico_Xdb_o = Xdb_1[f == ponto_especifico_f].values[0]

# Plotar um gráfico de barras simples
plt.figure()
plt.semilogx(f , Xdb_o , label="Phase Noise")
plt.semilogx(f_1 , Xdb_1 , label="Phase Noise + IRR")
plt.scatter(ponto_especifico_f, ponto_especifico_Xdb_o, color='black', marker='o', label=f'{ponto_especifico_Xdb_o:.2f} dBc/Hz')
plt.grid(visible=True)
plt.legend()
plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
plt.xlabel('Freq (Hz)')
plt.ylabel('Phase Noise [dBc/Hz]')
plt.show()
