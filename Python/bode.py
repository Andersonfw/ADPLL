import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

R = 100e3
C = 1e-3
# Defina a função de transferência em s
num = [1]  # Numerador
den = [1, R*C]  # Denominador
sys = signal.TransferFunction(num, den)

# Calcule a resposta em frequência
f, mag, phase = signal.bode(sys)

# Plote o mapa de Bode
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.semilogx(f, mag)  # Plota o ganho em escala logarítmica
ax1.set_ylabel('Ganho (dB)')
ax2.semilogx(f, phase)  # Plota a fase em escala logarítmica
ax2.set_xlabel('Frequência (Hz)')
ax2.set_ylabel('Fase (graus)')

# Encontre a frequência correspondente a -3dB
idx_3db = np.argmax(mag < max(mag - 3))

# Marque o ponto de -3dB no gráfico
ax1.axvline(x=f[idx_3db], ymin=-100, ymax=100, color='red', linestyle='--')
ax1.scatter(f[idx_3db], mag[idx_3db], color='red', marker='o')

# Adicionar o marcador com valor e posição específicos
posicao_marcador = (f[idx_3db], mag[idx_3db])
ax1.annotate(f'Frequeância {f[idx_3db]}Hz',  # Valor do marcador
             xy=posicao_marcador,  # Posição do marcador
             # xycoords='data',  # Coordenadas em relação aos dados do gráfico
             xytext=(f[idx_3db - 10], mag[idx_3db - 10]),  # Posição do texto do valor do marcador
             textcoords='offset points',  # Coordenadas do texto em relação ao marcador
             arrowprops=dict(facecolor='black',shrink=0.05))  # Estilo de seta do marcador
plt.show()
