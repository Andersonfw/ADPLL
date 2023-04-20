import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import control
R = 100e3
C = 1e-3
# Defina a função de transferência em s
num = [1]  # Numerador
den = [1, R*C]  # Denominador
den = [1536, 576, -86.4, 4.6, -0.1]
num = [7680, 4800, 1.205e4, 1465, 72.5, 1.5]
sys = signal.TransferFunction(num, den)

print(sys)


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
             arrowprops=dict(arrowstyle="->"))

# find the cross-over frequency and gain at cross-over
wc = np.interp(-180.0,np.flipud(phase),np.flipud(f))
Kcu = np.interp(wc,f,mag)

print('Crossover freq = ', wc, ' HZ')
print('Gain at crossover = ', Kcu)

# ax1,ax2 = plt.gcf().axes     # get subplot axes
#
# plt.sca(ax1)                 # magnitude plot
# plt.plot(plt.xlim(),[Kcu,Kcu],'r--')
# plt.plot([wc,wc],plt.ylim(),'r--')
# plt.title("Gain at Crossover = {0:.3g}".format(Kcu))
#
# plt.sca(ax2)                 # phase plot
# plt.plot(plt.xlim(),[-180,-180],'r--')
# plt.plot([wc,wc],plt.ylim(),'r--')
# plt.title("Crossover Frequency = {0:.3g} rad/sec".format(wc))

plt.show()