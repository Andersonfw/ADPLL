import matplotlib.pyplot as plt
import control
import numpy as np


def bode(name, sys, omega, margins=False):
    mag, phase, omega = control.bode(sys, omega, Plot=False)
    mag_dB = 20 * np.log10(mag)
    if margins:
        gm, pm, sm, wg, wp, ws = control.stability_margins(sys)

    plt.subplot(211)
    if margins:
        plt.hlines(0, omega[0], omega[-1], linestyle='--')
        plt.vlines([wp, wg], np.min(mag_dB), np.max(mag_dB), linestyle='--')
    plt.semilogx(omega, mag_dB)
    plt.xlabel('rad')
    plt.ylabel('dB')
    plt.grid()

    if margins:
        plt.title(
            name +
            ' bode pm: {:0.2f} deg @{:0.2f} rad/s gm: {:0.2f} @{:0.2f} rad/s'.
            format(pm, wp, gm, wg))
    else:
        plt.title(name + ' bode')
    plt.subplot(212)
    phase_deg = np.rad2deg(phase)
    plt.semilogx(omega, phase_deg)
    if margins:
        plt.vlines([wp, wg],
                   np.min(phase_deg),
                   np.max(phase_deg),
                   linestyle='--')
        plt.hlines(-180, omega[0], omega[-1], linestyle='--')
    plt.grid()

# Definir as funções de transferência H1(s) e H2(s)
H1 = control.TransferFunction([1], [1, 2, 1])
H2 = control.TransferFunction([1], [1, 3, 2])

bode("teste",H1, (1e-3, 1e3),margins=True)
# Criar uma figura com dois subplots (1 linha, 2 colunas)
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
#
# # Plotar a magnitude e fase de H1(s) na primeira subfigura
# control.bode(H1, omega_limits=(1e-3, 1e3), dB=True, deg=True, kargs="color=red")
# ax1.set_title('H1(s)')
#
# # Plotar a magnitude e fase de H2(s) na segunda subfigura
# # control.bode_plot(H2, omega_limits=(1e-3, 1e3), dB=True, deg=True, axmag=ax2, axphase=ax2)
# # ax2.set_title('H2(s)')
#
# # Adicionar rótulos aos eixos
# ax1.set_xlabel('Frequência (rad/s)')
# ax1.set_ylabel('Magnitude (dB)')
# ax2.set_xlabel('Frequência (rad/s)')
# ax2.set_ylabel('Fase (graus)')
#
# # Exibir a figura
plt.show()
