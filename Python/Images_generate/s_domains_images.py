import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
from scipy import signal

warnings.filterwarnings('ignore')

def fun_calc_psd(x , fs=1 , rbw=100e3 , fstep=None):
    '''
    Calculate power spectral density

    INPUT arguments
    x     :  data vector
    fs    :  sample rate [Hz]
    rbw   :  resolution bandwidth [Hz]
    fstep :  FFT frequency bin separation [Hz]
    OUTPUT
    XdB	: spectrum of x [dB]
    f	: frequency vector [Hz]
    '''

    if fstep is None:
        fstep = rbw / 1.62
    len_x = len(x)
    nwin = round(fs * 1.62 / rbw)
    nfft = round(fs / fstep)
    if nwin > len_x:
        nwin = len_x
        rbw = fs * 1.62 / nwin
    num_segments = 8
    # nwin = math.floor(len(x) / num_segments)
    fftstr = (f'len(x)={len_x:.2f}, rbw={rbw / 1e3:.2f}kHz, fstep={fstep / 1e3:.2f}kHz, nfft={nfft:d}, nwin={nwin:d}')
    print(f'Calculating the PSD: {fftstr} ...')
    f , X = signal.welch(x , fs=fs , window=signal.windows.blackman(nwin) , nperseg=nwin , nfft=nfft ,scaling='density')
    X *= (np.sinc(f / fs)) ** 2  # correct for ZOH
    XdB = 10 * np.log10(X)
    XdB_sig = np.max(XdB)
    print(f'Signal PSD peak = {XdB_sig:.2f} dB, 10log(rbw) = {10 * np.log10(rbw):.1f}')
    return XdB , f


'''
            Exibir duas funções de transferência em uma mesma figura
            Parameters
            ----------
            name1 : Nome da legenda da função 1
            name2 : Nome da legenda da função 2
            sys1 : Função de transferênci da função 1
            sys1 : Função de transferênci da função 2
            w :   Lista de frequências para ser usada pela resposta em frequência
            margem : (opcional) Plotar a margem e ganho de fase
            offset: (opcional) Offset de magnitude
'''
def format_plot(name1, name2, sys1, sys2,  omega, margins=False, offset=0):
    mag1, phase1, omega1 = control.bode(sys1, omega, Hz=False, Plot=False)
    mag2, phase2, omega2 = control.bode(sys2, omega, Hz=False, Plot=False)

    mag_dB1 = 20 * np.log10(mag1) + offset
    mag_dB2 = 20 * np.log10(mag2) + offset
    if margins:
        gm1, pm1, sm1, wg1, wp1, ws1 = control.stability_margins(sys1)
        gm2, pm2, sm2, wg2, wp2, ws2 = control.stability_margins(sys2)
        wp1 = wp1 / (2 * np.pi)
        wg1 = wg1 / (2 * np.pi)
        wp2 = wp2 / (2 * np.pi)
        wg2 = wg2 / (2 * np.pi)

    omega1 = omega1 / (2 * np.pi)
    omega2 = omega2 / (2 * np.pi)

    #mag_dB1 -= 20 * np.log10(1/omega1) - 10 * np.log10(omega1[0] / omega1)

    plt.subplot(211)
    """     if margins:
        plt.hlines(0, omega1[0], omega1[-1], linestyle='--')
        plt.vlines([wp1, wg1], np.min(mag_dB1), np.max(mag_dB1), linestyle='--') 
    """

    plt.semilogx(omega1, mag_dB1, '-r',  label="{}".format(name1))
    plt.xlabel('Frequência (Hz)')
    step = 10
    # plt.yticks(np.arange(-180, -90, step))
    plt.ylabel('PN (dBc/Hz)')
    plt.legend()
    plt.grid()

    """     if margins:
            plt.title(
                ' bode pm: {:0.2f} deg @{:0.2f} Hz gm: {:0.2f} dB @{:0.2f} Hz'.
                format(pm1, wp1,  20 * np.log10(gm1), wg1))
        else: 
    """
   # plt.title(name1 + ' bode')
    """     plt.subplot(223)
        phase_deg = np.rad2deg(phase1)
        plt.semilogx(omega1, phase_deg, '-r',  label="{}".format(name1))
        if margins:
            plt.vlines([wp1, wg1],
                    np.min(phase_deg),
                    np.max(phase_deg),
                    linestyle='--')
            plt.hlines(-180, omega1[0], omega1[-1], linestyle='--')
        plt.legend()
        plt.grid()
    """
    plt.subplot(212)
    """     if margins:
        plt.hlines(0, omega2[0], omega2[-1], linestyle='--')
        plt.vlines([wp2, wg2], np.min(mag_dB2), np.max(mag_dB2), linestyle='--') 
    """
    plt.semilogx(omega2, mag_dB2, '-b',  label="{}".format(name2))
    plt.xlabel('Frequência (Hz)')
    plt.ylabel('Magnitude (dB)')
    plt.legend()
    plt.grid()
    plt.tight_layout()

    """  
    if margins:
        plt.title(
            ' bode pm: {:0.2f} deg @{:0.2f} Hz gm: {:0.2f} dB @{:0.2f} Hz'.
            format(pm2, wp2,  20 * np.log10(gm2), wg2))
    else:
    """    
    #plt.title(name2 + ' bode')
    """ 
    plt.subplot(224)
    phase_deg = np.rad2deg(phase2)
    plt.semilogx(omega2, phase_deg, '-b',  label="{}".format(name2)) 
    if margins:
    plt.vlines([wp2, wg2],
                np.min(phase_deg),
                np.max(phase_deg),
                linestyle='--')
    plt.hlines(-180, omega2[0], omega2[-1], linestyle='--') 
    plt.legend()
    plt.grid() 
    """
'''
            Extrai os parêmetros do numerador e denominador de uma função de transferência em S
            Parameters
            ----------
            H : Função de transferência no dominio de frequência S
'''
def extract_parameters(H):
    # Converter a string em uma expressão simbólica
    expressao = sp.sympify(H)
    # Obter os numeradores e denominadores como objetos simbólicos
    num_Hol, den_Hol = expressao.as_numer_denom()
    # print("numerador em s:", num_Hol)
    # print("denominador em s:", den_Hol)
    # Extrair coeficientes do numerador e denominador de Hs
    num_Hol = sp.Poly(num_Hol, s).all_coeffs()
    den_Hol = sp.Poly(den_Hol, s).all_coeffs()
    # print("numerador:", num_Hol)
    # print("denominador:", den_Hol)
    # Transforma em array para aplicar a transformada de laplace
    num = np.array(num_Hol, dtype=float)
    den = np.array(den_Hol, dtype=float)
    return num, den


'''
            Calcula o filtro IRR de acordo com os coeficientes recebidos
            Parameters
            ----------
            lp : array de coeficientes
'''


def IRR_filter(lp):
    IRR_order = len(lp)
    IRR = 1
    for i in range(IRR_order - 1):
        IRR *= (1 + s / fr) / (1 + (s / (lp[i] * fr)))
    IRR *= (1 + s / fr) / (1 + (s / (lp[IRR_order - 1] * fr)))

    return IRR


'''
            Plota o gráfico de Bode na função tranferência no dominío S 
            Parameters
            ----------
            H : (opcional) Função de transferência no dominío de frequência S
            numerador : (opcional) array dos parâmetros do numerador da Função de transferência sem "S"
            denominador : (opcional) array dos parâmetros do denominador da Função de transferência sem "S"
            irr : (opcional) array de coeficientes do filtro IRR
            w : (opcional) Lista de frequências para ser usada pela resposta em frequência
            margem : (opcional) Plotar a margem e ganho de fase
            prints : (opcional) Printar informações
'''


def h_to_bode(H=None, numerador=None, denominador=None, irr=None, freq=None, margem=False, prints=False, plot=False):
    if irr:
        IRR = IRR_filter(irr)
    else:
        IRR = 1
    if H:
        num_irr, den_irr = extract_parameters(IRR)
        IRR = control.tf(num_irr, den_irr)
        num_H, den_H = extract_parameters(H)
        tf = control.tf(num_H, den_H)
    else:
        num_irr, den_irr = extract_parameters(IRR)
        IRR = control.tf(num_irr, den_irr)
        tf = control.tf(numerador, denominador)
    if prints:
        print("Função de transfêrencia sem filtro IRR: ", tf)
    tf = tf * IRR
    # tf = (tf /(1 + tf))
    # tf = (1 / (1 + tf))
    # tf = control.feedback(tf, 1)
    if prints:
        print("função de transferencia em Laplace com filtro IRR: ", tf)
    if freq is not None:
        mag, phase, f = control.bode(tf, omega=freq, Hz=True, dB=True, deg=True, margins=margem, plot=plot)
    else:
        mag, phase, f = control.bode(tf, Hz=True, dB=True, deg=True, margins=margem, plot=plot)

    return mag, phase, f, tf


'''
        DEFINIÇÕES
'''
# Símbolo 's' para a variável de Laplace
s = sp.symbols('s')
N = 4.8e9/26e6  # relação de f/f_r

# Livro Bogdan PG 137/152
# fr = 26e6  # Frequeência de referência
# a = 2 ** -7  # alpha value
# p = 2 ** -15  # rho value
# #   Coeficiêntes do filtro IIR
# l = [2**-3, 2**-3, 2**-3, 2**-4]

# Aassignment 4
# fr = 40e6  # Frequeência de referência
fr = 33e6  # Frequeência de referência
# a = 2 ** -7  # alpha value
a = 2 ** -6  # alpha value
# p = 2 ** -14  # rho value
p = 2 ** -15  # rho value
# Coeficiêntes do filtro IIR
# l = [2 ** -2, 2 ** -3, 2 ** -2, 2 ** -3]
l = [2 ** -2 , 2 ** -1 , 2 ** -1 , 2 ** -1]

# Open Loop Unit Gain
w1 = a * fr * ( 0.5 + 0.5 * np.sqrt(1 + (4 * p / a**2)))

print("Frequência de ganho unitarío é: ", w1)
# Função de tranferência de loop aberto
HOL = (a + p * fr / s) * (fr / s)
# HOL = a * (fr/s)

if __name__ == "__main__":
    '''
    loop BandWidth = alpha/2pi * FREF
    
    '''
    w = np.logspace(3, 8, 100000)  # List of frequencies in rad/sec to be used for frequency response ( 10^-1 até 10^3)
    margem = True
    # show the response in frequency of signal sampled
    # # Plotar bode da TF
    mag, phase, f, Hol = h_to_bode(H=HOL, freq=w, irr=None, margem=margem, prints=True, plot=False)
    # ax1, ax2 = plt.gcf().axes  # get subplot axes
    # ax1.set_title('Magnitude Open Loop Hol(s)')
    # ax2.set_title('Phase Open Loop Hol(s)')

    mag, phase, f, Hol_irr = h_to_bode(H=HOL, freq=w, irr=l, margem=margem, prints=True, plot=False)
    # plt.figure()
    format_plot("Open Loop Hol(s)" ,"Open Loop Hol(s) + IRR", Hol, Hol_irr, w, margins=True)

    plt.figure()
    # # Função de tranferência de loop fechado para referência
    Hcl_ref = N * (Hol / (1 + Hol))
    Hcl_ref_irr = N * (Hol_irr / (1 + Hol_irr))
    format_plot("Closed Loop Hcl reference(s)", "Closed Loop Hcl reference(s) + IRR", Hcl_ref, Hcl_ref_irr, w, margins=True)

    plt.figure()
    # # Função de tranferência de loop fechado para o TDC
    Hcl_TDC = (Hol / (1 + Hol))
    Hcl_TDC_irr = (Hol_irr / (1 + Hol_irr))
    format_plot("$H_{TDC}$", "$H_{TDC}$ + IRR", Hcl_TDC, Hcl_TDC_irr, w, margins=True)
    """ 
        plt.figure()
        t, y = control.step_response(Hcl_TDC)

        # Plotar a resposta ao degrau
        plt.plot(t, y)
        plt.title('Resposta ao Degrau')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Resposta')
        plt.grid()
        plt.show()

    """


    plt.figure()
    # # Função de tranferência de loop fechado para o DCO
    Hcl_DCO = (1 / (1 + Hol))
    Hcl_DCO_irr = (1 / (1 + Hol_irr))
    format_plot("$H_{DCO}$", "$H_{DCO}$ + IRR", Hcl_DCO, Hcl_DCO_irr, w, margins=True)

    plt.figure()
    num = [1]
    den = [1, 0,0,0]  # Denominador é s
    tf = control.TransferFunction(num, den)
    Hcl_ADPLL = Hcl_DCO * Hcl_TDC
    Hcl_ADPLL_irr = Hcl_DCO_irr * Hcl_TDC_irr * Hcl_ref_irr
    format_plot("Ruído de fase do DCO", "Ruído de fase do ADPLL + IRR", Hcl_ADPLL, Hcl_ADPLL_irr, w, margins=False, offset=0)
    # plt.figure()
    # control.bode(Hcl_DCO_irr, w, Hz=True, dB=True, deg=False, Plot=True, margins=False)
    """     mag, phase, f = control.bode(Hcl_ADPLL, omega=w, Hz=True, dB=False, deg=True, margins=margem, plot=False)
    Xdb_o , f = fun_calc_psd(mag , fr , 100e3 , 1e3)
    plt.figure()
    plt.semilogx(f , Xdb_o , label="Phase Noise")
    plt.grid(visible=True)
    plt.legend()
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Phase Noise [dBc/Hz]')
    """
    plt.tight_layout()
    plt.show()
