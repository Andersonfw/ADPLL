"""
Created on abril 26 18:02:04 2023

@author: Ânderson Felipe Weschenfelder
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


class LSB_BANK:
    def __init__(self, fc, nb, fr):
        self.fc = fc
        self.nb = nb
        self.fr = fr
        self.fmin = fc - fr / 2
        self.fmax = fc + fr / 2
        self.cmax = 1 / (L * (2 * np.pi * self.fmin) ** 2)
        self.cmin = 1 / (L * (2 * np.pi * self.fmax) ** 2)
        self.lsb = (self.cmax - self.cmin) / 2 ** nb
        self.freq_lsb = fr / 2 ** nb


def Init_DCO():
    '''
    Configuração inicial do DCO
    '''

    pvt_bank = LSB_BANK(F0, PVT_NB, FR_PVT)
    acq_bank = LSB_BANK(F0, ACQ_NB, FR_ACQ)
    trk_i_bank = LSB_BANK(F0, TRK_NB_I, FR_TRK_I)
    trk_f_bank = LSB_BANK(F0, TRK_NB_F, FR_TRK_F)
    global C0, pvt_lsb, acq_lsb, trk_i_lsb, trk_f_lsb, FREQ_RES_PVT, FREQ_RES_ACQ, FREQ_RES_TRK, FREQ_RES_TRK_F
    C0 = pvt_bank.cmin
    # C0=10e-4
    pvt_lsb = pvt_bank.lsb
    acq_lsb = acq_bank.lsb
    trk_i_lsb = trk_i_bank.lsb
    trk_f_lsb = trk_i_bank.lsb / (2 ** TRK_NB_F)
    FREQ_RES_PVT = pvt_bank.freq_lsb
    FREQ_RES_ACQ = acq_bank.freq_lsb
    FREQ_RES_TRK = trk_i_bank.freq_lsb
    FREQ_RES_TRK_F = trk_f_bank.freq_lsb
    print("PVT -----  fmin: ", pvt_bank.fmin, " fmax: ", pvt_bank.fmax, "cmax: ", pvt_bank.cmax, " cmim: ",
          pvt_bank.cmin, " LSB: ", pvt_bank.lsb)
    print("ACQ -----  fmin: ", acq_bank.fmin, " fmax: ", acq_bank.fmax, "cmax: ", acq_bank.cmax, " cmim: ",
          acq_bank.cmin, " LSB: ", acq_bank.lsb)
    print("TRK_I -----  fmin: ", trk_i_bank.fmin, " fmax: ", trk_i_bank.fmax, "cmax: ", trk_i_bank.cmax, " cmim: ",
          trk_i_bank.cmin, " LSB: ", trk_i_bank.lsb)
    print("TRK_F -----  fmin: ", trk_f_bank.fmin, " fmax: ", trk_f_bank.fmax, "cmax: ", trk_f_bank.cmax, " cmim: ",
          trk_f_bank.cmin, " LSB: ", trk_f_bank.lsb)
    return


def SET_DCO(pvt_OTW=255, acq_OTW=255, trk_i_OTW=64, trk_f_OTW=0):
    '''
    Ajusta a frequência do DCO conforme valor binario de cada bank capacitor

    INPUT arguments
    pvt_OTW     :  Valor inteiro do bank pvt
    acq_OTW     :  Valor inteiro do bank acquisition
    trk_i_OTW   :  Valor inteiro do bank trekking integer
    trk_f_OTW   :  Valor inteiro do bank trekking fractional

    OUTPUT
    f	: valor de frequência [Hz]
    '''

    global C0, pvt_lsb, acq_lsb, trk_i_lsb, trk_f_lsb
    pvt = 255 - int(pvt_OTW)
    acq = 255 - int(acq_OTW)
    trk_i = 64 - int(trk_i_OTW)
    trk_f = 32 - int(trk_f_OTW)
    f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb + trk_f_lsb * trk_f)))
    return f


def TDC(tR, t_ckv):
    global TDC_res, T0
    tR_Q = int((
                           tR - t_ckv) / TDC_res)  # Diferença de tempo entre a última borda de clock de CKV até a borda de REF. (FIG 2 Time-Domain Modeling of an RF All-Digital PLL)
    # delta_tR = int(((t_CKV - ntdc_init)/(n - n_init)) / TDC_res)
    error = 1 - (tR_Q * TDC_res) / T0
    return error


def plot_DCO_signal():
    '''
    plotar sinal de saída do DCO e comparar ao sinal de inicio

    INPUT arguments
    none
    '''
    global len_simulation, dco_init_time, dco_freq_init, f_CKV, T0, OVERSAMPLE, fs
    fs = OVERSAMPLE * f_CKV
    t = np.arange(0, 10 * T0, 1 / fs)
    # x_noise = np.sin(2 * np.pi * 1/(tc + jitter + wander) * t)
    x = np.sin(2 * np.pi * f_CKV * t)
    plt.figure()
    plt.plot(dco_init_time[:len_simulation] / 1e-9, dco_freq_init[:len_simulation], label="DCO Inicial")
    plt.plot(t[:len_simulation] / 1e-9, x[:len_simulation], label="DCO out")
    plt.grid(visible=True)
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (V)')
    plt.show()


def plot_histogram_noise(lenght):
    '''
    plotar histograma do ruído Wander e jitter

    INPUT arguments
    lenght     :  Quantidade de pontos aleatórios
    '''
    global Jt_noise, Wt_noise
    jitter_noise = np.random.randn(lenght) * Jt_noise
    wander_noise = np.random.randn(lenght) * Wt_noise
    plt.subplot(121)
    plt.hist(jitter_noise, bins=50, label="Normal distribution of the Jitter noise")
    plt.legend()
    plt.grid(visible=True)
    plt.subplot(122)
    plt.hist(wander_noise, bins=50, color="r", label="Normal distribution of the wander noise")
    plt.legend()
    plt.grid(visible=True)
    plt.show()


def fun_calc_psd(x, fs=1, rbw=100e3, fstep=None):
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
    fftstr = f'len(x)={len_x:.2f}, rbw={rbw / 1e3:.2f}kHz, fstep={fstep / 1e3:.2f}kHz, nfft={nfft:d}, nwin={nwin:d}'
    print(f'Calculating the PSD: {fftstr} ...')
    f, X = signal.welch(x, fs=fs, window=signal.windows.blackman(nwin), nperseg=nwin, nfft=nfft, scaling='density')
    X *= (np.sinc(f / fs)) ** 2  # correct for ZOH
    XdB = 10 * np.log10(X)
    XdB_sig = np.max(XdB)
    print(f'Signal PSD peak = {XdB_sig:.2f} dB, 10log(rbw) = {10 * np.log10(rbw):.1f}')
    return XdB, f


'''
        DEFINIÇÕES GERAIS
'''
F0 = 2045e6  # frequência central de ajuste do DCO
FCW = 72  # 76.9230  # Frequency command word
FREF = 26e6  # Frequência de referência
FDCO = FREF * FCW  # Frequência desajda na sáida do DCO
FREF_edge = 1 / FREF  # tempo de borda de FREF
FDCO_edge = 1 / FDCO  # tempo de borda de F0

'''
        BANK CAPACITOR
'''
FR_PVT = 500e6  # range de frequência em PVT mode
FR_ACQ = 100e6  # range de frequência em acquisition mode
FR_TRK_I = 2e6  # range de frequência em Trekking integer mode
FR_TRK_F = 2e6  # range de frequência em Trekking fractional mode
FREQ_RES_PVT = 0
FREQ_RES_ACQ = 0
FREQ_RES_TRK = 0
FREQ_RES_TRK_F = 0
L = 1e-9  # Indutor utilizado
C0 = 0  # valor de capacitância inicial
PVT_NB = 8  # número de bits em PVT mode
ACQ_NB = 8  # número de bits em acquisition mode
TRK_NB_I = 6  # número de bits Trekking integer mode
TRK_NB_F = 5  # número de bits Trekking fractional mode
pvt_lsb = 0  # valor do LSB em PVT mode
acq_lsb = 0  # valor do LSB em acquisition mode
trk_i_lsb = 0  # valor do LSB em Trekking integer mode
trk_f_lsb = 0  # valor do LSB em Trekking fractional mode

'''
        TDC
'''
TDC_res = 15e-12
TDC_chains = 40

'''
        RUÍDO
'''
Wt_noise = 12e-15  # Wander noise time
Wf_noise = 1 / Wt_noise  # Wander noise frequency
Jt_noise = 111e-15  # jitter noise time
Jf_noise = 1 / Jt_noise  # jitter noise frequency

'''
        VARIÁVEIS DE CONTROLE DA SIMULAÇÃO
'''
TIME = 2000  # simulação de X bordas de FREF
OVERSAMPLE = 100  # oversample de frequência para discretizar a frequência do DCO
len_simulation = 6 * OVERSAMPLE  # plotar 6 períodos do DCO

'''
        VARIÁVEIS DE CONTROLE
'''
RR_k = 0  # reference phase
RV_n = 0  # variable phase with index n
RV_k = 0  # variable phase with index k
t_CKV = 0  # period of DCO output
t_R = 0  # time reference
TDEV_I = 0  # time deviation integer
TDEV_F = 0  # time deviation fractional
last_jitter = 0  # last value of jitter
last_error = 0  # last value os error

pvt_bank_calib = False
acq_bank_calib = False
trk_bank_calib = False
OTW_pvt = 128  # initial value of pvt bank
OTW_acq = 128  # initial value of acq bank
OTW_trk = 32  # initial value of trk integer bank
OTW_trk_f = 0  # initial value of trk fractional bank
phase_dif = 0  # phase difference
prev_phase = 0  # new phase difference
count = 0  # counter of k index
k = 1  # index k
n = 0  # index n
freqs = np.zeros(TIME)  # array of different DCO output values

'''
        Main
'''
if __name__ == "__main__":
    Init_DCO()
    f_CKV = SET_DCO(128, 128, 32, 0)
    T0 = 1 / f_CKV
    fs = OVERSAMPLE * f_CKV
    print("frequência inicial do DCO é: ", f_CKV / 1e6, "MHz")
    dco_init_time = np.arange(0, 10 * T0, 1 / fs)
    dco_freq_init = np.sin(2 * np.pi * f_CKV * dco_init_time)
    # plot_histogram_noise(10000)
    # plot_DCO_signal()

    for k in range(1, TIME):
        RR_k += FCW
        t_R = k * FREF_edge
        ntdc_init = t_CKV
        n_init = n
        while t_CKV < t_R:
            n += 1
            # delta_f = KDCO * OTW
            delta_f = f_CKV - FDCO
            TDEV_I = delta_f / (FDCO * (FDCO + delta_f))
            jitter = 0
            wander = 0
            jitter = np.random.randn() * Jt_noise
            wander = np.random.randn() * Wt_noise
            last_t_CKV = t_CKV  # aramazena o valor anterior de t_CKV
            t_CKV = n * T0 + jitter + wander - last_jitter  # - TDEV_I
            last_jitter = jitter
            RV_n += 1
            # if trk_bank_calib:
            #     count += 1
            # if count == 4:
            #     count = 0
            #     error_f = phase_error % 1 #+ last_error
            #     last_error = error_f % 1
            #     OTW_trk_f = error_f * 2** TRK_NB_F
        RV_k = RV_n
        # error_TDC = TDC(t_R, last_t_CKV)
        # delta_tR = int((t_R - last_t_CKV) / TDC_res)  # Diferença de tempo entre a última borda de clock de CKV até a borda de REF. (FIG 2 Time-Domain Modeling of an RF All-Digital PLL)
        # delta_tR = int(((t_CKV - ntdc_init)/(n - n_init)) / TDC_res)
        # error_fractional = 1 - (delta_tR * TDC_res) / T0
        error_fractional = TDC(t_R, last_t_CKV)
        phase_error = RR_k - RV_k + error_fractional
        if not pvt_bank_calib:
            OTW_prev = OTW_pvt + (int(phase_error) * 2 ** -2)
            OTW_pvt = OTW_prev
            if k >= 158:
                print("debug")
            if k == 150:
                pvt_bank_calib = True
            num = abs((FREF * FCW) - f_CKV)
            if num <= FREQ_RES_PVT:
                count += 1
                if count == 5:
                    pvt_bank_calib = True
                    count = 0
            else:
                count = 0

        elif not acq_bank_calib:
            OTW_prev = OTW_acq + (int(phase_error) * 2 ** -4)
            OTW_acq = OTW_prev
            # if k == 400:
            #     acq_bank_calib = True
            #     trk_bank_calib = True
            num = abs((FREF * FCW) - f_CKV)
            if num <= FREQ_RES_ACQ:
                count += 1
                if count == 5:
                    acq_bank_calib = True
                    trk_bank_calib = True
            else:
                count = 0

        elif trk_bank_calib:
            OTW_prev = OTW_trk + (int(phase_error) * 2 ** -13)
            OTW_trk = OTW_prev
        f_CKV = SET_DCO(OTW_pvt, OTW_acq, OTW_trk, OTW_trk_f)
        T0 = 1 / f_CKV
        freqs[k] = f_CKV
        # if f_CKV < (FREF * FCW):
        #     print("freq menor")
    print("freq ajustada: ", f_CKV / 1e6, "MHz E a desejada era de :", (FREF * FCW) / 1e6, "MHz diferença de :",
          (f_CKV - (FREF * FCW)) / 1e3, "kHz")

    plt.figure()
    plt.plot(np.arange(1, TIME, 1), freqs[1:TIME])
    plt.grid(visible=True)

    plot_DCO_signal()

    #
    # Xdb_o, f = fun_calc_psd((x_or), fs, 1e3, 10e3)
    # Xdb, f = fun_calc_psd((x), fs, 1e3, 10e3)
    #
    # plt.figure()
    # # plt.plot(f / 1e6, Xdb_o, label="Original")
    # plt.semilogx(f, Xdb_o, label="Original")
    # # plt.plot(f / 1e6, Xdb, label="ERROR")
    # plt.semilogx(f, Xdb, label="With Noise")
    # plt.grid(visible=True)
    # plt.legend()
    # plt.xlabel('Freq (MHz)')
    # plt.ylabel('Amplitude (dB)')
    # # plt.show()
'''
For each f_ref cycle
    Accumulate FCW:
    Generate f_ref clock transition:
    While t [n] < t [k]
    ∑∆ -Modulator
    DCO
    Generate f_ckv clock transition:
    End
    TDC – calculate the timing difference
    Phase detector
    Loop filter
    Generate the current oscillator tuning word
End
'''
