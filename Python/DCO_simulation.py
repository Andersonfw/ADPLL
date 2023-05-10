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


def diff_freq(f):
    '''
        Calcula a diferença de frequêcia do valor desejado ao atual

        INPUT arguments
        f     :  Valor da frequência atual
        '''
    return abs((FREF * FCW) - f)


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
PVT_NB = 8  # número de bits em PVT mode
ACQ_NB = 8  # número de bits em acquisition mode
TRK_NB_I = 6  # número de bits Trekking integer mode
TRK_NB_F = 5  # número de bits Trekking fractional mode
FR_PVT = 500e6  # range de frequência em PVT mode
FR_ACQ = 100e6  # range de frequência em acquisition mode
FR_TRK_I = 2e6  # range de frequência em Trekking integer mode
FR_TRK_F = FR_TRK_I/ 2 ** TRK_NB_I # range de frequência em Trekking fractional mode
FREQ_RES_PVT = 0
FREQ_RES_ACQ = 0
FREQ_RES_TRK = 0
FREQ_RES_TRK_F = 0
L = 1e-9  # Indutor utilizado
C0 = 0  # valor de capacitância inicial
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
        LOOP FILTER
'''
Kp_PVT = 2 ** -2
Kp_ACQ = 2 ** -5
Kp_TRK = 2 ** -5
Ki_TRK = 2 ** -11

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
t_CKV = [0]  # period of DCO output
t_R = 0  # time reference
TDEV_I = 0  # time deviation integer
TDEV_F = 0  # time deviation fractional
jitter = [0]
wander = [0]
error_fractional = np.zeros(TIME)

pvt_bank_calib = False
acq_bank_calib = False
trk_bank_calib = False
OTW_pvt = 0  # initial value of pvt bank
OTW_acq = 128  # initial value of acq bank
OTW_trk = 32  # initial value of trk integer bank
OTW_trk_f = 0  # initial value of trk fractional bank
phase_dif = 0  # phase difference
prev_phase = 0  # new phase difference
count = 0  # counter of k index
k = 1  # index k
n = 0  # index n
freqs = np.zeros(TIME)  # array of different DCO output values
NTW = np.zeros(TIME)  # normalize tuning word
'''
        Main
'''
if __name__ == "__main__":
    Init_DCO()
    f_CKV = SET_DCO(255, 255, 64, 0)
    T0 = 1 / f_CKV
    fs = OVERSAMPLE * f_CKV
    print("frequência inicial do DCO é: ", f_CKV / 1e6, "MHz")
    dco_init_time = np.arange(0, 10 * T0, 1 / fs)
    dco_freq_init = np.sin(2 * np.pi * f_CKV * dco_init_time)
    # plot_histogram_noise(10000)
    # plot_DCO_signal()

    KDCO = FREQ_RES_PVT
    OTW = OTW_pvt
    F_start_DCO = f_CKV
    phase_error = np.zeros(TIME)

    div = 0
    NumberSamples = TIME * FCW
    BusSize = 5  # bits
    Fraction = 21  # usable 0 to 1
    FractionInternal = 2 ** BusSize * Fraction
    AccumulatorBits = 21  # bits
    AccumulatorSize = 2 ** AccumulatorBits

    C1 = np.zeros(NumberSamples)  # Carry out of the first accumulator
    C2 = np.zeros(NumberSamples)  # Carry out of the 2nd accumulator
    C3 = np.zeros(NumberSamples)  # Carry out of the 3nd accumulator
    U1 = np.zeros(NumberSamples)  # output of the 1st accum
    U2 = np.zeros(NumberSamples)  # output of the 2nd accum
    U3 = np.zeros(NumberSamples)  # output of the 3rd accum
    Yout1 = np.zeros(NumberSamples)  # output to the divider for 1 stage SDM
    Yout2 = np.zeros(NumberSamples)  # output to the divider for 2 stage SDM
    Yout3 = np.zeros(NumberSamples)  # output to the divider for 3 stage SDM
    # out = np.zeros(NumberSamples)
    index = 2
    for k in range(1, TIME):
        RR_k += FCW  # reference phase accumulator
        t_R = k * FREF_edge
        while t_CKV[n] < t_R:
            n += 1
            RV_n = n  # variable phase accumulator
            # delta_f = f_CKV - (F_start_DCO + (KDCO * (OTW)))
            # delta_f = f_CKV - FDCO
            # TDEV_I = delta_f / (f_CKV * (f_CKV + delta_f))
            if trk_bank_calib:
                div+=1
                if div == 4:
                    index +=1
                    div = 0
                    FractionInternal = 2 ** BusSize * error_fractional[k - 1]
                    U1[index] = FractionInternal + U1[index - 1]
                    U2[index] = U1[index - 1] + U2[index - 1]
                    # U2[index] = U1[index-1] + U2[index-1]
                    U3[index] = U2[index - 1] + U3[index - 1]
                    if U1[index] > AccumulatorSize:
                        C1[index] = 1  # carry 1
                        U1[index] -= AccumulatorSize
                    if U2[index] > AccumulatorSize:
                        C2[index] = 1  # carry 2
                        U2[index] -= AccumulatorSize
                    if U3[index] > AccumulatorSize:
                        C3[index] = 1  # carry 3
                        U3[index] -= AccumulatorSize
                    out = C1[index - 3] + C2[index - 2] - C2[index - 3] + C3[index - 1] - 2 * C3[index - 2] + C3[index - 3]
                    OTW_trk_f +=out
                    # f_CKV = SET_DCO(OTW_pvt, OTW_acq, OTW_trk, OTW_trk_f)
                    # T0 = 1 / f_CKV

            jitter.append(np.random.randn() * Jt_noise)
            wander.append(np.random.randn() * Wt_noise)
            t_CKV.append(n * T0 + jitter[n] + wander[n] - jitter[n - 1])  # - TDEV_I

        RV_k = RV_n  # variable phase accumulator
        error_fractional[k] = TDC(t_R, t_CKV[n - 1])  # TDC
        phase_error[k] = (RR_k - RV_k + error_fractional[k])  # Phase detector

        ##################### PVT MODE #################################################
        if not pvt_bank_calib:
            NTW[k] = OTW_pvt + (int(phase_error[k]) * Kp_PVT)  # calcula o novo valor de NTW como inteiro
            OTW_pvt = NTW[k]  # ajusta o novo valor de controle dos capacitores do PVT bank
            if diff_freq(f_CKV) <= FREQ_RES_PVT:
                count += 1
                if count == 5:
                    pvt_bank_calib = True
                    count = 0
            else:
                count = 0
        ##################### ACQUISITION MODE #########################################
        elif not acq_bank_calib:
            NTW[k] = OTW_acq + (int(phase_error[k]) * Kp_ACQ)  # calcula o novo valor de NTW como inteiro
            OTW_acq = NTW[k]  # ajusta o novo valor de controle dos capacitores do ACQ bank
            if diff_freq(f_CKV) <= FREQ_RES_ACQ:
                count += 1
                if count == 5:
                    acq_bank_calib = True
                    trk_bank_calib = True
            else:
                count = 0

        ##################### TREKING MODE ################################################

        elif trk_bank_calib:
            NTW[k] = OTW_trk + ((int(phase_error[k])) * Kp_TRK)  # calcula o novo valor de NTW como inteiro
            # NTW[k] = phase_error[k] * Kp_TRK - Kp_TRK * phase_error[k - 1] + Ki_TRK * phase_error[k - 1] + NTW[k - 1]
            OTW_trk = NTW[k]  # ajusta o novo valor de controle dos capacitores do TRK bank
        f_CKV = SET_DCO(OTW_pvt, OTW_acq, OTW_trk, OTW_trk_f)
        T0 = 1 / f_CKV
        freqs[k] = f_CKV  # insere o valor de frequência ajustado no index k

    print("freq ajustada: ", f_CKV / 1e6,
          "MHz E a desejada era de :", (FREF * FCW) / 1e6,
          "MHz diferença de :",(f_CKV - (FREF * FCW)) / 1e3, "kHz")

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
