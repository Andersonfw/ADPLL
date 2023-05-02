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


def Init_DCO():
    pvt_bank = LSB_BANK(F0, PVT_NB, FR_PVT)
    acq_bank = LSB_BANK(F0, ACQ_NB, FR_ACQ)
    trk_i_bank = LSB_BANK(F0, TRK_NB_I, FR_TRK_I)
    trk_f_bank = LSB_BANK(F0, TRK_NB_F, FR_TRK_F)
    global C0, pvt_lsb, acq_lsb, trk_i_lsb, trk_f_lsb
    C0 = pvt_bank.cmin
    # C0=10e-4
    pvt_lsb = pvt_bank.lsb
    acq_lsb = acq_bank.lsb
    trk_i_lsb = trk_i_bank.lsb
    trk_f_lsb = trk_i_bank.lsb / (2 ** TRK_NB_F)
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
    global C0, pvt_lsb, acq_lsb, trk_i_lsb, trk_f_lsb
    pvt = 255 - int(pvt_OTW)
    acq = 255 - int(acq_OTW)
    trk_i = 64 - int(trk_i_OTW)
    trk_f = 32 - int(trk_f_OTW)
    f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb + trk_f_lsb * trk_f)))
    return f

def TDC (T_ckv, res_TDC, ):

    return 0


def fun_calc_psd(x, fs=1, rbw=100e3, fstep=None):
    # Calculate power spectral density
    #
    # INPUT arguments
    # x     :  data vector
    # fs    :  sample rate [Hz]
    # rbw   :  resolution bandwidth [Hz]
    # fstep :  FFT frequency bin separation [Hz]
    # OUTPUT
    # XdB	: spectrum of x [dB]
    # f	: frequency vector [Hz]

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
FCW = 76#76.9230  # Frequency command word
FREF = 26e6  # Frequência de referência
FDCO = FREF * FCW  # Frequência desajda na sáida do DCO
FREF_edge = 1 / FREF  # tempo de borda de FREF
FDCO_edge = 1 / FDCO  # tempo de borda de F0

FR_PVT = 500e6  # range de frequência em PVT mode
FR_ACQ = 100e6  # range de frequência em acquisition mode
FR_TRK_I = 2e6  # range de frequência em Trekking integer mode
FR_TRK_F = 2e6  # range de frequência em Trekking fractional mode
L = 1e-9  # Indutor utilizado

PVT_NB = 8  # número de bits em PVT mode
ACQ_NB = 8  # número de bits em acquisition mode
TRK_NB_I = 6  # número de bits Trekking integer mode
TRK_NB_F = 5  # número de bits Trekking fractional mode
OVERSAMPLE = 100  # oversample de frequência para discretizar a frequência do DCO

TDC_res = 15e-12
TDC_chains = 40

Wt_noise = 12e-15  # Wander noise time
Wf_noise = 1 / Wt_noise  # Wander noise frequency
Jt_noise = 111e-15  # jitter noise time
Jf_noise = 1 / Jt_noise  # jitter noise frequency

TIME = 2000  # simulação de X bordas de FREF

'''
        VARIÁVEIS GLOBAIS
'''
C0 = 0  # valor de capacitância inicial
pvt_lsb = 0  # valor do LSB em PVT mode
acq_lsb = 0  # valor do LSB em acquisition mode
trk_i_lsb = 0  # valor do LSB em Trekking integer mode
trk_f_lsb = 0  # valor do LSB em Trekking fractional mode

'''
        Main
'''
if __name__ == "__main__":
    Init_DCO()
    f_CKV = SET_DCO(128, 128, 32, 0)
    T0 = 1 / f_CKV
    fs = OVERSAMPLE * f_CKV
    print("freq: ", f_CKV)
    t_init = np.arange(0, 10 * T0, 1 / fs)
    # jitter = np.random.randn(len(t)) * Jt_noise
    # wander = np.random.randn(len(t)) * Wt_noise
    # plt.subplot(121)
    # plt.hist(jitter, bins=100, label="Normal distribution of the Jitter noise")
    # plt.legend()
    # plt.grid(visible=True)
    # plt.subplot(122)
    # plt.hist(watter, bins=100, color="r", label="Normal distribution of the wander noise")
    # plt.legend()
    # plt.grid(visible=True)
    # plt.show()

    # x = np.sin(2 * np.pi * 1/(tc + jitter + wander) * t)
    x_init = np.sin(2 * np.pi * f_CKV * t_init)
    len_simulation = 6 * OVERSAMPLE  # plotar 6 periodos do DCO
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
    #
    # plt.figure()
    # plt.plot(t[:len_simulation] / 1e-9, x[:len_simulation], label="V(t)")
    # plt.plot(t[:len_simulation] / 1e-9, x_or[:len_simulation], label="V(t) Original")
    # plt.grid(visible=True)
    # plt.legend()
    # plt.xlabel('Time (ns)')
    # plt.ylabel('Amplitude (V)')
    # plt.show()

    RR_k = 0
    RV_n = 0
    RV_k = 0
    t_CKV = 0
    t_R = 0
    TDEV_I = 0
    TDEV_F = 0
    last_jitter = 0
    last_error = 0

    pvt_bank_calib = False
    acq_bank_calib = False
    trk_bank_calib = False
    OTW_pvt = 128
    OTW_acq = 128
    OTW_trk = 32
    OTW_trk_f = 0
    phase_dif = 0
    prev_phase = 0
    count = 0
    k = 1
    n = 0
    freqs = np.zeros(TIME)

    for k in range(1, TIME):
        RR_k += FCW
        t_R = k * FREF_edge
        while t_CKV < t_R:
            n += 1
            delta_f = f_CKV - FDCO
            TDEV_I = delta_f / (FDCO * (FDCO + delta_f))
            jitter = 0
            wander = 0
            # jitter = np.random.randn() * Jt_noise
            # wander = np.random.randn() * Wt_noise
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
        error_TDC = TDC(k, n)
        delta_tR = t_CKV - t_R
        error_fractional = 1 - delta_tR / T0
        phase_error = RR_k - RV_k + error_fractional
        if not pvt_bank_calib:
            OTW_prev = OTW_pvt + (int(phase_error) * 2 ** -4)
            OTW_pvt = OTW_prev
            if k == 200:
                pvt_bank_calib = True

        elif not acq_bank_calib:
            OTW_prev = OTW_acq + (int(phase_error) * 2 ** -4)
            OTW_acq = OTW_prev
            if k == 400:
                acq_bank_calib = True
                trk_bank_calib = True

        elif trk_bank_calib:
            OTW_prev = OTW_trk + (int(phase_error) * 2 ** -13)
            OTW_trk = OTW_prev
        f_CKV = SET_DCO(OTW_pvt, OTW_acq, OTW_trk, OTW_trk_f)
        T0 = 1 / f_CKV
        freqs[k] = f_CKV
    print("freq ajustada: ", f_CKV, " E a desejada era de :", FREF * FCW, "diferença de :", f_CKV - (FREF * FCW))

    plt.figure()
    plt.plot(np.arange(1, TIME, 1), freqs[1:TIME])
    plt.grid(visible=True)

    fs = OVERSAMPLE * f_CKV
    print("freq: ", f_CKV)
    t = np.arange(0, 10 * T0, 1 / fs)
    # x = np.sin(2 * np.pi * 1/(tc + jitter + wander) * t)
    x_or = np.sin(2 * np.pi * f_CKV * t)

    plt.figure()
    plt.plot(t_init[:len_simulation] / 1e-9, x_init[:len_simulation], label="V(t) Inicial")
    plt.plot(t[:len_simulation] / 1e-9, x_or[:len_simulation], label="V(t)")
    plt.grid(visible=True)
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (V)')
    plt.show()

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
