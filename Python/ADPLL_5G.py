"""
Created on setembro 12 22:53:30 2023

@author: Ânderson Felipe Weschenfelder
"""
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import decimal
import locale
import datetime
import time
import os
import glob
import pandas as pd
from matplotlib.backend_bases import MouseButton
import random

class LSB_BANK:
    def __init__(self , fc , nb , fr):
        self.fc = fc
        self.nb = nb
        self.fr = fr
        self.fmin = fc - fr / 2
        self.fmax = fc + fr / 2
        self.cmax = 1 / (L * (2 * np.pi * self.fmin) ** 2)
        self.cmin = 1 / (L * (2 * np.pi * self.fmax) ** 2)
        self.lsb = (self.cmax - self.cmin) / 2 ** nb
        self.lsb_sup = self.lsb * (1 + MISMATCH_DCO)
        self.lsb_inf = self.lsb * (1 - MISMATCH_DCO)
        self.freq_lsb = fr / 2 ** nb

pvt_array = np.zeros(128)
acq_array = np.zeros(64)
trk_array = np.zeros(128)
def Init_DCO():
    '''
    Configuração inicial do DCO
    '''

    pvt_bank = LSB_BANK(F0_PVT , PVT_NB , FR_PVT)
    acq_bank = LSB_BANK(F0_ACQ , ACQ_NB , FR_ACQ)
    trk_i_bank = LSB_BANK(F0_TRK , TRK_NB_I , FR_TRK_I)
    trk_f_bank = LSB_BANK(F0_TRK , TRK_NB_F , FR_TRK_F)
    global C0 , pvt_lsb , acq_lsb , trk_i_lsb , trk_f_lsb , FREQ_RES_PVT , FREQ_RES_ACQ , FREQ_RES_TRK , FREQ_RES_TRK_F, pvt_lsb_inf, pvt_lsb_sup, acq_lsb_inf, acq_lsb_sup, trk_lsb_inf, trk_lsb_sup
    global pvt_array, acq_array, trk_array, MISMATCH_DCO
    C0 = pvt_bank.cmin
    pvt_lsb = pvt_bank.lsb
    pvt_lsb_sup = pvt_bank.lsb_sup
    pvt_lsb_inf = pvt_bank.lsb_inf
    acq_lsb = acq_bank.lsb
    acq_lsb_sup = acq_bank.lsb_sup
    acq_lsb_inf = acq_bank.lsb_inf
    trk_i_lsb = trk_i_bank.lsb
    trk_lsb_sup = trk_i_bank.lsb_sup
    trk_lsb_inf = trk_i_bank.lsb_inf
    trk_f_lsb = trk_i_bank.lsb / (2 ** TRK_NB_F)
    FREQ_RES_PVT = pvt_bank.freq_lsb
    FREQ_RES_ACQ = acq_bank.freq_lsb
    FREQ_RES_TRK = trk_i_bank.freq_lsb
    FREQ_RES_TRK_F = trk_f_bank.freq_lsb


    for i in range(128):
        var_pvt = np.random.normal(loc=0, scale = pvt_lsb * MISMATCH_DCO)
        pvt_array[i] = pvt_lsb + var_pvt
        var_trk = np.random.normal(loc=0, scale = trk_i_lsb * MISMATCH_DCO)
        trk_array[i] = trk_i_lsb + var_trk
        if i < 64:
            var_acq = np.random.normal(loc=0, scale = acq_lsb * MISMATCH_DCO)
            acq_array[i] = acq_lsb + var_acq
    # var_pvt = np.random.normal(loc=0, scale= pvt_bank.lsb * MISMATCH_DCO, size=2**PVT_NB)
    # pvt_array = np.full(2**PVT_NB, pvt_bank.lsb) + var_pvt
    # var_acq = np.random.normal(loc=0, scale= acq_bank.lsb * MISMATCH_DCO, size=2**ACQ_NB)
    # acq_array = np.full(2**ACQ_NB, acq_bank.lsb) + var_acq
    # var_trk = np.random.normal(loc=0, scale= trk_i_bank.lsb * MISMATCH_DCO, size=2**TRK_NB_I)
    # trk_array = np.full(2**TRK_NB_I, trk_i_bank.lsb) + var_trk



    print("PVT -----  fmin: " , pvt_bank.fmin , " fmax: " , pvt_bank.fmax , "cmax: " , pvt_bank.cmax , " cmim: " ,
          pvt_bank.cmin , " LSB: " , pvt_bank.lsb, " FREQ LSB: " , pvt_bank.freq_lsb, "MEAN: ", np.std(pvt_array))
    print("ACQ -----  fmin: " , acq_bank.fmin , " fmax: " , acq_bank.fmax , "cmax: " , acq_bank.cmax , " cmim: " ,
          acq_bank.cmin , " LSB: " , acq_bank.lsb, " FREQ LSB: " , acq_bank.freq_lsb, "MEAN: ", np.std(acq_array))
    print("TRK_I -----  fmin: " , trk_i_bank.fmin , " fmax: " , trk_i_bank.fmax , "cmax: " , trk_i_bank.cmax ,
          " cmim: " ,
          trk_i_bank.cmin , " LSB: " , trk_i_bank.lsb, " FREQ LSB: " , trk_i_bank.freq_lsb, "MEAN: ", np.std(trk_array))
    # print("TRK_F -----  fmin: ", trk_f_bank.fmin, " fmax: ", trk_f_bank.fmax, "cmax: ", trk_f_bank.cmax, " cmim: ",
    #       trk_f_bank.cmin, " LSB: ", trk_f_bank.lsb)
    return


def binaryValue(N_bits , OscilatorValue):
    ret = (2 ** N_bits - 1) - abs(int(OscilatorValue))
    if ret < 0:
        ret = 0
    return ret


def SET_DCO(pvt_OTW=255 , acq_OTW=255 , trk_i_OTW=64 , trk_f_OTW=0):
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

    global C0 , MISMATCH_DCO, pvt_lsb , acq_lsb , trk_i_lsb , trk_f_lsb , AVG_FCKV , avg_counter , freq_array , trk_bank_calib, pvt_lsb_inf, pvt_lsb_sup, acq_lsb_inf, acq_lsb_sup, trk_lsb_inf, trk_lsb_sup
    global pvt_bank_calib, acq_bank_calib, trk_bank_calib, desvio_padrao_pvt, desvio_padrao_acq, desvio_padrao_trk
    global pvt_array, acq_array, trk_array
    
    pvt = binaryValue(PVT_NB , pvt_OTW)
    acq = binaryValue(ACQ_NB , acq_OTW)
    trk_i = binaryValue(TRK_NB_I , trk_i_OTW)
    trk_f = 0  # binaryValue(TRK_NB_F, trk_f_OTW)


    if not pvt_bank_calib:
        desvio_padrao_pvt = np.random.normal(loc=0, scale = pvt_lsb * MISMATCH_DCO)
    
    elif not acq_bank_calib:
        desvio_padrao_acq = np.random.normal(loc=0, scale = acq_lsb * MISMATCH_DCO)

    elif trk_bank_calib:
        desvio_padrao_trk = np.random.normal(loc=0, scale = trk_i_lsb * MISMATCH_DCO)
    # f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * random.uniform(pvt_lsb_inf, pvt_lsb_sup) + acq * random.uniform(acq_lsb_inf, acq_lsb_sup) + trk_i * random.uniform(trk_lsb_inf, trk_lsb_sup) + trk_f_lsb * trk_f)))
    # f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb + trk_f_lsb * trk_f)))
    # f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * (pvt_lsb+desvio_padrao_pvt) + acq * (acq_lsb + desvio_padrao_acq) + trk_i * (trk_i_lsb + desvio_padrao_trk) + trk_f_lsb * trk_f)))
    f = 1 / (2 * np.pi * np.sqrt(L * (C0 + np.sum(pvt_array[:pvt]) + np.sum(acq_array[:acq]) + np.sum(trk_array[:trk_i]) )))

    freq_array[avg_counter] = f / TDC_DIVISOR   # accumulator of CKV after divisor to average in TDC
    avg_counter += 1
    if avg_counter == AVG_FCKV:
        avg_counter = 0
        
    return f
# def SET_DCO(pvt_OTW=255 , acq_OTW=255 , trk_i_OTW=64 , trk_f_OTW=0):
#     '''
#     Ajusta a frequência do DCO conforme valor binario de cada bank capacitor

#     INPUT arguments
#     pvt_OTW     :  Valor inteiro do bank pvt
#     acq_OTW     :  Valor inteiro do bank acquisition
#     trk_i_OTW   :  Valor inteiro do bank trekking integer
#     trk_f_OTW   :  Valor inteiro do bank trekking fractional

#     OUTPUT
#     f	: valor de frequência [Hz]
#     '''

#     global C0 , MISMATCH_DCO, pvt_lsb , acq_lsb , trk_i_lsb , trk_f_lsb , AVG_FCKV , avg_counter , freq_array , trk_bank_calib, pvt_lsb_inf, pvt_lsb_sup, acq_lsb_inf, acq_lsb_sup, trk_lsb_inf, trk_lsb_sup
#     global pvt_bank_calib, acq_bank_calib, trk_bank_calib, desvio_padrao_pvt, desvio_padrao_acq, desvio_padrao_trk
#     pvt = binaryValue(PVT_NB , pvt_OTW)
#     acq = binaryValue(ACQ_NB , acq_OTW)
#     trk_i = binaryValue(TRK_NB_I , trk_i_OTW)
#     trk_f = 0  # binaryValue(TRK_NB_F, trk_f_OTW)


#     if not pvt_bank_calib:
#         desvio_padrao_pvt = np.random.normal(loc=0, scale = pvt_lsb * MISMATCH_DCO)
    
#     elif not acq_bank_calib:
#         desvio_padrao_acq = np.random.normal(loc=0, scale = acq_lsb * MISMATCH_DCO)

#     elif trk_bank_calib:
#         desvio_padrao_trk = np.random.normal(loc=0, scale = trk_i_lsb * MISMATCH_DCO)
#     # f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * random.uniform(pvt_lsb_inf, pvt_lsb_sup) + acq * random.uniform(acq_lsb_inf, acq_lsb_sup) + trk_i * random.uniform(trk_lsb_inf, trk_lsb_sup) + trk_f_lsb * trk_f)))
#     #f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb + trk_f_lsb * trk_f)))
#     f = 1 / (2 * np.pi * np.sqrt(L * (C0 + pvt * (pvt_lsb+desvio_padrao_pvt) + acq * (acq_lsb + desvio_padrao_acq) + trk_i * (trk_i_lsb + desvio_padrao_trk) + trk_f_lsb * trk_f)))

#     freq_array[avg_counter] = f / TDC_DIVISOR   # accumulator of CKV after divisor to average in TDC
#     avg_counter += 1
#     if avg_counter == AVG_FCKV:
#         avg_counter = 0
        
#     return f

tdc_delay = []
def TDC(tR , t_ckv ,TCKV_accumulator):
    global TDC_res , T0 , TDC_chains , RV_n , RV_k , NUM_ZEROS, AVG_FCKV, TDC_BITS_RESOLUTION, TDC_BITS_RESOLUTION, TDC_MISMATCH
    # tR_Q = int((tR - t_ckv) / TDC_res)  # Diferença de tempo entre a última borda de clock de CKV até a borda de REF. (FIG 2 Time-Domain Modeling of an RF All-Digital PLL)
   
    '''
    Average CKV Clock and normalization to DCO period
    Ref: ADPLL Design for WiMAX (pg. 67)
    '''
    tdc_resolution = random.uniform(TDC_res*(1 + TDC_MISMATCH), TDC_res*(1 - TDC_MISMATCH))
    t = TCKV_accumulator * AVG_FCKV
    T_Q =  int(( TCKV_accumulator * AVG_FCKV)  / tdc_resolution)    # quantização do accumulador em relação a resolução do TDC
    K_TDC = int(( 1 / T_Q ) * 2**TDC_BITS_AVG_TCKV) # ganho do TDC
    '''
    Calculate the diference between edges and quantization to TDC resolution delay
    '''
    difftime = t_ckv - tR
    tR_Q = abs(int((t_ckv - tR) / tdc_resolution)) #quantização da diferença em relação a resolução do TDC
    Binary_error = tR_Q * ( K_TDC / 2**TDC_BITS_AVG_TCKV) * 2**TDC_BITS_RESOLUTION
    tdc_delay.append(T_Q)
    Float_error = (1 / 2**TDC_BITS_RESOLUTION) * Binary_error
    error = 1 - Float_error
    return error


def SDM_modulator(ntw_f):
    global BusSize , U1 , U2 , U3 , C1 , C2 , C3 , AccumulatorSize , out , f_CKV , T0 , OTW_trk , index

    FractionInternal = 2 ** BusSize * ntw_f
    U1[index] = FractionInternal + U1[index - 1]
    U2[index] = U1[index - 1] + U2[index - 1]
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
    out[index] = C1[index - 3] + C2[index - 2] - C2[index - 3] + C3[index - 1] - 2 * C3[index - 2] + C3[index - 3]
    # out[index] =  C1[index] + C2[index] - C2[index-1] + C3[index] - 2*C3[index-1] + C3[index-2]
    # NTW[k - 1] += out[index]
    # otw_prev = OTW_trk + out[index]
    otw_prev = OTW_trk + C1[index]
    # OTW_trk += out[index]
    if otw_prev > 64:
        otw_prev = 64
    elif otw_prev < 0:
        otw_prev = 0
    f_CKV = SET_DCO(OTW_pvt , OTW_acq , otw_prev , OTW_trk_f)
    T0 = 1 / f_CKV


def plot_DCO_signal():
    '''
    plotar sinal de saída do DCO e comparar ao sinal de inicio

    INPUT arguments
    none
    '''
    global len_simulation , dco_init_time , dco_freq_init , f_CKV , T0 , OVERSAMPLE , fs
    fs = OVERSAMPLE * f_CKV
    t = np.arange(0 , 10 * T0 , 1 / fs)
    # x_noise = np.sin(2 * np.pi * 1/(tc + jitter + wander) * t)
    x = np.sin(2 * np.pi * f_CKV * t)
    plt.figure()
    plt.plot(dco_init_time[:len_simulation] / 1e-9 , dco_freq_init[:len_simulation] , label="DCO Inicial")
    plt.plot(t[:len_simulation] / 1e-9 , x[:len_simulation] , label="DCO out")
    plt.grid(visible=True)
    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (V)')
    plt.show()

def plot_phaseNoiseMask():
    # Geração das frequências
    frequencies_1k_to_1M = np.linspace(1000, 1e6, 300)  # De 1k a 1M
    frequencies_2M_to_2M= np.linspace(1e6, 2e6, num=100)  # De 1M a 2M
    frequencies_2M_to_3M = np.linspace(2e6, 3e6, num=100)  # De 2M a 3M
    frequencies_above_3M = np.linspace(3e6, 5e9, num=300)  # Acima de 3M

    # Geração dos níveis de fase noise correspondentes
    phase_noise_1k_to_1M = np.full_like(frequencies_1k_to_1M, -66)
    phase_noise_1M_to_2M = np.full_like(frequencies_2M_to_2M, -98)
    phase_noise_2M_to_3M = np.full_like(frequencies_2M_to_3M, -108)
    phase_noise_above_3M = np.full_like(frequencies_above_3M, -108)

    # Concatenação dos dados
    phase_noise = np.concatenate((phase_noise_1k_to_1M, phase_noise_1M_to_2M, phase_noise_2M_to_3M, phase_noise_above_3M))
    freq = np.concatenate((frequencies_1k_to_1M, frequencies_2M_to_2M, frequencies_2M_to_3M, frequencies_above_3M))

    return phase_noise, freq
    


def plot_histogram_noise(lenght):
    '''
    plotar histograma do ruído Wander e jitter

    INPUT arguments
    lenght     :  Quantidade de pontos aleatórios
    '''
    global Jt_noise , Wt_noise
    jitter_noise = np.random.randn(lenght) * Jt_noise
    wander_noise = np.random.randn(lenght) * Wt_noise
    plt.subplot(121)
    plt.hist(jitter_noise , bins=100 , label="jitter") #"Jitter $\sigma_{\Delta t}$ = 73 fs ") 
    #label="Normal distribution of the Jitter noise")
    plt.legend()
    plt.grid(visible=True)
    plt.subplot(122)
    plt.hist(wander_noise , bins=100 , color="r" , label="Wander")#"Wander $\sigma_{\Delta T}$ = 4 fs ") 
    #label="Normal distribution of the wander noise")
    plt.legend()
    plt.grid(visible=True)
    plt.show()


def diff_freq(f):
    '''
        Calcula a diferença de frequêcia do valor desejado ao atual

        INPUT arguments
        f     :  Valor da frequência atual
        '''
    return abs((F_DESIRED * 2) - f)


def IRR_filter(k , ntw):
    global y1 , y2 , y3 , y_IRR , IRR_coef

    y1[k] = (1 - IRR_coef[0]) * y1[k - 1] + IRR_coef[0] * ntw
    y2[k] = (1 - IRR_coef[1]) * y2[k - 1] + IRR_coef[1] * y1[k]
    y3[k] = (1 - IRR_coef[2]) * y3[k - 1] + IRR_coef[2] * y2[k]
    y_IRR[k] = (1 - IRR_coef[3]) * y_IRR[k - 1] + IRR_coef[3] * y3[k]

    return y_IRR[k]


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


def saveresults(timestop , timediff , fout_n , desv_n , fout_T , desv_T , fout_SDM , desv_SDM , result , dfresult):
    '''
    Save results of simulation in a csv file
    '''
    global csvsaveresults , FCW , IRR , SDM , NOISE , TIME , FREF

    # ['starttime', "stoptime", "duration", "noise enable", "IRR enable", "SDM enable", "FCW", "times FREF",
    #  "Fout nor.", "desv nor.", "Fout total mean", "desv total mean", "Fout SDM mean", "desv SDM mean"]
    locale.setlocale(locale.LC_ALL , 'pt_BR.UTF-8')
    dfresult['stoptime'] = timestop
    dfresult['duration'] = locale.format_string('%.4f' , timediff)
    dfresult['noise enable'] = "Ativado" if NOISE else "Desativado"
    dfresult['IRR enable'] = "Ativado" if IRR else "Desativado"
    dfresult['SDM enable'] = "Ativado" if SDM else "Desativado"
    dfresult['FCW'] = locale.format_string('%.6f' , FCW)
    dfresult['times FREF'] = locale.format_string('%.3f' , TIME)
    dfresult['Fout required'] = locale.format_string('%.4f' , (F_DESIRED * 2) / 1e6)
    dfresult['Fout nor.'] = locale.format_string('%.4f' , fout_n)
    dfresult['desv nor.'] = locale.format_string('%.4f' , desv_n)
    dfresult['Fout total mean'] = locale.format_string('%.4f' , fout_T)
    dfresult['desv total mean'] = locale.format_string('%.4f' , desv_T)
    dfresult['Fout SDM mean'] = locale.format_string('%.4f' , fout_SDM)
    dfresult['desv SDM mean'] = locale.format_string('%.4f' , desv_SDM)

    testeTesting = pd.concat([result , dfresult] , ignore_index=True , axis=0)
    testeTesting.to_csv(csvsaveresults , index=False , sep=';')

def SaveCsvValues(name, x, y):
    '''
    Save results of phase Noise in a csv file
    '''
    global csvsavePhaseresults

    locale.setlocale(locale.LC_ALL , 'pt_BR.UTF-8')
    #phasenoise_locale = [locale.format_string('%.4f', valor, grouping=True) for valor in phasenoise]
    # freq_locale = [locale.format_string('%.4f', valor, grouping=True) for valor in freq]
    # df = pd.DataFrame({'PN': phasenoise_locale, 'freq': freq_locale})
    df = pd.DataFrame({'x': x, 'y': y})
    df.to_csv(name , index=False , sep=';',float_format='%.2f')

'''
        DEFINIÇÕES GERAIS
'''
FCW_I_bits = 8  # Frequency Command Word integer resolution
FCW_F_bits = 24 # Frequency Command Word Fractional resolution
DIVISION_OUTPUT = 2 # Divisor after DCO output

# F0 = 2045e6  # frequência central de ajuste do DCO
FREF = 26e6  # Frequência de referência
# F_DESIRED = 2e9
F_DESIRED = 2.48e9
# FCW = 76.923076927661896  #  Frequency command word
'''
Normalize the FCW integer + fractional in relationship the Nbits resolutions
FCW_FRAC = FREF/2^(Nbits)
'''
##########################################################
fcw_temp = (F_DESIRED * 2) / FREF
fcw_f_res = FREF/2**FCW_F_bits
fCW_f , fCW_i = math.modf(fcw_temp) 
fCW_f = int((fCW_f * FREF) / fcw_f_res)
fCW_f = (fCW_f * fcw_f_res )/ FREF
##########################################################
FCW =   fCW_i + fCW_f  # Frequency command word
FDCO = FREF * FCW  # Frequência desejada na saída do DCO
FREF_edge = 1 / FREF  # tempo de borda de FREF
FDCO_edge = 1 / FDCO  # tempo de borda de F0
NOISE = True
IRR = True
SDM = False
SAVE = False
ENGLISH = False

'''
        BANK CAPACITOR
'''
PVT_NB = 7  # número de bits em PVT mode
ACQ_NB = 7  # número de bits em acquisition mode
TRK_NB_I = 7  # número de bits Trekking integer mode
TRK_NB_F = 5  # número de bits Trekking fractional mode
FR_PVT = 1000e6  # range de frequência em PVT mode
""" FR_ACQ = 80e6  # range de frequência em acquisition mode
FR_TRK_I = 2e6  # range de frequência em Trekking integer mode """
FR_ACQ = 200e6  # range de frequência em acquisition mode
FR_TRK_I = 4e6  # range de frequência em Trekking integer mode
FR_TRK_F = FR_TRK_I / 2 ** TRK_NB_I  # range de frequência em Trekking fractional mode
F0_PVT = 4.8e9 # frequência central do PVT bank
F0_ACQ = 2.44e9 * 2  # frequência central do ACQ bank
F0_TRK = 4.8e9  # frequência central do TRK bank
FREQ_RES_PVT = 0
FREQ_RES_ACQ = 0
FREQ_RES_TRK = 0
FREQ_RES_TRK_F = 0
L = 1e-9  # Indutor utilizado
C0 = 0  # valor de capacitância inicial
MISMATCH_DCO = 0#0.01/100 # 0,01%
pvt_lsb = 0  # valor do LSB em PVT mode
pvt_lsb_sup = 0
pvt_lsb_inf = 0
acq_lsb = 0  # valor do LSB em acquisition mode
acq_lsb_sup = 0
acq_lsb_inf = 0
trk_i_lsb = 0  # valor do LSB em Trekking integer mode
trk_lsb_sup = 0
trk_lsb_inf = 0
trk_f_lsb = 0  # valor do LSB em Trekking fractional mode
desvio_padrao_pvt = 0 
desvio_padrao_acq = 0
desvio_padrao_trk = 0
'''
        TDC
'''
TDC_res = 16e-12  # delay of each  inverter
TDC_chains = 28  # number of inverters
AVG_FCKV = 128  # 128  # number of edges to average actual frequency
NUM_ZEROS = 0
TDC_BITS_RESOLUTION = 16  #Bits resolution of TDC
TDC_BITS_AVG_TCKV = 15    #Bits resolution to average the clock CKV
TDC_DIVISOR = DIVISION_OUTPUT
TDC_MISMATCH = 1/100  #5%
'''
        LOOP FILTER
        
The traditional rule of thumb is to select the offset frequency where (reference
phase noise+TDC phase noise) and DCO phase noise meet as the PLL bandwidth
Besides the bandwidth, for type-II PLL, we have another value which is the damping factor
ζ. For our case, a over-damped type-II PLL (ζ ≥ 1) is preferred since it avoids a significant
gain peaking 7
Calculating the TDC noise and DCO noise, the values are equal ate 70kHz of distance the fv to DCO, thus the 
Loop Filter bandwidth is dimensioned to thhis frequency
'''
Kp_PVT = 2 ** -2
Kp_ACQ = 2 ** -5
Kp_TRK = 2 ** -5
Ki_TRK = 2 ** -12
w_n = np.sqrt(Ki_TRK) * FREF
damping_factor = 0.5 * (Kp_TRK /  np.sqrt(Ki_TRK))
print(f"LOOP FILTER --- Wn={w_n}rad/s and Damping Factor={damping_factor}")
MOD_ARITH = 2 ** 8

'''
        NOISE
'''
noise_floor = -150  # -150  # noise floor [dBc)
L_j = 10 ** (noise_floor / 10)  # noise level
f_desired = F_DESIRED * 2 # F0  # desired frequency
t_required = 1 / f_desired  # period of frequency
Thermal_noise = -111  # 6dB acima do desjado para dobro de freq.  # Up converted Thermal noise with deltaf frequency offset [dBc]
L_w = 10 ** (Thermal_noise / 10)  # noise level
deltaf = 1e6  # offset frequency

j_noise = (t_required / (2 * np.pi)) * np.sqrt(L_j * f_desired)  # Jitter noise standard deviation
W_noise = deltaf / f_desired * np.sqrt(t_required) * np.sqrt(L_w)  # Wander noise standard deviation (including the 1/f noise)
# Converte o número em um decimal
j_decimal = decimal.Decimal(j_noise)
w_decimal = decimal.Decimal(W_noise)
# Arredonda o número com uma precisão de 15 dígitos
j_decimal = j_decimal.quantize(decimal.Decimal('1e-15') , rounding=decimal.ROUND_HALF_EVEN)
w_decimal = w_decimal.quantize(decimal.Decimal('1e-15') , rounding=decimal.ROUND_HALF_EVEN)
# Converte o número arredondado de volta para notação científica
Jt_noise = float('{:e}'.format(j_decimal))
Wt_noise = float('{:e}'.format(w_decimal))
print("jitter noise" , Jt_noise)
print("Wander noise" , Wt_noise)
# Wt_noise = 12e-15  # Wander noise time
# Jt_noise = 111e-15  # jitter noise time

'''
        VARIÁVEIS DE CONTROLE DA SIMULAÇÃO
'''
TIME = 50000  # simulação de X bordas de FREF
OVERSAMPLE = 100  # over sample de frequência para discretizar a frequência do DCO
len_simulation = 6 * OVERSAMPLE  # plotar 6 períodos do DCO

'''
        IRR FILTER
'''
y1 = np.zeros(TIME)
y2 = np.zeros(TIME)
y3 = np.zeros(TIME)
y_IRR = np.zeros(TIME)
IRR_coef = [2 ** -2 , 2 ** -2 , 2 ** -2 , 2 ** -2]
# IRR_coef = [2 ** -2 , 2 ** -1 , 2 ** -1 , 2 ** -1]
f_BW_Irr = IRR_coef[0]/ 2 * np.pi * FREF

'''
        SIGMA DELTA MODULATOR
'''
SDM_DIVISION = 4  # Ratio of division DCO clock
count_div = 0
NumberSamples = TIME * int(FCW / SDM_DIVISION + 1)
BusSize = 5  # bits
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
out = np.zeros(NumberSamples)  # output to the SDM
index = 2  # index to the SDM arrays

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
jitter = [0]  # jitter noise
wander = [0]  # wander noise
error_fractional = np.zeros(TIME)  # fractional error from TDC
phase_error = np.zeros(TIME)  # phase detector
fractional_error_trk = []  # fractional error to plot
fractional_error_trk_IRR = []  # fractional error with IRR aplied to plot

pvt_bank_calib = False
acq_bank_calib = False
trk_bank_calib = False
OTW_pvt = 0  # initial value of pvt bank
OTW_acq = 2 **ACQ_NB /2  -1# initial value of acq bank
OTW_trk = 2 **TRK_NB_I /2  # initial value of trk integer bank
OTW_trk_f = 0  # initial value of trk fractional bank
phase_dif = 0  # phase difference
prev_phase = 0  # new phase difference
count = 0  # counter times when the frequency deviation is within the resolution capacitor bank value
k = 1  # index k
n = 0  # index n
avg_counter = 0  # counter of number f_ckv edges to calculate the average
freq_array = np.zeros(AVG_FCKV)  # array to the values of fckv
freqs = np.zeros(TIME)  # array of different DCO output values
NTW = np.zeros(TIME)  # normalize tuning word
OTW = np.zeros(TIME)  # oscilator tuning word

'''
        Files
'''
csvsaveresults = "resultssimulations.csv"
dfresult = pd.DataFrame()

################################################################################################################
'''
        Main
'''
if __name__ == "__main__":
    starTime = datetime.datetime.now()
    print("iniciando simulação do ADPLL em: " , starTime.strftime("%H:%M:%S") , "\r\nFrequência de referência: " ,
          FREF / 1e6 , "MHz" ,
          "\r\nFrequência desejada: " , F_DESIRED / 1e6 , "MHz" , "\r\nBordas de referência: " , TIME , "vezes \r\nFCW: " ,
          FCW ,
          "\r\nRuido: " , "Ativado" if NOISE else "Desativado" ,
          "\r\nSDM: " , "Ativado" if SDM else "Desativado" ,
          "\r\nFiltro IRR: " , "Ativado" if IRR else "Desativado")
    Init_DCO()
    f_CKV = SET_DCO(OTW_pvt , OTW_acq , OTW_trk , 0)
    T0 = 1 / f_CKV
    fs = OVERSAMPLE * f_CKV
    print("frequência inicial do DCO é: " , f_CKV / 1e6 , "MHz")

    ################  PLOT DO RUIDO ##################################
    #plot_histogram_noise(100000)
    #plot_DCO_signal()

    ################  CRIA CABEÇALHO DO CSV SE NÃO EXISTE ##################################
    if not os.path.exists(csvsaveresults):
        columns = ['starttime' , "stoptime" , "duration" , "noise enable" , "IRR enable" , "SDM enable" , "FCW" ,
                   'Fout required' , "times FREF" ,
                   "Fout nor." , "desv nor." , "Fout total mean" , "desv total mean" , "Fout SDM mean" ,
                   "desv SDM mean"]
        simulationResults = pd.DataFrame(columns=columns)
    else:
        simulationResults = pd.read_csv(csvsaveresults , delimiter=';')
    dfresult['starttime'] = [starTime.strftime("%H:%M:%S")]
    freqmeanall = []
    freqmeanSDM = []
    K = 0
    for k in range(1 , TIME):
        K += FCW  # reference phase accumulator
        RR_k += FCW
        t_R = k * FREF_edge
        U1[index] = 0
        U2[index] = 0
        U3[index] = 0
        while t_CKV[n] < t_R:
            n += 1
            RV_n = n  # variable phase accumulator
            jitter.append(np.random.randn() * Jt_noise)
            wander.append(np.random.randn() * Wt_noise)
            if trk_bank_calib:
                count_div += 1
                if count_div == SDM_DIVISION:
                    index += 1
                    count_div = 0
                    NTW_f = OTW_trk % 1  # NTW[k - 1] % 1  OTW_trk % 1  # get fractional part from NTW
                    if SDM:
                        SDM_modulator(NTW_f)
                        freqmeanSDM.append(f_CKV)
                    freqmeanall.append(f_CKV)
                if NOISE:
                     t_CKV.append(t_CKV[n - 1] + T0 + jitter[n] + wander[n] - jitter[n - 1])  # - TDEV_I
                else:
                    t_CKV.append(t_CKV[n - 1] + T0)
            else:
                # if NOISE:
                #     t_CKV.append(t_CKV[n - 1] + T0 + jitter[n] + wander[n] - jitter[n - 1])
                # else:
                t_CKV.append(t_CKV[n - 1] + T0)
        if trk_bank_calib:
            error_fractional[k] = TDC(t_R*TDC_DIVISOR , t_CKV[n - 1]*TDC_DIVISOR ,1 / np.sum(freq_array))

        RV_k = RV_n  # variable phase accumulator
        erro_esperado = f_CKV / FREF
        erro_esperado = FCW - erro_esperado
        phase_error[k] = (RR_k - RV_k + error_fractional[k])  # Phase detector
        if int(phase_error[k]) < 0:
            phase_error[k] = phase_error[k - 1]
           # print("if int(phase_error[k]) < 0:", phase_error[k])
            pass
        #######################################################################################################################

        ##################### PVT MODE #################################################
        if not pvt_bank_calib:
            NTW[k] = (int(phase_error[k]) * Kp_PVT)  # calcula o novo valor de NTW como inteiro
            OTW[k] = NTW[k] * (FREF / FREQ_RES_PVT)  # gain normalization of TRK mode
            OTW_pvt = OTW[k]
            if diff_freq(f_CKV) <= FREQ_RES_PVT:
                count += 1
                if count == 5:#5:
                    pvt_bank_calib = True
                    print("Frequência na saída do PVT bank: " , f_CKV / 1e6 , "MHz em " , k ,
                          "bordas de FREF e valor do banco ajustado em: " , OTW_pvt)
                    count = 0
            else:
                count = 0
        #######################################################################################################################
        
        ##################### ACQUISITION MODE #########################################
        elif not acq_bank_calib:
            NTW[k] = (int(phase_error[k]) * Kp_ACQ)  # calcula o novo valor de NTW como inteiro
            OTW[k] = NTW[k] * (FREF / FREQ_RES_ACQ) # gain normalization of TRK mode
            OTW_acq = OTW[k]  # ajusta o novo valor de controle dos capacitores do ACQ bank
            if diff_freq(f_CKV) <= FREQ_RES_ACQ:
                count += 1
                if count == 5:
                    acq_bank_calib = True
                    trk_bank_calib = True
                    OTW_acq = OTW[k - 1]
                    print("Frequência na saída do ACQ bank: " , f_CKV / 1e6 , "MHz em " , k ,
                          "bordas de FREF e valor do banco ajustado em: " , OTW_acq)
                    OTW_trk = 0
                    OFFSET_ERROR_ACQ = int(phase_error[k])  # compensão do erro inteiro para ser apenas considerado o erro fracionário
            else:
                count = 0
        #######################################################################################################################

        ##################### TREKING MODE ################################################
        elif trk_bank_calib:
            # NTW[k] = OTW_trk + ((int(phase_error[k] - phase_error[k - 1])) * Kp_TRK)
            if phase_error[k] < (OFFSET_ERROR_ACQ - OFFSET_ERROR_ACQ / 2):
                print(" phase_error[k] < (OFFSET_ERROR_ACQ - OFFSET_ERROR_ACQ / 2): ", phase_error[k], k, OFFSET_ERROR_ACQ)
                # pass
                break
            else:
                phase_error[k] = abs(phase_error[k]) - OFFSET_ERROR_ACQ
            fractional_error_trk.append(phase_error[k])

            if IRR:
                phase_error[k] = IRR_filter(k , phase_error[k])  # aplica o filtro IRR
                fractional_error_trk_IRR.append(phase_error[k])
            NTW[k] = (phase_error[k]) * Kp_TRK - Kp_TRK * (phase_error[k - 1]) + Ki_TRK * (phase_error[k - 1]) + NTW[
                k - 1]  # calcula o novo valor de NTW
            # NTW[k] = (phase_error[k]) * Kp_TRK
            OTW[k] = NTW[k] * (FREF / FREQ_RES_TRK) # gain normalization of TRK mode
            if OTW[k] > 2**TRK_NB_I:
                OTW[k] = 2**TRK_NB_I
            elif OTW[k] < 0:
                OTW[k] = 0
            OTW_trk = OTW[k]  # calcula o novo valor de NTW como inteiro
        #######################################################################################################################
        f_CKV = SET_DCO(OTW_pvt , OTW_acq , OTW_trk , OTW_trk_f)
        last_To = T0
        T0 = 1 / f_CKV
        freqs[k] = f_CKV/DIVISION_OUTPUT  # insere o valor de frequência ajustado no index k
        freqmeanall.append(f_CKV)
        #######################################################################################################################
    
    #####   PRINTS DE INFORMAÇÃO DA SIMULAÇÃO ##############
    print("freq ajustada: " , f_CKV / 1e6 ,
          "MHz E a desejada era de :" , (F_DESIRED * 2) / 1e6 ,
          "MHz diferença de :" , (f_CKV - (F_DESIRED * 2)) / 1e3 , "kHz e valor do banco ajustado em: " , OTW_trk)

    SDMfreq = 0
    SDMDesv = 0
    Last_100_edges_CKV = np.array(freqmeanall[len(freqmeanall)-100:])
    if SDM:
        SDMfreq = np.mean(freqmeanSDM) / 1e6
        SDMDesv = (np.mean(freqmeanSDM) - (F_DESIRED * 2)) / 1e3
        print("freq ajustada considerando a média SDM: " , SDMfreq ,
              "MHz E a desejada era de :" , (F_DESIRED * 2) / 1e6 ,
              "MHz diferença de :" , SDMDesv , "kHz")

    print("freq ajustada considerando a média das ultimas 100 Bordas: " , np.mean(Last_100_edges_CKV) / 1e6 ,
          "MHz E a desejada era de :" , (F_DESIRED * 2) / 1e6 ,
          "MHz diferença de :" , (np.mean(Last_100_edges_CKV) - (F_DESIRED * 2)) / 1e3 , "kHz")
    
    print("ULTIMOS 100 CLOCKS:  Max freq:", np.max(Last_100_edges_CKV), "MHz  MIN freq:", np.min(Last_100_edges_CKV), "MHz")

    print("Overflow TDC: " , NUM_ZEROS)
    print("number repetition of n index:" , n)
    print("ratio of n/k:" , n / k)

    stopTime = datetime.datetime.now()
    diftime = stopTime - starTime
    print("Encerando a simulação em: " , stopTime.strftime("%H:%M:%S"))
    print("Duração da simulação: " , diftime.total_seconds())

    if SAVE:
        saveresults(timestop=stopTime.strftime("%H:%M:%S") , timediff=diftime.total_seconds() , fout_n=f_CKV / 1e6 ,
                    desv_n=(f_CKV - F_DESIRED*DIVISION_OUTPUT) / 1e3 ,
                    fout_T=np.mean(freqmeanall) / 1e6 , desv_T=(np.mean(freqmeanall) - F_DESIRED*2) / 1e3 ,
                    fout_SDM=SDMfreq , desv_SDM=SDMDesv , result=simulationResults , dfresult=dfresult)
    ################################################################################################################
    
    ################ ERRO DE FASE EM RAD/S #################
    # print("Calculando o erro de fase em rad/s")
    # tckv = np.array(t_CKV[len(t_CKV) - 500000:])
    # phase = np.zeros(len(tckv))
    # tref = 1 / F_DESIRED
    # for i in range(len(tckv) - 2):
    #     diff = tref - ((tckv[i + 1] - tckv[i]) * DIVISION_OUTPUT)
    #     phase[i] = diff * 2*np.pi * F_DESIRED

    # print("Max error: ",np.max(phase), " rad/s")

    print("Calculando o erro de fase em rad/s")
    tckv = np.array(t_CKV[len(t_CKV) - 500000:])
    phase = np.zeros(len(tckv))
    tref = 1 / (F_DESIRED * DIVISION_OUTPUT)
    for i in range(len(tckv) - 2):
        diff = tref - ((tckv[i + 1] - tckv[i]) )
        phase[i] = diff * 2*np.pi * F_DESIRED

    print("Max error: ",np.max(phase), " rad/s")
    ########################################################

    ###############   PLOT DOS VALORES DE FREQUÊNCIA E/OU ERROS ######################
    if ENGLISH:
        label1 = "Reference Clock Cycles"
        label2 = 'Frequency out of DCO (Hz)'
    else:
        label1 = "Ciclos de clock de referência"
        label2 = 'Frequência de saída do DCO (Hz)'
    plt.figure()
    plt.plot(freqs[1:500]/1e9, '-r', label="DCO")
    # plt.plot(np.arange(1, TIME, 1), OTW[1:TIME], '-b')
    # # plt.plot(np.arange(0, len(fractional_error_trk), 1), fractional_error_trk, '.')
    # # plt.stem(np.arange(0, len(fractional_error_trk), 1), fractional_error_trk, linefmt='r', markerfmt='.', basefmt="-b")
    # # plt.plot(np.arange(0, len(fractional_error_trk_IRR), 1), fractional_error_trk_IRR, '-b')
    plt.legend()
    plt.xlabel(label1)
    plt.ylabel(label2)
    plt.grid(visible=True)
    # plt.show()
    ##################################################################################

    plt.figure()
    plt.plot(tdc_delay, '-r', label="TDC")
    # plt.plot(np.arange(1, TIME, 1), OTW[1:TIME], '-b')
    # # plt.plot(np.arange(0, len(fractional_error_trk), 1), fractional_error_trk, '.')
    # # plt.stem(np.arange(0, len(fractional_error_trk), 1), fractional_error_trk, linefmt='r', markerfmt='.', basefmt="-b")
    # # plt.plot(np.arange(0, len(fractional_error_trk_IRR), 1), fractional_error_trk_IRR, '-b')
    plt.legend()
    plt.xlabel(label1)
    plt.ylabel('TDC delays')
    plt.grid(visible=True)

    #############  SALVA EM UM CSV OS DADOS DO ARRAY DE phase #########################
    print("Save CSV")
    df = pd.DataFrame([phase])
    nome_arquivo_csv = 'phaseErrorinRad.csv'
    # df.to_csv(nome_arquivo_csv , index=False , header=False)
    ##################################################################################

    ##############  APLICA UM FILTRO NOS VALORES DE PHASE ############################
    b = 1
    a = np.array([1 , -1])
    x = signal.lfilter([b], a,  phase - np.mean(phase))
    ##################################################################################

    ###############  CALCULA O PHASE NOISE DO DCO ####################################
    print("cálculo da PSD")
    print(len(x))
    Xdb_o , f = fun_calc_psd(x , F_DESIRED * DIVISION_OUTPUT , 100e3 , 1e3)
    mask_phase_noise, freq  = plot_phaseNoiseMask() # obter a mascara de phase noise
    marker = 1e6  # Substitua pelo valor específico de frequência desejado
    indice = np.where(f == marker)[0][0]
    marker_dB = Xdb_o[indice]
    SaveCsvValues("phasenoise_PSD_without_IRR.csv",x=f, y=Xdb_o)
    plt.figure()
    if ENGLISH:
        label1 = "Phase Noise"
    else:
        label1 = "Ruído de fase"
    plt.semilogx(f , Xdb_o , label=label1)
    plt.scatter(marker, marker_dB, color='black', marker='o', label=f'{marker_dB:.2f} dBc/Hz  @1 MHz')
    plt.semilogx(freq, mask_phase_noise, label='Phase Noise MASK')
    plt.grid(visible=True)
    plt.legend()
    plt.yticks([-160, -150, -140, -130, -120, -110, -100, -90, -80, -70])
    plt.xlabel('Freq. (Hz)')
    plt.ylabel(label1 + ' [dBc/Hz]')
    # plt.show()
    ##################################################################################

    stopTime = datetime.datetime.now()
    diftime = stopTime - starTime
    print("Encerando a simulação em: " , stopTime.strftime("%H:%M:%S"))
    print("Duração da simulação: " , diftime.total_seconds())
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