"""
Created on abril 26 18:02:04 2023

@author: Ânderson Felipe Weschenfelder
"""

import numpy as np
import matplotlib.pyplot as plt


class C_LSB:
    def __init__(self,fc, nb, fr, L):
        self.fc = fc
        self.nb = nb
        self.fr = fr
        self.fmin,  self.fmax, self.cmin, self.cmax, self.lsb = deltaC_calc(fc, nb, fr, L)
def deltaC_calc (fc, nb, fr, L):
    fmin = fc - fr/2
    fmax = fc + fr/2
    c_max = 1 / (L * (2 * np.pi * fmin) ** 2)
    c_min = 1 / (L * (2 * np.pi * fmax) ** 2)
    C_LSB = (c_max - c_min) / 2 ** nb

    # print("fmin: ", fmin, " fmax: ", fmax, "cmax: ", c_max, " cmim: ", c_min, " LSB: ", C_LSB)
    return fmin, fmax, c_min, c_max, C_LSB

def DCO (pvt_bin,acq_bin, trk_i_bin, trk_f_bin):
    pvt = int(pvt_bin, 2)
    acq = int(acq_bin, 2)
    trk_i = int(trk_i_bin, 2)
    trk_f = int(trk_f_bin, 2)

    C0 = 4.809222e-12
    L = 1e-9
    FR_PVT = 500e6
    FR_ACQ = 100e6
    FR_TRK = 2e6
    F0 = 2045e6

    PVT_NB = 8
    ACQ_NB = 8
    TRK_NB = 6

    pvt_lsb = C_LSB(F0, PVT_NB, FR_PVT, L)
    print("fmin: ", pvt_lsb.fmin, " fmax: ", pvt_lsb.fmax, "cmax: ", pvt_lsb.cmax, " cmim: ", pvt_lsb.cmin, " LSB: ", pvt_lsb.lsb)
    acq_lsb = deltaC_calc(F0, ACQ_NB, FR_ACQ, L)
    trk_i_lsb = deltaC_calc(F0,TRK_NB , FR_TRK, L)

    fc = 1/(2*np.pi*np.sqrt(L*(C0 + pvt * pvt_lsb + acq * acq_lsb + trk_i * trk_i_lsb)))
    tc = 1/fc
    print("freq: ",fc)
    t = np.arange(0, 10 * tc, tc/20)

    x = np.sin(2 * np.pi * fc * t)
    len_simulation = 6*20

    plt.figure()
    plt.plot(t[:len_simulation] / 1e-9, x[:len_simulation], label="V(t)")
    plt.grid(visible=True)
    plt.legend()
    plt.show()




DCO("10010000","10000000","100000","101")
