"""
Created on setembro 12 17:46:23 2023

@author: Ânderson Felipe Weschenfelder
"""
import numpy as np


DELTA_F = 70e3
TDC_RES = 16e-12
FRES_DCO = 32e3
FREF = 26e6
TREF = 1/FREF
FDCO = 4.8e9/2
TDCO = int(1/FDCO / TDC_RES) * TDC_RES



noise_DCO_time = (1/12) * ((FRES_DCO/DELTA_F)**2) * TREF * (np.sinc(DELTA_F/FREF)**2)
noise_DCO_DB = 10*np.log10(noise_DCO_time)

noise_TDC_time = ( ( (2*np.pi)**2) / 12) * (TDC_RES/TDCO)**2 * TREF
noise_TDC_DB = 10*np.log10(noise_TDC_time)


print("TDC   ---- time:", noise_TDC_time, "  ", noise_TDC_DB, "DB")
print("DCO   ---- time:", noise_DCO_time, "  ", noise_DCO_DB, "DB")
