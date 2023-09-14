"""
Created on setembro 12 17:46:23 2023

@author: Ânderson Felipe Weschenfelder
"""
import numpy as np


FREF = 26e6
TREF = 1/FREF
FDCO = 2.4e9
TDCO = 1/FDCO
DELTA_F = 3.5e6
TDC_RES = 12e-12
FRES_DCO = 31e3


noise_DCO_time = (1/12) * ((FRES_DCO/DELTA_F)**2) * TREF * (np.sinc(DELTA_F/FREF)**2)
noise_DCO_DB = 10*np.log10(noise_DCO_time)

noise_TDC_time = ( ( (2*np.pi)**2) / 12) * (TDC_RES/TDCO)**2 * TREF
noise_TDC_DB = 10*np.log10(noise_TDC_time)


print("TDC   ---- time:", noise_TDC_time, "  ", noise_TDC_DB, "DB")
print("DCO   ---- time:", noise_DCO_time, "  ", noise_DCO_DB, "DB")
