"""
Created on maio 07 17:21:52 2023

@author: Ânderson Felipe Weschenfelder
"""
import numpy as np
import matplotlib.pyplot as plt

FCWf = 0.3438
error = 1/12

cycles = 19

out1 = np.zeros(cycles)
in1 = np.zeros(cycles)
e1 = np.zeros(cycles)

out2 = np.zeros(cycles)
in2 = np.zeros(cycles)
e2 = np.zeros(cycles)

outf = np.zeros(cycles)
cout = np.zeros(cycles)
in1[0] = FCWf
in2[0] = 0
media1 = 0
media2 = 0
mediaf = 0
for n in range(1, cycles):
    in1[n] += in1[n - 1]
    out1[n] = in1[n] + e1[n - 1] - error
    if out1[n] > 0:
        out1[n] = 1
    else:
        out1[n] = 0
    e1[n] = in1[n] + e1[n - 1] - out1[n]

    in2[n] += e1[n]
    out2[n] = in2[n] + e2[n - 1] - error
    if out2[n] > 0:
        out2[n] = 1
    else:
        out2[n] = 0
    e2[n] = in2[n] + e2[n - 1] - out2[n]

    cout[n] = FCWf - (e2[n] - e2[n - 1]) ** 2
    outf[n] = out1[n] + out2[n - 1]
    media1 += out1[n]
    media2 += out2[n]
    mediaf += outf[n]

media1 = media1 / (cycles - 1)
media2 = media2 / (cycles - 1)
mediaf = mediaf / (cycles - 1)

print(out1)
print(media1)
print(media2)
print(mediaf)
print(cout.sum() / (cout.size - 1))
