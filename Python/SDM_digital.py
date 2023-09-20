# """
# Created on maio 07 17:21:52 2023
#
# @author: Ânderson Felipe Weschenfelder
# """
# import numpy as np
# import matplotlib.pyplot as plt
#
# FCWf = 0.3438
# error = 1/12
#
# cycles = 19
#
# out1 = np.zeros(cycles)
# in1 = np.zeros(cycles)
# e1 = np.zeros(cycles)
#
# out2 = np.zeros(cycles)
# in2 = np.zeros(cycles)
# e2 = np.zeros(cycles)
#
# outf = np.zeros(cycles)
# cout = np.zeros(cycles)
# in1[0] = FCWf
# in2[0] = 0
# media1 = 0
# media2 = 0
# mediaf = 0
# for n in range(1, cycles):
#     in1[n] += in1[n - 1]
#     out1[n] = in1[n] + e1[n - 1] - error
#     if out1[n] > 0:
#         out1[n] = 1
#     else:
#         out1[n] = 0
#     e1[n] = in1[n] + e1[n - 1] - out1[n]
#
#     in2[n] += e1[n]
#     out2[n] = in2[n] + e2[n - 1] - error
#     if out2[n] > 0:
#         out2[n] = 1
#     else:
#         out2[n] = 0
#     e2[n] = in2[n] + e2[n - 1] - out2[n]
#
#     cout[n] = FCWf - (e2[n] - e2[n - 1]) ** 2
#     outf[n] = out1[n] + out2[n - 1]
#     media1 += out1[n]
#     media2 += out2[n]
#     mediaf += outf[n]
#
# media1 = media1 / (cycles - 1)
# media2 = media2 / (cycles - 1)
# mediaf = mediaf / (cycles - 1)
#
# print(out1)
# print(media1)
# print(media2)
# print(mediaf)
# print(cout.sum() / (cout.size - 1))

import numpy as np
import matplotlib.pyplot as plt
'''
%MASH ((Multi Stage Noise Shaping)111 SDM
%this is a simulation of a MASH111 SDM where we look at the noise
%generated after a lengthy simulation time. The inputs are the fraction we
%are trying to synthesize, and the number of bits for the accumulators.
%The values in the accumulators are all positive, and the fraction is only
%valid from 0 to 1
%This has no dither it is just the basic modulator
%The setup is:
%                 --------> + ----------> + ----> Output
%                |         diff         diff^2
%                |          |             |
%              carry1     carry2        carry3
%Fraction=>Accumulator1=>Accumulator2=>Accumulator3
%Carry1 from accumulator 1 is added to the dirivative of carry2 and the
%2nd derivative of carry3 to make the output.
%1/29/2013 fixed the 3 stage differentiation
%1/29/2013 added better plot and description
'''

# NumberSamples=2^16;
NumberSamples = 50
BusSize = 21  # bits
Fraction = 0.0205  # usable 0 to 1
FractionInternal = 2**BusSize * Fraction
AccumulatorBits = 21  # bits
AccumulatorSize = 2**AccumulatorBits

C1 = np.zeros(NumberSamples)    # Carry out of the first accumulator
C2 = np.zeros(NumberSamples)    # Carry out of the 2nd accumulator
C3 = np.zeros(NumberSamples)    # Carry out of the 3nd accumulator
U1 = np.zeros(NumberSamples)    # output of the 1st accum
U2 = np.zeros(NumberSamples)    # output of the 2nd accum
U3 = np.zeros(NumberSamples)    # output of the 3rd accum
Yout1 = np.zeros(NumberSamples) # output to the divider for 1 stage SDM
Yout2 = np.zeros(NumberSamples) # output to the divider for 2 stage SDM
Yout3 = np.zeros(NumberSamples) # output to the divider for 3 stage SDM
out = np.zeros(NumberSamples)

for index in range(0, NumberSamples):
    U1[index] = FractionInternal + U1[index-1]
    U2[index] = U1[index - 1] + U2[index - 1]
    # U2[index] = U1[index-1] + U2[index-1]
    U3[index] = U2[index-1] + U3[index-1]
    if U1[index] > AccumulatorSize:
        C1[index] = 1  # carry 1
        U1[index] -= AccumulatorSize
    if U2[index] > AccumulatorSize:
        C2[index] = 1  # carry 2
        U2[index] -= AccumulatorSize
    if U3[index] > AccumulatorSize:
        C3[index] = 1  # carry 3
        U3[index] -= AccumulatorSize

    # The output is the overflow from acc 1, plus the diff of the overflow from
    # acc 2, plus the 2nd derivative of the overflow from acc 3
    Yout3[index] = C1[index] + C2[index] - C2[index-1] + C3[index] - 2*C3[index-1] + C3[index-2]  # output to the divider - 3 stages
    # Yout3[index] = C1[index - 2] + C2[index - 1] - C2[index - 2] + C3[index] - 2 * C3[index - 1] + C3[index - 2]
    Yout2[index] = C1[index] + C2[index] - C2[index-1]  # output to the divider - 2 stages
    Yout1[index] = C1[index]  # output to the divider - 1 stage
    out[index] = C1[index - 3] + C2[index - 2] - C2[index - 3] + C3[index - 1] - 2 * C3[index - 2] + C3[index - 3]
MeanFrac = np.mean(Yout2)
Meanout = np.mean(out)
print(f"\nMeanFracMASH = {MeanFrac:.4f}\n")
print(f"\nMeanout = {Meanout:.4f}\n")
# note how the close in noise improves as the order increases
fig1, ax1 = plt.subplots()
SignalFreq1 = 20*np.log10(np.abs(np.fft.fft(Yout1)))
ax1.plot(np.fft.fftshift(SignalFreq1) - np.max(SignalFreq1), 'g')
SignalFreq2 = 20*np.log10(np.abs(np.fft.fft(Yout2)))
SignalFreq3 = 20*np.log10(np.abs(np.fft.fft(Yout3)))
ax1.plot(np.fft.fftshift(SignalFreq2) - np.max(SignalFreq2), 'r')
ax1.plot(np.fft.fftshift(SignalFreq3) - np.max(SignalFreq3), 'b')
ax1.legend(['1 stage', '2 stage', '3 stage'])
ax1.set_title('MASH111 SDM Noise')
ax1.set_ylim([-150, 0])
ax1.set_xlim([0, NumberSamples])
ax1.grid(True)

fig2, ax2 = plt.subplots()
# ax2.plot(out, 'r')
# ax2.plot(Yout3, 'b')
ax2.plot(C2, 'r')
ax2.plot(C1, 'g')
ax2.plot(C2, 'r')
ax2.plot(C3, 'b')
ax2.set_title('Yout3')
ax2.grid(True)

plt.show()
