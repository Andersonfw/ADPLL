
import numpy as np
def SDM(word_lenght, fractional_value):
    BusSize = word_lenght  # bits
    Fraction = fractional_value  # usable 0 to 1
    FractionInternal = 2 ** BusSize * Fraction
    AccumulatorBits = word_lenght  # bits
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
    out = np.zeros(NumberSamples)

    for index in range(0, NumberSamples):
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

        # The output is the overflow from acc 1, plus the diff of the overflow from
        # acc 2, plus the 2nd derivative of the overflow from acc 3
        Yout3[index] = C1[index] + C2[index] - C2[index - 1] + C3[index] - 2 * C3[index - 1] + C3[
            index - 2]  # output to the divider - 3 stages
        # Yout3[index] = C1[index - 2] + C2[index - 1] - C2[index - 2] + C3[index] - 2 * C3[index - 1] + C3[index - 2]
        Yout2[index] = C1[index] + C2[index] - C2[index - 1]  # output to the divider - 2 stages
        Yout1[index] = C1[index]  # output to the divider - 1 stage
        out[index] = C1[index - 3] + C2[index - 2] - C2[index - 3] + C3[index - 1] - 2 * C3[index - 2] + C3[index - 3]