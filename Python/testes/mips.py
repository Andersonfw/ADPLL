import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
from scipy import signal
import pandas as pd


class struction_set:
    def __init__(self):
        self.opcode = 0
        self.op1 = 0
        self.op2 = 0
        self.op3 = 0
        self.temp1 = 0
        self.temp2 = 0
        self.temp3 = 0
        self.valida = 0

def get_instruction():

    return 0

def read_instruction():

    return 0

def ULA(operation = None):

    if operation == ADD:

        return
    
    if operation == ADDI:

        return
    
    if operation == SUB:

        return
    
    if operation == SUBI:

        return
    
    if operation == BEQ:

        return
    
    if operation == J:

        return
    
    return 0

def memory_acess():

    return 0

def write_memory():

    return 0


N_REGISTRADORES = 32        # nº de registradores

REGISTERS = np.zeros(N_REGISTRADORES)



'''
        Main
'''
if __name__ == "__main__":

    print("Simulador MIPS")