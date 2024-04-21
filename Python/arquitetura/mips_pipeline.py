import pandas as pd
import matplotlib.pyplot as plt
import os


FILE = 'code.txt'
valid_opcodes = {"add", "addi", "sub", "subi", "halt", "noop", "beq"}
variables = dict()

class Instruction:
    def __init__(self,num, label=None, opcode=None, op1=None, op2=None, op3=None):
        self.label = label if label is not None else 0
        self.opcode = opcode if opcode is not None else 0
        if op1 is not None:
            self.op1 = int(op1) if  op1.startswith("-") or op1.isdigit() else op1
        if op2 is not None:
            self.op2 = int(op2) if op2.isdigit() else op2
        if op3 is not None:
            self.op3 = int(op3) if op3.isdigit() else op3
        self.num = num
        self.UlaResult = None
        self.DataMem = None
        


class Pipeline:
    def __init__(self):
        self.registers = [0] * 32  # Banco de registradores R0 a R31
        self.memory = []  # Memória de programa
        self.pc = 0  # Contador de programa
        self.pipeline_registers = [None] * 5  # Registradores do pipeline (estágios)
        self.pipeline_registers_valid = [False] * 5  # Validação dos registradores do pipeline
        self.stop = False
        # self.pipeline

    def load_program(self, filename):
        halt_detect = False
        count = 0
        print("-----------------------------------------------------------")
        print("Reading File: ", filename)
        with open(filename, 'r') as file:
            for line in file:
                parts = line.strip().split(',')
                print(parts)

                if parts[0] not in valid_opcodes:
                     if parts[1] == 'halt':
                         halt_detect = True
                         variables[parts[1]] = count
                     if halt_detect == True and parts[1] != 'halt':
                         variables[parts[0]] = count #int(parts[2])
                     else:
                         variables[parts[0]] = count #int(parts[2])
                        
                instruction = Instruction(count, *parts)
                self.memory.append(instruction)
                count += 1
        print("End File")
        print("-----------------------------------------------------------")

    def fetch(self):
        if self.pc < (variables['halt'] + 1):
            self.pipeline_registers[0] = self.memory[self.pc]
            self.pipeline_registers_valid[0] = True
            instruction = self.memory[self.pc]
            self.pc = instruction.num + 1

    def decode(self):
        if self.pipeline_registers_valid[0]:
            self.pipeline_registers[1] = self.pipeline_registers[0]
            self.pipeline_registers_valid[1] = True
        else:
            self.pipeline_registers_valid[1] = False

    def execute(self):
        if self.pipeline_registers_valid[1]:
            instruction = self.pipeline_registers[1]
            opcode = instruction.opcode
            self.pipeline_registers[2] = self.pipeline_registers[1]
            self.pipeline_registers_valid[2] = True
            if opcode == 'add':
                self.UlaResult = self.registers[int(instruction.op2)] + self.registers[int(instruction.op3)]
                # self.registers[int(instruction.op1)] = self.registers[int(instruction.op2)] + self.registers[int(instruction.op3)]
            elif opcode == 'addi':
               # self.registers[int(instruction.op2)] = self.registers[int(instruction.op1)] + meu_dicionario[instruction.op3]#int(instruction.op3)
                self.UlaResult = self.registers[int(instruction.op1)] + variables[instruction.op3]
            elif opcode == 'sub':
                # self.registers[int(instruction.op1)] = self.registers[int(instruction.op2)] - self.registers[int(instruction.op3)]
                self.UlaResult = self.registers[int(instruction.op2)] - self.registers[int(instruction.op3)]
            elif opcode == 'subi':
                # self.registers[int(instruction.op1)] = self.registers[int(instruction.op2)] - int(instruction.op3)
                self.UlaResult = self.registers[int(instruction.op1)] - variables[instruction.op3]
            elif opcode == 'beq':
                self.UlaResult = self.registers[int(instruction.op2)] - self.registers[int(instruction.op1)]
                if self.UlaResult == 0:
                    self.pc = variables[instruction.op3]
                    # Implementação do branch equal
            elif opcode == 'j':
                pass  # Implementação do jump
            elif opcode == 'noop':
                pass
            elif opcode == 'halt':
                print("halt")
                self.stop = True
            else:
                self.pipeline_registers_valid[2] = False

    def memory_access(self):
        if self.pipeline_registers_valid[2]:
            instruction = self.pipeline_registers[2]
            self.pipeline_registers_valid[3] = True
            self.pipeline_registers[3] = self.pipeline_registers[2]
            opcode = instruction.opcode
            if opcode == 'addi' or opcode == 'subi':
                self.DataMem = self.memory[self.UlaResult]

            pass  # Implementação do acesso à memória
        else:
            self.pipeline_registers_valid[3] = False

    def write_back(self):
        if self.pipeline_registers_valid[3]:
            instruction = self.pipeline_registers[3]
            self.pipeline_registers_valid[4] = True
            opcode = instruction.opcode
            if opcode == 'addi' or opcode == 'subi':
                self.registers[int(instruction.op2)] = self.DataMem.op1
                print("Registers: R0:", self.registers[0], " R1:", self.registers[1], " R2:", self.registers[2])
            elif opcode == 'add'or opcode == 'sub':
                self.registers[int(instruction.op3)] = self.UlaResult
                print("Registers: R0:", self.registers[0], " R1:", self.registers[1], " R2:", self.registers[2])
            pass  # Implementação da escrita de volta
        else:
            self.pipeline_registers_valid[4] = False
            


    def run(self):
        counter_clock = 0
        while True:
            self.write_back()
            self.memory_access()
            self.execute()
            self.decode()
            self.fetch()
            pip1.append(self.pipeline_registers_valid[0])
            pip2.append(self.pipeline_registers_valid[1])
            pip3.append(self.pipeline_registers_valid[2])
            pip4.append(self.pipeline_registers_valid[3])
            pip5.append(self.pipeline_registers_valid[4])
            clock.append(counter_clock)
            counter_clock +=1
            # if self.pc >= meu_dicionario['halt'] and all(not v for v in self.pipeline_registers_valid):
            if self.stop == True:
                break



pip1 = []
pip2 = []
pip3 = []
pip4 = []
pip5 = []
clock = []

pipeline_simulator = Pipeline()
pipeline_simulator.load_program(FILE)
pipeline_simulator.run()
plt.figure()
plt.plot(clock, pip1, label="Fetch")
plt.plot(clock, pip2, label="Decode")
plt.plot(clock, pip3, label="Execute")
plt.plot(clock, pip4, label="Data Memory")
plt.plot(clock, pip5, label="Write Back")
plt.legend()
plt.show(), 