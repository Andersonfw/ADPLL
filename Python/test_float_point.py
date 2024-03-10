import struct

def decimal_para_float24(decimal):
    # Converte o número decimal para um float de 32 bits
    float32 = struct.pack('f', decimal)
    # Extrai os bits relevantes
    bits = struct.unpack('I', float32)[0]  # I representa um inteiro não assinado de 32 bits
    # Ajusta os bits para representar um float de 24 bits
    bits = (bits >> 8) & 0xFFFFFF  # Desloca 8 bits para a direita e mascara para manter apenas os 24 bits menos significativos
    # Converte os bits de volta para bytes
    float24 = struct.pack('I', bits)
   # return struct.pack('I', bits)
    return float24

# Exemplo de uso
numero_decimal = 3.14159
float24 = decimal_para_float24(numero_decimal)
print("Float de 24 bits:", repr(float24))
print("Float de 24 bits em hexadecimal:", float24.hex())


def decimal_para_float_binario(decimal, num_bits):
    if num_bits not in [32, 24, 16]:
        print("Erro: O número de bits deve ser 32, 25 ou 16.")
        return None
    formato = {32: 'f', 24: 'e', 16: 'e'}
    try:
        # Converte o número decimal para float de 32 bits
        float32 = struct.pack(formato[num_bits], decimal)
        # Converte o float para binário
        float_binario = ''.join(f'{byte:08b}' for byte in float32)
        print(f"Número decimal {decimal} convertido para {num_bits} bits em binário float: {float_binario}")
        return float_binario
    except struct.error as e:
        print(f"Erro: Não foi possível converter o número {decimal} para {num_bits} bits em binário float.")
        return None

def float_binario_para_decimal(float_binario, num_bits):
    if num_bits not in [32, 24, 16]:
        print("Erro: O número de bits deve ser 32, 25 ou 16.")
        return None
    formato = {32: 'f', 24: 'e', 16: 'e'}
    try:
        # Converte o binário de volta para float
        float_bytes = bytes(int(float_binario[i:i+8], 2) for i in range(0, len(float_binario), 8))
        float_decimal = struct.unpack(formato[num_bits], float_bytes)[0]
        print(f"Binário float {float_binario} convertido de volta para decimal: {float_decimal}")
        return float_decimal
    except struct.error as e:
        print(f"Erro: Não foi possível converter o binário {float_binario} para decimal.")
        return None

# Exemplo de uso
numero_decimal = 3.14159
num_bits = 32
float_binario = decimal_para_float_binario(numero_decimal, num_bits)
float_decimal_novo = float_binario_para_decimal(float_binario, num_bits)
