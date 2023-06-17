def int_to_complemento2(valor, num_bits):
    # Verificar se o valor está dentro do intervalo representável
    limite_superior = 2 ** (num_bits - 1) - 1
    limite_inferior = -2 ** (num_bits - 1)
    # if valor < limite_inferior or valor > limite_superior:
    #     raise ValueError("Valor fora do intervalo representável.")

    # Converter para representação binária
    binario = bin(valor & int("1" * num_bits, 2))[2:]

    # Preencher com zeros à esquerda se necessário
    binario = binario.zfill(num_bits)

    return binario


def complemento2_to_int(binario):
    # Verificar se o número é negativo (bit mais significativo é 1)
    if binario[0] == "1":
        # Aplicar complemento de 2 invertendo todos os bits
        invertido = "".join("1" if bit == "0" else "0" for bit in binario)

        # Adicionar 1 ao resultado
        complemento2 = bin(int(invertido, 2) + 1)

        # Converter para valor inteiro negativo
        valor = int(complemento2, 2) * -1
    else:
        # Converter para valor inteiro positivo
        valor = int(binario, 2)

    return valor


valor = 228
num_bits = 8

complemento2 = complemento2_to_int(int_to_complemento2(valor, num_bits))
# complemento2 = complemento2_to_int(complemento2)
print("Valor em complemento de 2:", complemento2)