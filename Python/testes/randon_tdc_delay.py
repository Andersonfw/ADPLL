import numpy as np
import matplotlib.pyplot as plt

# Valor original da variável (em ps)
valor_original = 14e-18

# Variação permitida (+/- 5%)
variacao_percentual = 0.05

# Gerar 10000 amostras aleatórias com distribuição normal
desvio_padrao = valor_original * variacao_percentual
amostras = np.random.normal(loc=0, scale=desvio_padrao, size=10000)
vetor_original = np.full(10000, valor_original)
vetor_variado = vetor_original + amostras

# Plotar o histograma
plt.hist(vetor_variado, bins=30, density=True, alpha=0.7, color='blue')
plt.title('Histograma da Variação da Variável')
plt.xlabel('Valor (ps)')
plt.ylabel('Densidade de Probabilidade')
plt.grid(True)
plt.show()
