import numpy as np
import matplotlib.pyplot as plt

# Geração das frequências
frequencies_1k_to_1M = np.logspace(3, 6, num=300)  # De 1k a 1M
frequencies_1M_to_10M = np.logspace(6, 7, num=200)  # De 1M a 10M
frequencies_above_10M = np.logspace(7, 9, num=300)  # Acima de 10M

# Geração dos níveis de fase noise correspondentes
phase_noise_1k_to_1M = np.full_like(frequencies_1k_to_1M, -90)
phase_noise_1M_to_10M = np.full_like(frequencies_1M_to_10M, -130)
phase_noise_above_10M = np.full_like(frequencies_above_10M, -140)

# Concatenação dos dados
frequencies = np.concatenate((frequencies_1k_to_1M, frequencies_1M_to_10M, frequencies_above_10M))
phase_noise = np.concatenate((phase_noise_1k_to_1M, phase_noise_1M_to_10M, phase_noise_above_10M))

# Plot
plt.semilogx(frequencies, phase_noise, label='Phase Noise')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase Noise (dBc/Hz)')
plt.title('Phase Noise Mask')
plt.grid(True)
plt.legend()

plt.ylim(-160, -60)
plt.show()
