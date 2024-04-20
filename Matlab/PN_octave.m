# Define filename (replace with your actual path)
nome_arquivo_csv = 'C:\Users\ander\OneDrive - Associacao Antonio Vieira\UNISINOS\TCC\Python\phaseErrorinRad.csv';

# Frequencies
frequencies_1k_to_3_5M = linspace(1000, 1e6, 300);
frequencies_3_5M_to_10M = linspace(1e6, 2e6, 200);
frequencies_above_10M = linspace(2e6, 5e9, 300);

# Phase noise
phase_noise_1k_to_1M = -66 * ones(size(frequencies_1k_to_3_5M));
phase_noise_1M_to_10M = -98 * ones(size(frequencies_3_5M_to_10M));
phase_noise_above_10M = -108 * ones(size(frequencies_above_10M));

# Combine frequencies and phase noise
freq = [frequencies_1k_to_3_5M, frequencies_3_5M_to_10M, frequencies_above_10M];
phase_noise = [phase_noise_1k_to_1M, phase_noise_1M_to_10M, phase_noise_above_10M];

# Read data from CSV
x = dlmread(nome_arquivo_csv);

# Sampling frequency
Fs = 26e6;
Ts = 1/Fs;

# Filter phase (remove mean)
phase = filter(1, [1 -1], x - mean(x));

# Calculate statistics
media = mean(x);
max_phase = max(phase);
min_phase = min(phase);

# Number of segments and window length
num_segments = 16;
window_length = floor(length(phase) / num_segments);

# Calculate PSD using Welch's method
[PSDphase, f] = pwelch(phase, window_length, [], [], Fs);
PSDphase = PSDphase .* (sinc(f/Fs)).^2;  % Correct for ZOH
PSDphase = 10 * log10(PSDphase);

# Output filename (replace with your desired path)
filename = 'C:\Users\ander\OneDrive - Associacao Antonio Vieira\UNISINOS\TCC\Python\output.csv';

# Save filtered phase to CSV
#dlmwrite(filename, phase);

# Calculate PSD with custom function (assuming fun_calc_psd is defined elsewhere)
[custom_PSDphase, f1] = PSD_octave(phase, 2.4e9, 50e3, 1e3);
# Find maximum PSD
psd_max = max(PSDphase);

# Plot results
figure;
hold on;
##h(1) = semilogx(f, PSDphase, 'b');
h(1) = semilogx(f1, custom_PSDphase, 'b');
#h(2) = semilogx(freq, phase_noise, 'r');
grid on;
set(gca, 'fontsize', 15, 'fontweight', 'bold');
xlabel('Frequency [Hz]');
ylabel('Phase Noise [dBc/Hz]');
legend('Matlab Model', 'Phase Noise Mask');
set(h(1), 'LineWidth', 2);
hold off;



# Additional plot (optional)
## x = blackman(5265);
## figure;
## plot(x, 'b');
## grid on;

