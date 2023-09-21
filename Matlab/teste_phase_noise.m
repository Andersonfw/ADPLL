
nome_arquivo_csv = 'C:\Users\ander\OneDrive - Associacao Antonio Vieira\UNISINOS\TCC\Python/phaseErrorinRad.csv';

% Frequências de 1k a 3,5M
frequencies_1k_to_3_5M = linspace(1000, 3.5e6, 300);

% Frequências de 3,5M a 10M
frequencies_3_5M_to_10M = linspace(3.5e6, 10e6, 200);

% Frequências acima de 10M
frequencies_above_10M = linspace(10e6, 2.2e9, 300);

% Geração dos níveis de fase noise correspondentes
phase_noise_1k_to_1M = -89 * ones(size(frequencies_1k_to_3_5M));
phase_noise_1M_to_10M = -124 * ones(size(frequencies_3_5M_to_10M));
phase_noise_above_10M = -132 * ones(size(frequencies_above_10M));

% Concatenação dos dados
freq = [frequencies_1k_to_3_5M, frequencies_3_5M_to_10M, frequencies_above_10M];
phase_noise = [phase_noise_1k_to_1M, phase_noise_1M_to_10M, phase_noise_above_10M];


x = readmatrix(nome_arquivo_csv);

Fs = 26e6   
Fs = 2e9 *10
Ts = 1/Fs
kv = Fs
% x = x/Fs;
% phase = filter (Ts *2* pi*kv ,[1 -1], x - mean ( x ));
phase = filter (1 ,[1 -1], x - mean ( x ));
% phase = x;
media = mean(x)
max(phase)
min(phase)
num_segments = 16
window_length = floor ( length ( phase )/ num_segments );
[ PSDphase ,f] = pwelch(phase , window_length ,[] ,[] ,Fs);%,'twosided');
PSDphase = PSDphase .* (sinc(f/Fs)).^2;	% correct for ZOH
PSDphase = 10*log10(PSDphase);


% Nome do arquivo CSV
filename = 'C:\Users\ander\OneDrive - Associacao Antonio Vieira\UNISINOS\TCC\Python/phase_filtered.csv';

% Salvar o array no arquivo CSV
writematrix(phase, filename)
[ PSDphase ,f] = fun_calc_psd(phase, 2e9, 2e3, 700);
% [ PSDphase ,f] = fun_calc_psd(x, 2e9, 100e3, 1e3);

psdmax =  max(PSDphase)
% PSDphase = PSDphase - psdmax;
figure
%h= semilogx ( f, 10* log10 ( PSDphase ), 'b')
h= semilogx ( f, ( PSDphase ), 'b', freq, ( phase_noise ), 'r')
% semilogx ( )
grid on
set(gca , 'fontsize', 15, 'fontweight', 'bold')
xlabel ('Frequency [Hz]')
ylabel ('Phase Noise [dBc/Hz]')
grid on
set(h(1) , 'LineWidth', 2)
%axis ([1e4 0.5e7 -150 -104])
legend ('Matlab Model ', "Phase Noise Mask")


% x = blackman(5265);
% figure
% %h= semilogx ( f, 10* log10 ( PSDphase ), 'b')
% plot ( x,'b')
% grid on

