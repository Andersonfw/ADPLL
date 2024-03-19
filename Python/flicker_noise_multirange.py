import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import control
import sympy as sp
import warnings
from scipy import signal

def fun_calc_psd(x , fs=1 , rbw=100e3 , fstep=None):
    '''
    Calculate power spectral density

    INPUT arguments
    x     :  data vector
    fs    :  sample rate [Hz]
    rbw   :  resolution bandwidth [Hz]
    fstep :  FFT frequency bin separation [Hz]
    OUTPUT
    XdB	: spectrum of x [dB]
    f	: frequency vector [Hz]
    '''

    if fstep is None:
        fstep = rbw / 1.62
    len_x = len(x)
    nwin = round(fs * 1.62 / rbw)
    nfft = round(fs / fstep)
    if nwin > len_x:
        nwin = len_x
        rbw = fs * 1.62 / nwin
    num_segments = 8
    # nwin = math.floor(len(x) / num_segments)
    fftstr = (f'len(x)={len_x:.2f}, rbw={rbw / 1e3:.2f}kHz, fstep={fstep / 1e3:.2f}kHz, nfft={nfft:d}, nwin={nwin:d}')
    print(f'Calculating the PSD: {fftstr} ...')
    f , X = signal.welch(x , fs=fs , window=signal.windows.blackman(nwin) , nperseg=nwin , nfft=nfft ,scaling='density')
    X *= (np.sinc(f / fs)) ** 2  # correct for ZOH
    XdB = 10 * np.log10(X)
    XdB_sig = np.max(XdB)
    print(f'Signal PSD peak = {XdB_sig:.2f} dB, 10log(rbw) = {10 * np.log10(rbw):.1f}')
    return XdB , f

def extract_parameters(H):
    # Converter a string em uma expressão simbólica
    expressao = sp.sympify(H)
    # Obter os numeradores e denominadores como objetos simbólicos
    num_Hol, den_Hol = expressao.as_numer_denom()
    # print("numerador em s:", num_Hol)
    # print("denominador em s:", den_Hol)
    # Extrair coeficientes do numerador e denominador de Hs
    num_Hol = sp.Poly(num_Hol, s).all_coeffs()
    den_Hol = sp.Poly(den_Hol, s).all_coeffs()
    # print("numerador:", num_Hol)
    # print("denominador:", den_Hol)
    # Transforma em array para aplicar a transformada de laplace
    num = np.array(num_Hol, dtype=float)
    den = np.array(den_Hol, dtype=float)
    return num, den


# Definição de 's' como um símbolo
s = sp.symbols('s')

fr = 33e6  # Frequeência de referência
# a = 2 ** -7  # alpha value
a = 2 ** -6  # alpha value
# p = 2 ** -14  # rho value
p = 2 ** -15  # rho value
# Definição da função de transferência em S
z = 1 + s/fr
A = 10**(10/20) 
H_s = (a + p * fr / s) * (fr / s)
ak = []
fc_k = [100, 1e3, 10e3, 100e3, 1e6]
h_flicker = []
tf_f = []
bode = []
tf_all = 1
for i in range(len(fc_k)):
    ak.append(2 *np.pi * fc_k[i]/fr)
    h_flicker.append((ak[i]* A**(-i+1)*z)/(z-(1-ak[i])))
    num_H, den_H = extract_parameters(h_flicker[i])
    tf_f.append(control.tf(num_H, den_H))
    bode.append(control.bode(tf_f[i], Hz=True, dB=True, deg=True, plot=True))
    tf_all += tf_f[i]

# Transformada de Fourier da função de transferência em S
omega = sp.symbols('omega', real=True)
H_jw = H_s.subs(s, sp.I * omega)

# Valores de frequência angular para avaliação
frequencias = np.logspace(1, 3, 1000)  # Escolha suas frequências de interesse

# Avaliação da função de transferência em frequências angulares
#H_jw_valores = [H_jw.evalf(subs={omega: w}) for w in frequencias]

# Plot da magnitude da função de transferência no domínio da frequência
#magnitude = [abs(valor) for valor in H_jw_valores]
y = []

y.append(0)
noise = np.random.randn(len(frequencias))
for i in range(1 , len(frequencias)-1):
    ak = frequencias[i+1]/fr
    yk = (1 - ak)*y[i-1] + ak*A**(-i+1)*noise[i]
    y.append(yk) 
y.append(y[i-1])
plt.plot(frequencias, y)
plt.xlabel('Frequência Angular (rad/s)')
plt.ylabel('Magnitude')
plt.title('Magnitude da função de transferência no domínio da frequência')
plt.grid(True)

plt.figure()
bode.append(control.bode(tf_all, Hz=True, dB=True, deg=True, plot=True))
# Xdb_o , f = fun_calc_psd(y , fr ,  1, 1e3)
# plt.figure()
# plt.semilogx(f , Xdb_o , label="Phase Noise")
# plt.grid(visible=True)
# plt.legend()
# plt.xlabel('Freq (Hz)')
# plt.ylabel('Phase Noise [dBc/Hz]')
plt.show()
