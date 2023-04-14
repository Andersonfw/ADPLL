import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import csv
from pprint import pprint
from statistics import mean
import pandas as pd

Shunt = 3
Vin = 5.05
names = [
#	'PA_TX',
#	'PA_SHUTDOWN',
#	'PA_BY_PASS',
	'PA_RX'
]

estados =[ 
#	'Inicialização',
#	'Leitura da Tensão e RTC',
#	'Leitura da Temperatura',
#	'Salva na FLASH',
	'Transmissão',
	'Recepção'
]

medidas = [
	'Etapas',
	'Corrente [mA]',
	'Tempo [ms]',
	'Energia [mJ]',
]

tab = {}

for name in names:
	#Monta a tabela
	for j in medidas:
		tab[j] = []
	#Monta a primeira coluna
	tab['Etapas'] = estados
	
	mX = []
	mY = []
	for i in estados:
		mX.append([])
		mY.append([])
	
	with open(name+'.csv', 'rt') as f:
		reader = csv.reader(f,delimiter=';')	# usa o separador como ';'
		for row in reader: 						# percorre linha por linha
			key = estados.index(row[2])
			mX[key].append(float(row[0]))
			mY[key].append(float(row[1]))


	for i in range(len(estados)):
		x = mX[i]
		y = mY[i]
		print(len(x))
		if len(x):

			#passo a escala para ms e mW
			#x = (Vin-np.array(x) )*np.array(x) *1000/Shunt
			#y = (Vin-np.array(y) )*np.array(y) *1000/Shunt
		
			#Converte para corrente e aplica mA
			x = np.array(x) *1000/Shunt
			y = np.array(y) *1000/Shunt
			z = y*Vin

			CorrenteMedia = np.mean(y)
			PotenciaMedia = np.mean(z)
			dT = max(x) - min(x)
			Energia = PotenciaMedia*dT
			I_trapz = trapz(z,x)

			tab[medidas[1]].append(round(CorrenteMedia,2))
			tab[medidas[2]].append(round(dT,2))
			tab[medidas[3]].append(round(Energia,2))

			plt.title(name)
			plt.ylabel('Corrente [mA]')
			plt.xlabel('Tempo [ms]')
			plt.grid(which='both', axis='both')
			plt.grid(True)
			plt.margins(0.1, 0.1)
			plt.plot(x, y,label= estados[i])
				#+': '+str(round(CorrenteMedia,2))+'mA | Tempo = '+
				#str(round(dT,2))+'ms')
				
			#plt.plot([min(x),max(x)],[media,media],label='Media = '+str(round(media,2))+'mA')
			plt.legend()
	plt.savefig(name+'.pdf')
	plt.show()
	
	df = pd.DataFrame(tab)
	tabLatex = df.to_latex(
		index=False,
		label='tab:'+name,
		caption='Medições em modo '+name.replace('_',' ')+' de transmissão.',
		bold_rows = True,
		)
	with open('Tabelas/tabela'+name+'.tex', 'w' , encoding='utf-8') as f:
		f.write(tabLatex)

exit()

'''
names = [
	'integral.csv',
	'integral-sem-envio.csv',
]

for name in names:
	x = []
	y = []
	with open(name, 'rt') as f:
		reader = csv.reader(f,delimiter=';')	# usa o separador como ';'
		for row in reader: 						# percorre linha por linha
			x.append(float(row[0]))
			y.append(float(row[1]))

	I_trapz = trapz(y,x)
	
	dT = max(x) - min(x)
	media = I_trapz/dT
	
	#passo a escala para mW
	I_trapz *= 1000
	print(round (I_trapz,2),'mW')
	media *= 1000
	print(round (media,2),'mW')

	#passo a escala para ms e mW
	x = np.array(x)*1000
	y = np.array(y)*1000

	plt.title('Potência dissipada')
	plt.ylabel('Potência [mW]')
	plt.xlabel('Tempo [ms]')
	plt.grid(which='both', axis='both')
	plt.grid(True)
	plt.margins(0.1, 0.1)
	plt.plot(x, y,label=
		'Trabalho = '+str(round(I_trapz,2))+'mJ | Tempo = '+
		str(round(max(x),2))+'ms')
		
	plt.plot([min(x),max(x)],[media,media],label='Media = '+str(round(media,2))+'mW')
	plt.legend()
	plt.savefig(name+'.pdf')
	plt.show()

'''


name = 'supercapCarga.csv'
x = []
y = []
with open(name, 'rt') as f:
	reader = csv.reader(f,delimiter=';')	# usa o separador como ';'
	for row in reader: 						# percorre linha por linha
		x.append(float(row[0]))
		y.append(float(row[1]))

ymax = max(y)
ymin = min(y)
dy = ymax - ymin
y1T = dy*0.63
y4T = dy*0.98
for i,v in enumerate(y):
	if v < y1T: x1T = x[i]
	if v < y4T: x4T = x[i]
print('x1T',x1T)
print('x4T',x4T)
print('ymax',ymax)
print('ymin',ymin)

L1Tx=[x[0],x1T,x1T]
L1Ty=[y1T,y1T,0]
L4Tx=[x[0],x4T,x4T]
L4Ty=[y4T,y4T,0]

plt.title('Carga do Capacitor')
plt.ylabel('Tensão [V]')
plt.xlabel('Tempo [s]')
plt.grid(which='both', axis='both')
plt.grid(True)
plt.margins(0.1, 0.1)
plt.plot(x, y,label='Vcap')
plt.plot(L1Tx,L1Ty,label='0.63Vs = '+str(round(y1T,1))+'V | 1τ = '+str(round(x1T,1))+'s')
plt.plot(L4Tx,L4Ty,label='0.98Vs = '+str(round(y4T,1))+'V | 4τ = '+str(round(x4T,1))+'s')
plt.legend()
plt.savefig(name+'.pdf')
plt.show()