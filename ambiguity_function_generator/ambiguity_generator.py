#Ambiguity function generator
#author: Devin Huyghebaert
#Date: Jan. 23, 2017

import numpy as np
import cmath as math
import matplotlib.pyplot as plt
#import ps_routines as cp
import random as rand
import scipy

#800 ranges
nrang = 100

#array for storing the code to be analyzed
b = np.zeros((10000),dtype=complex )

#array for storing ambiguity function
b_ambig = np.zeros((nrang*2+10001),dtype=complex)+0.1

#records how many of each phase
phase_1 = 0
phase_2 = 0

#commented lines for generation of pseudo-random code
for num in range(nrang*2+10001):
	random_num = rand.randint(0,1)
	if random_num == 0:
		b_ambig[num] = (1)
		phase_1 = phase_1+1
	else:
		b_ambig[num] = (-1)
		phase_2 = phase_2+1

print phase_1
print phase_2

#output to file if code needs to be stored
#b_ambig[nrang:nrang+10000].tofile("pseudo_random_code_test.txt")

#read in code to be tested
#test_sig = scipy.fromfile(open("pseudo_random_code_test_8.txt"),dtype=scipy.complex64)

test_sig=b_ambig[nrang:nrang+10000]

#assign code to be stored in array
b = test_sig

#take the conjugate of the code
b = np.conj(b)

#generate array to store information for plotting
c = np.zeros( (200,nrang),dtype=float)
c_test_fft = np.zeros((10000,nrang*2),dtype=float)

#constants for the fft
#80000 length code
N=10000
#800000 samples/second
T=1.0/100000.0
x = np.linspace(0.0,N*T,N)

for num in range(nrang*2):
	print num
	#calculate the fft of the code multiplied by the conjugate at a range, add a small number to prevent taking log of 0
	c_test_fft[:,num] = 20*np.log10(abs(np.fft.fft(b_ambig[(num):(num+10000)]*b))+0.000000000001)
	c_test_fft[:,num] = np.fft.fftshift(c_test_fft[:,num])

fft_freq = np.fft.fftfreq(N,T)

plt.imshow(c_test_fft[4961:5040,:].T,vmax=np.amax(c_test_fft),vmin=np.amax(c_test_fft)-35,origin='lower',extent=[np.amin(fft_freq)/200,np.amax(fft_freq)/200,(nrang*-1.5),(nrang*1.5)],aspect='auto')
plt.xlabel('Doppler (Hz)')
plt.title('pseudo-random code test')
plt.colorbar()
plt.ylabel('Range (km)')
plt.savefig('ambig_func_prcode_test.png')
plt.close()
