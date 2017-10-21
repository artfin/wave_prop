import matplotlib.pyplot as plt
import numpy as np

NPoints = 100 
T = 2250 # sampling time

x = np.linspace( 0, NPoints * T, NPoints )
#print(x)

signal = np.sin( x )
#print(signal)

signal_f = abs(np.fft.fft(signal))**2
print(signal_f)

freqs = 2 * np.pi * np.linspace( 0.0, 1.0 / (2.0 * T), NPoints/2 )
#print(freqs)

plt.plot( freqs, signal_f[0:NPoints/2] )
plt.show()
