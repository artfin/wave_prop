import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

# Number of samplepoints
N = 600
# sample spacing
T = 1.0 / 800.0

x = np.linspace(0.0, N * T, N)
y = 10 + np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0 * T), N/2)

#plt.subplot(2, 1, 1)
#plt.plot(xf, 2.0/N * np.abs(yf[0:N/2]))
#plt.subplot(2, 1, 2)
plt.plot(xf[1:], 2.0/N * np.abs(yf[0:N/2])[1:])
plt.show()
