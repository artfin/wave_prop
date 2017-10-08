import numpy as np

x = np.linspace( -1, 1, 4)
print('x: {0}'.format(x))

p = np.fft.fftfreq( 4, d = 1.5 )
print('p: {0}'.format(p))

p = np.fft.fftshift(p)
print('p: {0}'.format(p))
