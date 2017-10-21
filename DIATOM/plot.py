import matplotlib.pyplot as plt

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    freqs = []
    ints = []

    for line in lines:
        data = line.split()
        freqs.append( float(data[0]) )
        ints.append( float(data[1]) )

    return freqs, ints

freqs, ints = read_file( 'spectrum.txt' )
 
plt.plot( freqs, ints, marker = '^', color = 'k' )
plt.show()
