import numpy as np
import matplotlib.pyplot as plt
import pylab

fileIn = 'fort.100'
data = np.genfromtxt(fileIn);

print np.mean(data)
print np.std(data)**2

plt.figure()
plt.hist( data)
plt.show()
