import numpy as np
import matplotlib.pyplot as plt
fname = 'approx_sin.txt'
final =np.loadtxt(fname)

x=np.linspace(0, 2*3.14, 100)



plt.plot(x,final)
plt.show()
