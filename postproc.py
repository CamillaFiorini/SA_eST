import numpy as np
import matplotlib.pyplot as plt

var = np.loadtxt('results/double_CD/s_rho.dat', delimiter='\t')
var2 = np.loadtxt('results/state_CD/s_rho.dat', delimiter='\t')
var3 = np.loadtxt('results/sens_CD/s_rho.dat', delimiter='\t')
x = np.linspace(0,1,1000)
plt.plot(x, var, label = "double_CD")
plt.plot(x, var2, label = "state_CD")
plt.plot(x, var3, label = "sens_CD")
plt.legend()
plt.show()
