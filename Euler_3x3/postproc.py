import numpy as np
import matplotlib.pyplot as plt

var = np.loadtxt('shock_pos_sens_uc/s_rho.dat', delimiter='\t')
var2 = np.loadtxt('shock_pos_sens/s_rho.dat', delimiter='\t')
var3 = np.loadtxt('shock_pos_sens_diff/s_rho.dat', delimiter='\t')
x = np.linspace(0,1,1000)
print(var-var2)
plt.plot(x, var, '*')
plt.plot(x, var2)
plt.plot(x, var3)
plt.show()
