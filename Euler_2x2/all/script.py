import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

test_case = "results/isolated_shock/sens_CD/"
u = np.loadtxt(test_case+"u.dat")
s_u = np.loadtxt(test_case+"s_u.dat")
tau = np.loadtxt(test_case+"tau.dat")
s_tau = np.loadtxt(test_case+"s_tau.dat")

test_case_ref = "results/isolated_shock/double_CD/"
s_u_ref = np.loadtxt(test_case_ref+"s_u.dat")
s_tau_ref = np.loadtxt(test_case_ref+"s_tau.dat")

test_case2 = "results/isolated_shock/state_CD/"
u_2 = np.loadtxt(test_case2+"u.dat")
s_u_2 = np.loadtxt(test_case2+"s_u.dat")
tau_2 = np.loadtxt(test_case2+"tau.dat")
s_tau_2 = np.loadtxt(test_case2+"s_tau.dat")

test_case3 = "results/isolated_shock/"
u_3 = np.loadtxt(test_case3+"u.dat")
s_u_3 = np.loadtxt(test_case3+"s_u.dat")
tau_3 = np.loadtxt(test_case3+"tau.dat")
s_tau_3 = np.loadtxt(test_case3+"s_tau.dat")


N_x = np.size(u_3[0, :])
N_t = np.size(u_3[:, 0])
dx = 1./N_x

xx = np.linspace(0.5*dx,1-0.5*dx,N_x);
fig, ax = plt.subplots(1,2)
#
#ax[0].clear()
#ax[1].clear()
#ax[0].plot(xx, s_tau[-1,:])
#ax[0].plot(xx, s_tau_ref[-1,:])
#ax[0].set_ylim([np.min(np.min(s_tau))-.1, np.max(np.max(s_tau))+.1])
#ax[1].plot(xx, s_u[-1,:])
#ax[1].plot(xx, s_u_ref[-1,:])
#ax[1].set_ylim([np.min(np.min(s_u))-.1, np.max(np.max(s_u))+.1])
#plt.show()

def animate(i):
    ax[0].clear()
    ax[1].clear()
    #ax[0].plot(xx, s_tau[i,:], label='sens_CD')
    #ax[0].plot(xx, s_tau_2[i,:], label='state_CD')
    ax[0].plot(xx, s_tau_3[i,:], label='diff sampl')
    ax[0].set_xlim([0,1])
    ax[0].legend()
    #ax[0].plot(xx, s_tau_ref[i,:])
    ax[0].set_ylim([np.min(np.min(s_tau))-.1, np.max(np.max(s_tau))+.1])
    #ax[1].plot(xx, s_u[i,:], label='sens_CD')
    #ax[1].plot(xx, s_u_2[i,:], label='state_CD')
    ax[1].plot(xx, s_u_3[i,:], label='diff sampl')
    #ax[1].plot(xx, s_u_ref[i,:])
    ax[1].legend()
    ax[1].set_ylim([np.min(np.min(s_u))-.1, np.max(np.max(s_u))+.1])
    ax[1].set_xlim([0,1])

ani = FuncAnimation(fig, animate, interval=100)
plt.show()
