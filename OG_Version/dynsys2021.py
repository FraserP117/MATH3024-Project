import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

params = { 'figure.figsize': (8,8),
          'axes.labelsize': 40,
          'lines.linewidth': 2,
          #'text.fontsize': 12,
          #'lines.color': 'r',
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          # 'legend.fontsize': 10,
          # 'title.fontsize': 12,
          # 'text.usetex': False,
          # 'font': 'Helvetica',
          #'mathtext.bf': 'helvetica:bold',
          # 'xtick.major.pad': 6,
          # 'ytick.major.pad': 6,
          # 'xtick.major.size': 5,
          # 'ytick.major.size': 5,
          # 'xtick.minor.size': 3,      # minor tick size in points
          # 'xtick.major.width': 1.,    # major tick width in points
          # 'xtick.minor.width': 1.,    # minor tick width in points
          # 'ytick.minor.size': 3,      # minor tick size in points
          # 'ytick.major.width': 1.,    # major tick width in points
          # 'ytick.minor.width': 1.,    # minor tick width in points
          # 'tick.labelsize': 'small'
           }

plt.rcParams.update(params)



def dynamics(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68):
	xdot = np.zeros(3)
	fx = b_c*x[0] + 0.5*(a_c-b_c)*(np.abs(x[0]+1.0)-np.abs(x[0]-1.0))
	xdot[0] = alpha_c*(x[1]-x[0]-fx)
	xdot[1] = x[0]-x[1]+x[2]
	xdot[2] = -beta_c*x[1]
	return xdot

x_init = np.random.uniform(0.1,0.5,3)

t_init=0; t_final = 100; t_step = 0.01
tpoints = np.arange(t_init,t_final,t_step)

transient=int(0.8*len(tpoints))

alpha_c=10.0; beta_c=14.87; a_c=-1.27; b_c=-0.68; sigma=0.03
y = odeint(dynamics, x_init, tpoints,args=(alpha_c,beta_c,a_c,b_c), hmax = 0.01)

plt.figure()
plt.plot(y[:,0],y[:,1],'k')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.xlim([-4,4])
plt.ylim([-2,2])
plt.tight_layout()
