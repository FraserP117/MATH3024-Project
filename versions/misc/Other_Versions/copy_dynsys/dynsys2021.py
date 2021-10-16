import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

params = { 'figure.figsize': (8,8),
          'axes.labelsize': 40,
          'lines.linewidth': 2,
          'font.size': 12,
          'lines.color': 'r',
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'legend.fontsize': 10,
          'text.usetex': False,
          'font.family':['sans-serif'],
          'font.sans-serif': ['DejaVu Sans',
                    'Bitstream Vera Sans',
                    'Computer Modern Sans Serif',
                    'Lucida Grande',
                    'Verdana',
                    'Geneva',
                    'Lucid',
                    'Arial',
                    'Helvetica',
                    'Avant Garde',
                    'sans-serif'],
          'mathtext.bf': 'helvetica:bold',
          'xtick.major.pad': 6,
          'ytick.major.pad': 6,
          'xtick.major.size': 5,
          'ytick.major.size': 5,
          'xtick.minor.size': 3,      # minor tick size in points
          'xtick.major.width': 1.,    # major tick width in points
          'xtick.minor.width': 1.,    # minor tick width in points
	      'xtick.labelsize': 'small',
          'ytick.minor.size': 3,      # minor tick size in points
          'ytick.major.width': 1.,    # major tick width in points
          'ytick.minor.width': 1.,    # minor tick width in points
	      'ytick.labelsize':'small'
           }

plt.rcParams.update(params)


# Implementation of the system dynamics
def dynamics(x,t,alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68):
    xdot = np.zeros(3)
    # xdot[0] = alpha_c*(x[1]-x[0]-fx)
    # xdot[1] = x[0]-x[1]+x[2]
    # xdot[2] = -beta_c*x[1]
    # print(f"type(x): {type(x)}")
    # print(f"x: {x}")
    # x = str(x)
    # x_list = [int(e) for e in str(x) if e.isdigit()]
    fx = b_c*x[0] + 0.5*(a_c-b_c)*(np.abs(x[0]+1.0)-np.abs(x[0]-1.0))


    # linear = b_c*x_list[0]
    # diff = 0.5*(a_c-b_c)
    # abs_1 = np.abs(x_list[0]+1.0)
    # abs_2 = np.abs(x_list[0]-1.0)
    # abs_mult = abs_1 * abs_2
    # end_prod = diff*abs_mult
    # fx = linear + end_prod

    # fx = b_c*x[0] + 0.5*(a_c-b_c)*(np.abs(x[0]+1.0)-np.abs(x[0]-1.0))
    f_x_1 = x_list[1]
    f_x_2 = x_list[0]
    fx_fx = float(fx)
    xdot[0] = alpha_c*(f_x_1 - f_x_2 - fx_fx)
    # xdot[0] = alpha_c*(float(x[1])-float(x[0])-float(fx))
    # print(f"\n\nx_list: {x_list}")
    # print(f"type(x_list): {type(x_list)}\n\n")
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]
    return xdot

    # linear = b_c*float(x[0])
    # diff = 0.5*(a_c-b_c)
    # abs_1 = np.abs(float(x[0])+1.0)
    # abs_2 = np.abs(float(x[0])-1.0)
    # abs_mult = abs_1 * abs_2
    # end_prod = diff*abs_mult
    # fx = linear + end_prod
    # # fx = b_c*x[0] + 0.5*(a_c-b_c)*(np.abs(x[0]+1.0)-np.abs(x[0]-1.0))
    # print("ONE")
    # f_x_1 = float(x[1])
    # print("TWO")
    # f_x_2 = float(x[0])
    # fx_fx = float(fx)
    # xdot[0] = alpha_c*(f_x_1 - f_x_2 - fx_fx)
    # # xdot[0] = alpha_c*(float(x[1])-float(x[0])-float(fx))
    # xdot[1] = float(x[0])-float(x[1])+float(x[2])
    # xdot[2] = -beta_c*float(x[1])
    # return xdot

x_init = np.random.uniform(0.1,0.5,3)
# print(f"\n\nx_init: {x_init}")
# print(f"\n\ntype(x_init): {type(x_init)}\n\n")
t_init = 0; t_final = 100; t_step = 0.01
tpoints = np.arange(t_init, t_final, t_step)

transient=int(0.8*len(tpoints))

alpha_c=10.0; beta_c=14.87; a_c=-1.27; b_c=-0.68; sigma=0.03
# y = odeint(dynamics, x_init, tpoints,args=(alpha_c,beta_c,a_c,b_c), hmax = 0.01)
y = odeint(dynamics, x_init, tpoints, args=(alpha_c,beta_c,a_c,b_c), full_output = 1, hmax = 0.01) #, tfirst = True)
'''
odeint(
    func, y0, t, args=(), Dfun=None, col_deriv=0,
    full_output=0, ml=None, mu=None, rtol=None,
    atol=None, tcrit=None, h0=0.0, hmax=0.0,
    hmin=0.0, ixpr=0, mxstep=0, mxhnil=0, mxordn=12,
    mxords=5, printmessg=0, tfirst=False
)
'''
plt.figure()

# print(f"\n\ny[0]: {y[0]}")
# print(f"y[1]: {y[1]}")
# print(f"len(y): {len(y)}")
# print(f"type(y): {type(y)}\n\n")

# plt.plot(y[:,0],y[:,1],'k') # TypeError: tuple indices must be integers or slices, not tuple
# plt.plot(y[:0],y[:1],'k') # TypeError: tuple indices must be integers or slices, not tuple

# print(f"np.array(y[:0]): {np.array(y[:0])}")
# print(f"len(np.array(y[:0])): {len(np.array(y[:0]))}")

# -------------------

# print(f"\n\nnp.array(y[:1][0]): {np.array(y[:1][0])}")
# print(f"len(np.array(y[:1][0])): {len(np.array(y[:1][0]))}\n\n")
#
# print(f"\n\ntpoints: {tpoints}")
# print(f"len(tpoints): {len(tpoints)}\n\n")

plt.plot(tpoints, np.array(y[:1][0]),'k')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.xlim([-4,4])
plt.ylim([-2,2])
plt.tight_layout()
plt.show()
