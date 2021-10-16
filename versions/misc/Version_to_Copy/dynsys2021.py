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
          'legend.title_fontsize': 12,
          'text.usetex': False,
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
          'xtick.labelsize': 'small', # this needed to be added
          'ytick.minor.size': 3,      # minor tick size in points
          'ytick.major.width': 1.,    # major tick width in points
          'ytick.minor.width': 1.,    # minor tick width in points
          'ytick.labelsize':'small' # this needed to be added
           }

# update params
plt.rcParams.update(params)

# implement the system dynamics
def dynamics(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68):

    xdot = np.zeros(3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0]+1.0) - np.abs(x[0]-1.0))

    # the vector of dynamical variables
    xdot[0] = alpha_c * (x[1] - x[0] - fx)
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    return xdot

def dynamics_question_b(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()]

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()]

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime) + sigma*(x_list[0] - x_list_prime[0]) # coupple to x
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2]
    xdot_prime[2] = -beta_c*x_list_prime[1]

    # define the following errors:
    ex = xdot[0] - xdot_prime[0]
    ey = xdot[1] - xdot_prime[1]
    ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    ey_dot = ex - ey + ez
    ez_dot = -beta_c*ey
    e_dot = [ex_dot, ey_dot, ez_dot]

    # return the error dynamics
    return e_dot

def dynamics_question_c_yy(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()]

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()]

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime)
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2] + sigma*(x_list[1] - x_list_prime[1]) # coupple to y
    xdot_prime[2] = -beta_c*x_list_prime[1]

    # define the following errors:
    ex = xdot[0] - xdot_prime[0]
    ey = xdot[1] - xdot_prime[1]
    ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    # TO DO
    '''
    ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    ey_dot = ex - ey + ez
    ez_dot = -beta_c*ey
    e_dot = [ex_dot, ey_dot, ez_dot]
    '''

    # return the error dynamics
    return e_dot

def dynamics_question_c_zz(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()]

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()]

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime)
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2]
    xdot_prime[2] = -beta_c*x_list_prime[1] + sigma*(x_list[2] - x_list_prime[2]) # coupple to z

    # define the following errors:
    ex = xdot[0] - xdot_prime[0]
    ey = xdot[1] - xdot_prime[1]
    ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    # TO DO
    '''
    ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    ey_dot = ex - ey + ez
    ez_dot = -beta_c*ey
    e_dot = [ex_dot, ey_dot, ez_dot]
    '''

    # return the error dynamics
    return e_dot

def dynamics_question_d_xx(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68):
    xdot = np.zeros(3)

    x_list = [int(e) for e in str(x) if e.isdigit()]

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0]-x_list[1]+x_list[2]
    xdot[2] = -beta_c*x_list[1]

    return xdot

def dynamics_question_e_x_driver(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    '''
    Pecora-Carrol Couppling for an x-x couppled system.
    '''
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()] # change

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()] # change

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    # fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    # xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime)
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2]
    xdot_prime[2] = -beta_c*x_list_prime[1]

    # define the following errors:
    # ex = xdot[0] - xdot_prime[0]
    ey = xdot[1] - xdot_prime[1]
    ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    # ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    ey_dot = ex - ey + ez
    ez_dot = -beta_c*ey
    # e_dot = [ex_dot, ey_dot, ez_dot]
    e_dot = [xdot[0], ey_dot, ez_dot]

    # return the error dynamics
    return e_dot

def dynamics_question_e_y_driver(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    '''
    Pecora-Carrol Couppling for an x-x couppled system.
    '''
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()] # change

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()] # change

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    # fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    # xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime)
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2]
    xdot_prime[2] = -beta_c*x_list_prime[1]

    # define the following errors:
    ex = xdot[0] - xdot_prime[0]
    # ey = xdot[1] - xdot_prime[1]
    ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    # ey_dot = ex - ey + ez
    ez_dot = -beta_c*ey
    e_dot = [ex_dot, xdot[1], ez_dot]

    # return the error dynamics
    return e_dot

def dynamics_question_e_z_driver(x,t,alpha_c=10.0, beta_c=14.87, a_c=-1.27, b_c=-0.68, sigma=0.03):
    '''
    Pecora-Carrol Couppling for an x-x couppled system.
    '''
    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    x_list = [int(e) for e in str(x) if e.isdigit()] # change

    x_init_prime = np.random.uniform(0.1,0.5,3)
    x_list_prime = [int(e) for e in str(x_init_prime) if e.isdigit()] # change

    fx = b_c*x_list[0] + 0.5*(a_c-b_c)*(np.abs(x_list[0]+1.0)-np.abs(x_list[0]-1.0))
    xdot[0] = alpha_c*(x_list[1]-x_list[0]-fx)
    xdot[1] = x_list[0] - x_list[1] + x_list[2]
    xdot[2] = -beta_c*x_list[1]

    # fx_prime = b_c*x_list_prime[0] + 0.5*(a_c-b_c)*(np.abs(x_list_prime[0]+1.0)-np.abs(x_list_prime[0]-1.0))
    # xdot_prime[0] = alpha_c*(x_list_prime[1]-x_list_prime[0]-fx_prime)
    xdot_prime[1] = x_list_prime[0] - x_list_prime[1] + x_list_prime[2]
    xdot_prime[2] = -beta_c*x_list_prime[1]

    # define the following errors:
    ex = xdot[0] - xdot_prime[0]
    ey = xdot[1] - xdot_prime[1]
    # ez = xdot[2] - xdot_prime[2]

    # error dynamics:
    ex_dot = (-alpha_c - alpha_c*(fx-fx_prime) - 2*sigma)*ex + alpha_c*ey
    ey_dot = ex - ey + ez
    # ez_dot = -beta_c*ey
    e_dot = [ex_dot, ey_dot, xdot[2]]

    # return the error dynamics
    return e_dot

# x_init = np.random.uniform(0.1,0.5,6)
x_init = np.random.uniform(0.1,0.5,3)

# x_init = np.random.uniform(0.1,0.5,3)
# x_init_prime = np.random.uniform(0.1,0.5,3)
# x_inits = [x_init, x_init_prime]

t_init = 0; t_final = 100; t_step = 0.01
tpoints = np.arange(t_init, t_final, t_step) # discrete time intervals at which to numerically intergrate the systems

transient=int(0.8*len(tpoints))

alpha_c=10.0; beta_c=14.87; a_c=-1.27; b_c=-0.68; sigma=0.03
# y = odeint(dynamics, x_init, tpoints,args=(alpha_c,beta_c,a_c,b_c), hmax = 0.01)
# y = odeint(dynamics, x_init, tpoints, args=(alpha_c,beta_c,a_c,b_c,sigma), full_output = 1, hmax = 0.01)
# y = odeint(dynamics, x_inits, tpoints, args=(alpha_c,beta_c,a_c,b_c,sigma), full_output = 1, hmax = 0.01)
y = odeint(dynamics, x_init, tpoints, args=(alpha_c,beta_c,a_c,b_c,sigma), full_output = 1, hmax = 0.01) # the call to the code as is necessary for part b
'''
Returns y
array, shape (len(t), len(y0))
Array containing the value of y for each desired time in t, with the initial value y0 in the first row.
'''

'''
fig, axs = plt.subplots(3)
fig.suptitle('Vertically stacked subplots')
axs[0].plot(tpoints, np.array(y[0][0]),'k')
axs[1].plot(tpoints, np.array(y[0][1]),'k')
axs[2].plot(tpoints, np.array(y[0][2]),'k')
'''

plt.figure()
plt.plot(tpoints, np.array(y[0]),'k') # plot takes the time points
# plt.plot(y[:,0],y[:,1],'k')
# plt.plot(tpoints, y,'k')
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.xlim([-4,4])
plt.ylim([-2,2])
plt.tight_layout()
plt.show()
