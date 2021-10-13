import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.axes import Axes

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
def dynamics(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68):

    xdot = np.zeros(3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))

    # Constructing the vector of dynamicsl variables:
    xdot[0] = alpha_c * (x[1] - x[0] - fx)
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    return xdot

# # x-couppled system dynamics: 1.b
# def dynamics_x_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):
#
#     xdot = np.zeros(3)
#     xdot_prime = np.zeros(3)
#     edot = np.zeros(3)
#
#     # initalize the x prime vector
#     x_prime = np.random.uniform(0.1,0.5,3) # this shoud ot be here !!!!!!!!!!!!!!
#
#     fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
#     fx_prime = b_c * x_prime[0] + 0.5 * (a_c-b_c) * (np.abs(x_prime[0] + 1.0) - np.abs(x_prime[0] - 1.0))
#
#     # The first couppled system:
#     xdot[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x_prime[0] - x[0]) # x-coupling
#     xdot[1] = x[0] - x[1] + x[2]
#     xdot[2] = -beta_c * x[1]
#
#     # The second couppled system:
#     xdot_prime[0] = alpha_c * (x_prime[1] - x_prime[0] - fx_prime) + sigma * (x[0]- x_prime[0]) # x-coupling
#     xdot_prime[1] = x_prime[0] - x_prime[1] + x_prime[2]
#     xdot_prime[2] = -beta_c * x_prime[1]
#
#     '''
#     return the 6 D system
#     integrate the 6D system later
#
#     get rid of x_prime for the init, simply make the x argument 6 dimensional
#
#     the error dynamics are then calculated as the diferences between the correct indices of y - see below
#     '''
#
#     # The errors between the respective components in the couppled systems:
#     ex = x[0] - x_prime[0] # don't need this
#     ey = x[1] - x_prime[1] # don't need this
#     ez = x[2] - x_prime[2] # don't need this
#
#     # The error dynamics:
#     # ex = x_dot[0] - x_dot_prime[0]
#     # ey = x_dot[1] - x_dot_prime[1]
#     # ez = x_dot[2] - x_dot_prime[2]
#     edot[0] = xdot[0] - xdot_prime[0]
#     edot[1] = xdot[1] - xdot_prime[1]
#     edot[2] = xdot[2] - xdot_prime[2]
#
#     # The subsequent error dynamics:
#     # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
#     # edot[1] = ex - ey + ez
#     # edot[2] = -beta_c * ey
#
#     return edot



# x-couppled system dynamics: 1.b
def dynamics_x_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    x_dot_couppled_system = np.zeros(6)

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x_dot_couppled_system[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[3] - x[0]) # x-coupling
    x_dot_couppled_system[1] = x[0] - x[1] + x[2]
    x_dot_couppled_system[2] = -beta_c * x[1]
    x_dot_couppled_system[3] = alpha_c * (x[4] - x[3] - fx) + sigma * (x[0]- x[3]) # x-coupling
    x_dot_couppled_system[4] = x[3] - x[4] + x[5]
    x_dot_couppled_system[5] = -beta_c * x[3]

    return x_dot_couppled_system


    '''
    return the 6 D system
    integrate the 6D system later

    get rid of x_prime for the init, simply make the x argument 6 dimensional

    the error dynamics are then calculated as the diferences between the correct indices of y - see below
    '''

    # # The errors between the respective components in the couppled systems:
    # ex = x[0] - x_prime[0] # don't need this
    # ey = x[1] - x_prime[1] # don't need this
    # ez = x[2] - x_prime[2] # don't need this
    #
    # # The error dynamics:
    # # ex = x_dot[0] - x_dot_prime[0]
    # # ey = x_dot[1] - x_dot_prime[1]
    # # ez = x_dot[2] - x_dot_prime[2]
    # edot[0] = xdot[0] - xdot_prime[0]
    # edot[1] = xdot[1] - xdot_prime[1]
    # edot[2] = xdot[2] - xdot_prime[2]

    # The subsequent error dynamics:
    # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
    # edot[1] = ex - ey + ez
    # edot[2] = -beta_c * ey

    return edot


# y-couppled system dynamics: 1.c
def dynamics_y_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    edot = np.zeros(3)

    # initalize the x prime vector
    x_prime = np.random.uniform(0.1,0.5,3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x_prime[0] + 0.5 * (a_c-b_c) * (np.abs(x_prime[0] + 1.0) - np.abs(x_prime[0] - 1.0))

    # The first couppled system:
    xdot[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[1] - x_prime[1]) # y-coupling
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    # The second couppled system:
    xdot[0] = alpha_c * (x_prime[1] - x_prime[0] - fx_prime) + sigma * (x_prime[1] - x[1]) # y-coupling
    xdot[1] = x_prime[0] - x_prime[1] + x_prime[2]
    xdot[2] = -beta_c * x_prime[1]

    # The errors between the respective components in the couppled systems:
    ex = x[0] - x_prime[0] # don't need this
    ey = x[1] - x_prime[1] # don't need this
    ez = x[2] - x_prime[2] # don't need this

    # The subsequent error dynamics:
    # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
    # edot[1] = ex - ey + ez
    # edot[2] = -beta_c * ey
    edot[0] = xdot[0] - xdot_prime[0]
    edot[1] = xdot[1] - xdot_prime[1]
    edot[2] = xdot[2] - xdot_prime[2]

    return edot

# z-couppled system dynamics: 1.c
def dynamics_z_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    edot = np.zeros(3)

    # initalize the x prime vector
    x_prime = np.random.uniform(0.1,0.5,3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x_prime[0] + 0.5 * (a_c-b_c) * (np.abs(x_prime[0] + 1.0) - np.abs(x_prime[0] - 1.0))

    # The first couppled system:
    xdot[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[2] - x_prime[2]) # z-coupling
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    # The second couppled system:
    xdot[0] = alpha_c * (x_prime[1] - x_prime[0] - fx_prime) + sigma * (x_prime[2] - x[2]) # z-coupling
    xdot[1] = x_prime[0] - x_prime[1] + x_prime[2]
    xdot[2] = -beta_c * x_prime[1]

    # The errors between the respective components in the couppled systems:
    ex = x[0] - x_prime[0] # don't need this
    ey = x[1] - x_prime[1] # don't need this
    ez = x[2] - x_prime[2] # don't need this

    # The subsequent error dynamics:
    # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
    # edot[1] = ex - ey + ez
    # edot[2] = -beta_c * ey
    edot[0] = xdot[0] - xdot_prime[0]
    edot[1] = xdot[1] - xdot_prime[1]
    edot[2] = xdot[2] - xdot_prime[2]

    return edot

# x-couppled system dynamics: 1.b
def dynamics_x_coupling_stochastic(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    edot = np.zeros(3)

    # Gaussian white noise:
    # is this the same as xi? doe we simply call np.random.normal again?


    # samples from random normal:
    xi = np.random.normal(0, 1)

    # initalize the x prime vector
    x_prime = np.random.uniform(0.1,0.5,3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x_prime[0] + 0.5 * (a_c-b_c) * (np.abs(x_prime[0] + 1.0) - np.abs(x_prime[0] - 1.0))

    # The first couppled system:
    xdot[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[0] - x_prime[0]) # x-coupling
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    # The second couppled system:
    xdot[0] = alpha_c * (x_prime[1] - x_prime[0] - fx_prime) + sigma * (x_prime[0] - x[0]) # x-coupling
    xdot[1] = x_prime[0] - x_prime[1] + x_prime[2]
    xdot[2] = -beta_c * x_prime[1]

    # The errors between the respective components in the couppled systems:
    ex = x[0] - x_prime[0]
    ey = x[1] - x_prime[1]
    ez = x[2] - x_prime[2]

    # The subsequent error dynamics:
    # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
    # edot[1] = ex - ey + ez
    # edot[2] = -beta_c * ey
    edot[0] = xdot[0] - xdot_prime[0]
    edot[1] = xdot[1] - xdot_prime[1]
    edot[2] = xdot[2] - xdot_prime[2]

    return edot

# question 1.e.i
def dynamics_x_driver(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    xdot = np.zeros(3)
    xdot_prime = np.zeros(3)
    edot = np.zeros(3)

    # initalize the x prime vector
    x_prime = np.random.uniform(0.1,0.5,3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x_prime[0] + 0.5 * (a_c-b_c) * (np.abs(x_prime[0] + 1.0) - np.abs(x_prime[0] - 1.0))

    # The first couppled system:
    xdot[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[0] - x_prime[0]) # x-coupling
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    # The second couppled system:
    xdot[0] = alpha_c * (x_prime[1] - x_prime[0] - fx_prime) + sigma * (x_prime[0] - x[0]) # x-coupling
    xdot[1] = x_prime[0] - x_prime[1] + x_prime[2]
    xdot[2] = -beta_c * x_prime[1]

    # The errors between the respective components in the couppled systems:
    ex = x[0] - x_prime[0]
    ey = x[1] - x_prime[1]
    ez = x[2] - x_prime[2]

    # The subsequent error dynamics:
    # edot[0] = (-alpha_c - alpha_c * (fx - fx_prime) - 2 * sigma) * ex + alpha_c * ey
    # edot[1] = ex - ey + ez
    # edot[2] = -beta_c * ey

    edot[0] = xdot[0] - xdot_prime[0]
    edot[1] = xdot[1] - xdot_prime[1]
    edot[2] = xdot[2] - xdot_prime[2]

    return edot

# def init(dynamics, sigma):
#     # randomly initialize the state vector
#     x_init = np.random.uniform(0.1,0.5,3) # add 6 here
#
#     # initalise the time, final time and step size
#     t_init = 0; t_final = 100; t_step = 0.01
#
#     # construct the range of time values
#     tpoints = np.arange(t_init, t_final, t_step)
#
#     # ????????????????????????????
#     transient = int(0.8 * len(tpoints))
#
#     # assign values to all parameters for the dynamics
#     alpha_c = 10.0; beta_c = 14.87; a_c = -1.27; b_c = -0.68; sigma = sigma
#
#     # yield trajectories by integrating the dynamics in time for the number of points specified above
#     # y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c), full_output = 1, hmax = 0.01)
#     y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c, sigma), full_output = 1, hmax = 0.01)
#
#     return y, tpoints

def error_dynamics(y, tpoints):
    e_x = y[0] - y[3]
    e_y = y[1] - y[4]
    e_z = y[2] - y[5]

# def init(dynamics, sigma):
#     # randomly initialize the state vector
#     x_init = np.random.uniform(0.1,0.5,6) # make this 6D
#
#     # initalise the time, final time and step size
#     t_init = 0; t_final = 100; t_step = 0.01
#
#     # construct the range of time values
#     tpoints = np.arange(t_init, t_final, t_step)
#
#     # ????????????????????????????
#     transient = int(0.8 * len(tpoints))
#
#     # assign values to all parameters for the dynamics
#     alpha_c = 10.0; beta_c = 14.87; a_c = -1.27; b_c = -0.68; # sigma = sigma
#
#     # yield trajectories by integrating the dynamics in time for the number of points specified above
#     # y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c), full_output = 1, hmax = 0.01)
#     y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c), full_output = 1, hmax = 0.01)
#
#     # calculate the error dynamics via y
#     '''
#     y = [x, y, z, x_prime, y_prime, z_prime]
#     e_x = y[0] - y[3]
#     e_y = y[1] - y[4]
#     e_z = y[2] - y[5]
#     '''
#
#     return y, tpoints

def init(dynamics, sigma):
    # randomly initialize the state vector
    x_init = np.random.uniform(0.1,0.5,6) # make this 6D

    # initalise the time, final time and step size
    t_init = 0; t_final = 100; t_step = 0.01

    # construct the range of time values
    tpoints = np.arange(t_init, t_final, t_step)

    # ????????????????????????????
    transient = int(0.8 * len(tpoints))

    # assign values to all parameters for the dynamics
    alpha_c = 10.0; beta_c = 14.87; a_c = -1.27; b_c = -0.68; # sigma = sigma

    # yield trajectories by integrating the dynamics in time for the number of points specified above
    # y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c), full_output = 1, hmax = 0.01)
    y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c, sigma), full_output = 1, hmax = 0.01)

    print(y[:1][0])
    print(y[:1][0].shape)

    '''
    Odeint Returns
    -------
    y : array, shape (len(t), len(y0))
        Array containing the value of y for each desired time in t,
        with the initial value `y0` in the first row.
    '''

    # # calculate the error dynamics via y
    # # y = [x, y, z, x_prime, y_prime, z_prime]
    # edot[0] = y[0] - y[3]
    # edot[0] = y[1] - y[4]
    # edot[0] = y[2] - y[5]

    return y, tpoints
    # return edot, tpoints




def display_dynamics(y, tpoints, title, color):
    # create the figure; onto which we will plot the above trajectries
    plt.figure(facecolor = 'grey')

    # plt.plot(tpoints, np.array(y[:1][0]),'k', color = "green")
    plt.plot(tpoints, np.array(y[:1][0]), color = color, linewidth = 2)
    plt.title(title)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.xlim([0,20])
    plt.ylim([-5000,5500])
    plt.tight_layout()
    plt.show()

# def display_component_error_dynamics(y, tpoints, title, color):
#     # create the figure; onto which we will plot the above trajectries
#     plt.figure(facecolor = 'grey')
#
#     # plt.plot(tpoints, np.array(y[:1][0]),'k', color = "green")
#     plt.plot(tpoints, np.array(y[:1][0]), color = color)
#     plt.title(title)
#     plt.xlabel(r"$x$")
#     plt.ylabel(r"$y$")
#     plt.xlim([0,1.5])
#     plt.ylim([-6,6])
#     plt.tight_layout()
#     plt.show()

# def display_system_dynamics(X, tpoints, title):
#     # create the figure; onto which we will plot the above trajectries
#
#     plt.subplot(1, 3, 1)
#     plt.plot(tpoints, np.array(X[0][:1][0]), color = "blue")
#
#     plt.subplot(1, 3, 2)
#     plt.plot(tpoints, np.array(X[1][:1][0]), color = "green")
#
#     plt.subplot(1, 3, 3)
#     plt.plot(tpoints, np.array(X[2][:1][0]), color = "red")
#
#     plt.tight_layout()
#
#     plt.show()

    # # plt.plot(tpoints, np.array(y[:1][0]),'k', color = "green")
    # plt.plot(tpoints, np.array(y[:1][0]), color = "red")
    # plt.title(title)
    # plt.xlabel(r"$x$")
    # plt.ylabel(r"$y$")
    # plt.xlim([0,100])
    # plt.ylim([-6,6])
    # plt.tight_layout()
    # plt.show()


if __name__ == '__main__':
    # # The origional system dynamics
    # y, tpoints = init(dynamics, sigma = 0.5)
    # display_dynamics(y, tpoints, "System Dynamics", "green")

    # x-coupled dynamics:
    y, tpoints = init(dynamics_x_coupling, sigma = 0.5)
    for x in y[:1][0]:
        print(x)
    display_dynamics(y, tpoints, "System Dynamics: x-coupling", "blue")

    # # y-coupled dynamics:
    # y, tpoints = init(dynamics_y_coupling, sigma = 0.5)
    # # display_dynamics(y, tpoints, "System Dynamics: y-coupling", "red")
    # display_component_error_dynamics(y, tpoints, "System Dynamics: y-coupling", "red")

    # # z-coupled dynamics:
    # y, tpoints = init(dynamics_z_coupling)
    # display_dynamics(y, tpoints, "System Dynamics: z-coupling", "red")

    # # x-coupled dynamics:
    # y, tpoints = init(dynamics_x_driver, 5.0)
    # display_dynamics(y, tpoints, "System Dynamics: x-driving", "blue")

    # # whole system dynamics:
    # x, tpoints = init(dynamics_x_coupling)
    # y, tpoints = init(dynamics_y_coupling)
    # z, tpoints = init(dynamics_z_coupling)
    # X = [x, y, z]
    # display_system_dynamics(X, tpoints, "System Dynamics")

# if __name__ == '__main__':
#     # # The origional system dynamics
#     # y, tpoints = init(dynamics, sigma = 0.5)
#     # display_dynamics(y, tpoints, "System Dynamics", "green")
#
#     # x-coupled dynamics:
#     y, tpoints = init(dynamics_x_coupling, sigma = 0.5)
#     # display_dynamics(y, tpoints, "System Dynamics: x-coupling", "blue")
#     display_dynamics(y, tpoints, "System Dynamics: x-coupling", "blue")
#
#     # # y-coupled dynamics:
#     # y, tpoints = init(dynamics_y_coupling, sigma = 0.5)
#     # # display_dynamics(y, tpoints, "System Dynamics: y-coupling", "red")
#     # display_component_error_dynamics(y, tpoints, "System Dynamics: y-coupling", "red")
#
#     # # z-coupled dynamics:
#     # y, tpoints = init(dynamics_z_coupling)
#     # display_dynamics(y, tpoints, "System Dynamics: z-coupling", "red")
#
#     # # x-coupled dynamics:
#     # y, tpoints = init(dynamics_x_driver, 5.0)
#     # display_dynamics(y, tpoints, "System Dynamics: x-driving", "blue")
#
#     # # whole system dynamics:
#     # x, tpoints = init(dynamics_x_coupling)
#     # y, tpoints = init(dynamics_y_coupling)
#     # z, tpoints = init(dynamics_z_coupling)
#     # X = [x, y, z]
#     # display_system_dynamics(X, tpoints, "System Dynamics")
