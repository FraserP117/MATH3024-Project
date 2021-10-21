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

'''
NOTES:

return the 6 D system
integrate the 6D system later

get rid of x_prime for the init, simply make the x argument 6 dimensional

the error dynamics are then calculated as the diferences between the correct indices of y - see below
'''

'''
MAKE SURE TO GRAPH THE ERRORS IN THE X VECTOR
'''


# Implementation of the system dynamics
def dynamics(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68):

    xdot = np.zeros(3)

    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))

    # Constructing the vector of dynamicsl variables:
    xdot[0] = alpha_c * (x[1] - x[0] - fx)
    xdot[1] = x[0] - x[1] + x[2]
    xdot[2] = -beta_c * x[1]

    return xdot

# x-couppled system dynamics: 1.b
def dynamics_x_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 5.5):

    # init the 6D vector for the coupled dynamics
    x_dot_x_couppled_system = np.zeros(6)

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x_dot_x_couppled_system[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[3] - x[0]) # x-coupling
    x_dot_x_couppled_system[1] = x[0] - x[1] + x[2]
    x_dot_x_couppled_system[2] = -beta_c * x[1]
    x_dot_x_couppled_system[3] = alpha_c * (x[4] - x[3] - fx_prime) + sigma * (x[0]- x[3]) # x-coupling
    x_dot_x_couppled_system[4] = x[3] - x[4] + x[5]
    x_dot_x_couppled_system[5] = -beta_c * x[4]

    return x_dot_x_couppled_system

# y-couppled system dynamics: 1.c
def dynamics_y_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    x_dot_y_couppled_system = np.zeros(6)

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x_dot_y_couppled_system[0] = alpha_c * (x[1] - x[0] - fx)
    x_dot_y_couppled_system[1] = x[0] - x[1] + x[2] + sigma * (x[4] - x[1]) # y-coupling
    x_dot_y_couppled_system[2] = -beta_c * x[1]
    x_dot_y_couppled_system[3] = alpha_c * (x[4] - x[3] - fx_prime)
    x_dot_y_couppled_system[4] = x[3] - x[4] + x[5] + sigma * (x[1]- x[4]) # y-coupling
    x_dot_y_couppled_system[5] = -beta_c * x[4]

    return x_dot_y_couppled_system

# z-couppled system dynamics: 1.c
def dynamics_z_coupling(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    x_dot_z_couppled_system = np.zeros(6)

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x_dot_z_couppled_system[0] = alpha_c * (x[1] - x[0] - fx)
    x_dot_z_couppled_system[1] = x[0] - x[1] + x[2]
    x_dot_z_couppled_system[2] = -beta_c * x[1] + sigma * (x[5] - x[2]) # z-coupling
    x_dot_z_couppled_system[3] = alpha_c * (x[4] - x[3] - fx_prime)
    x_dot_z_couppled_system[4] = x[3] - x[4] + x[5]
    x_dot_z_couppled_system[5] = -beta_c * x[4] + sigma * (x[2]- x[5]) # z-coupling

    return x_dot_z_couppled_system

# x-couppled stochastic system dynamics: 1.d
def dynamics_x_coupling_stochastic(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.5):

    # init the 6D vector for the coupled dynamics
    x_dot_x_couppled_system = np.zeros(6)

    # samples from random normal:
    xi = np.random.normal((1, 2))

    # noise intensity:
    D = 0.000001

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # # The 6D vector x-couppled dynamics
    # x_dot_x_couppled_system[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[3] - x[0] + np.sqrt(D)*(xi[0] - xi[1])) # x-coupling and dynamic noise
    # x_dot_x_couppled_system[1] = x[0] - x[1] + x[2]
    # x_dot_x_couppled_system[2] = -beta_c * x[1]
    # x_dot_x_couppled_system[3] = alpha_c * (x[4] - x[3] - fx_prime) + sigma * (x[0]- x[3] + np.sqrt(D)*(xi[0] - xi[1])) # x-coupling and dynamic noise
    # x_dot_x_couppled_system[4] = x[3] - x[4] + x[5]
    # x_dot_x_couppled_system[5] = -beta_c * x[4]

    # The 6D vector x-couppled dynamics
    x_dot_x_couppled_system[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[3] - x[0])
    x_dot_x_couppled_system[1] = x[0] - x[1] + x[2]
    x_dot_x_couppled_system[2] = -beta_c * x[1]
    x_dot_x_couppled_system[3] = alpha_c * (x[4] - x[3] - fx_prime) + sigma * (x[0]- x[3])
    x_dot_x_couppled_system[4] = x[3] - x[4] + x[5]
    x_dot_x_couppled_system[5] = -beta_c * x[4]

    return x_dot_x_couppled_system

# x-couppled pecora-carrol system: 1.e.i
def dynamics_x_driver(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    x_couppled_pecora_carrol_system = np.zeros(6) # x_couppled_pecora_carrol_system = [x, y, z, x', y', z']

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    # fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x_couppled_pecora_carrol_system[0] = alpha_c * (x[1] - x[0] - fx)
    x_couppled_pecora_carrol_system[1] = x[0] - x[1] + x[2]
    x_couppled_pecora_carrol_system[2] = -beta_c * x[1]
    x_couppled_pecora_carrol_system[3] = x[0] - x[4] + x[5]
    x_couppled_pecora_carrol_system[4] = -beta_c * x[4] # -beta * y'

    return x_couppled_pecora_carrol_system

# y-couppled pecora-carrol system: 1.e.ii
def dynamics_y_driver(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    y_couppled_pecora_carrol_system = np.zeros(6) # y_couppled_pecora_carrol_system = [x, y, z, x', y', z']

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    y_couppled_pecora_carrol_system[0] = alpha_c * (x[1] - x[0] - fx)
    y_couppled_pecora_carrol_system[1] = x[0] - x[1] + x[2]
    y_couppled_pecora_carrol_system[2] = -beta_c * x[1]
    y_couppled_pecora_carrol_system[3] = alpha_c * (x[1] - x[3] - fx_prime)
    y_couppled_pecora_carrol_system[4] = -beta_c * x[4]

    return y_couppled_pecora_carrol_system


# z-couppled pecora-carrol system: 1.e.iii
def dynamics_z_driver(x, t, alpha_c = 10.0, beta_c = 14.87, a_c = -1.27, b_c = -0.68, sigma = 0.03):

    # init the 6D vector for the coupled dynamics
    # y_couppled_pecora_carrol_system = np.zeros(6) # y_couppled_pecora_carrol_system = [x, y, z, x', y', z']
    y_couppled_pecora_carrol_system = np.zeros(6) # y_couppled_pecora_carrol_system = [x, y, z, x', y', z']

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    y_couppled_pecora_carrol_system[0] = alpha_c * (x[1] - x[0] - fx)
    y_couppled_pecora_carrol_system[1] = x[0] - x[1] + x[2]
    y_couppled_pecora_carrol_system[2] = -beta_c * x[1]
    y_couppled_pecora_carrol_system[3] = alpha_c * (x[4] - x[3] - fx_prime)
    # y_couppled_pecora_carrol_system[4] = x[3] - x[4] + x[2]
    y_couppled_pecora_carrol_system[4] = -beta_c * x[4]  # duplicated component

    '''
    should we duplicate the z component in the OG and the OG' system, and hence have a 6D vector?
    '''

    return y_couppled_pecora_carrol_system

def error_dynamics(y):
    X = []
    for x in y[:1][0]:
        error = np.zeros(3)
        error[0] = x[0] - x[3]
        error[1] = x[1] - x[4]
        error[2] = x[2] - x[5]
        X.append(error)

    A = np.array(X)

    return A

def component_dynamics(dynamics, index):
    component = []
    for i in dynamics:
        component.append(i[index])
    comp = np.array(component)

    return comp

def init(dynamics):
    # randomly initialize the state vector
    x_init = np.random.uniform(0.1,0.5,3) # make this 6D

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
    y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c), full_output = 1, hmax = 0.01)

    return y, tpoints

def init_system(dynamics, sigma):
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

    return y, tpoints

def init_system_stochastic(dynamics, sigma):
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
    # y = odeint(dynamics, x_init, tpoints, args=(alpha_c, beta_c, a_c, b_c, sigma), full_output = 1, hmax = 0.01)

    # samples from random normal:
    xi = np.random.normal((1, 2))

    # noise intensity:
    D = 0.000001

    x = x_init + t_step*dynamics + np.sqrt(D*t_step)*(xi[0] - xi[1])

    for t in tpoints:
        x = x + t_step*dynamics + np.sqrt(D*t_step)*(xi[0] - xi[1])

    return y, tpoints

def display_dynamics(y, tpoints, title, color, is_errors):
    if is_errors:
        # create the figure; onto which we will plot the above trajectries
        plt.figure(facecolor = 'grey')

        # plt.plot(tpoints, np.array(y[:1][0]),'k', color = "green")
        plt.plot(tpoints, y, color = color, linewidth = 2)
        plt.title(title, fontsize = 12)
        plt.xlabel(r"$time$", fontsize = 12)
        plt.ylabel(r"$y$", fontsize = 12)
        # plt.xlim([0,60])
        # plt.ylim([-15,15])
        plt.xlim([0,100])
        plt.ylim([-10,10])
        plt.tight_layout()
        plt.show()
    else:
        # create the figure; onto which we will plot the above trajectries
        plt.figure(facecolor = 'grey')

        # plt.plot(tpoints, np.array(y[:1][0]),'k', color = "green")
        plt.plot(tpoints, np.array(y[:1][0]), color = color, linewidth = 2)
        plt.title(title, fontsize = 12)
        plt.xlabel(r"$time$", fontsize = 12)
        plt.ylabel(r"$y$", fontsize = 12)
        # plt.xlim([0,60])
        # plt.ylim([-15,15])
        plt.xlim([0,100])
        plt.ylim([-10,10])
        plt.tight_layout()
        plt.show()

def show_final_dynamics(sigma):

    # # The OG system dynamics (uncoupled)
    # y, tpoints = init(dynamics)
    # display_dynamics(y, tpoints, f"System Dynamics for sigma = {sigma}", "blue", False)

    # x-coupled dynamics:
    y, tpoints = init_system(dynamics_x_coupling, sigma = sigma)
    X = error_dynamics(y)
    display_dynamics(y, tpoints, f"x-coupling dynamics for sigma = {sigma}", "blue", is_errors = False)
    display_dynamics(X, tpoints, f"x-coupling error dynamics for sigma = {sigma}", "red", is_errors = True)

    # y-coupled dynamics:
    y, tpoints = init_system(dynamics_y_coupling, sigma = sigma)
    X = error_dynamics(y)
    display_dynamics(y, tpoints, f"y-coupling dynamics for sigma = {sigma}", "lime", is_errors = False)
    display_dynamics(X, tpoints, f"y-coupling error dynamics for sigma = {sigma}", "red", is_errors = True)

    # z-coupled dynamics:
    y, tpoints = init_system(dynamics_z_coupling, sigma = sigma)
    X = error_dynamics(y)
    display_dynamics(y, tpoints, f"z-coupling dynamics for sigma = {sigma}", "purple", is_errors = False)
    display_dynamics(X, tpoints, f"z-coupling error dynamics for sigma = {sigma}", "red", is_errors = True)


if __name__ == '__main__':
    # sigma = 0.19 # only exponential divergence for some rounds anything less than this and no exponential divergence

    sigma = 4.0
    print(f"coupling strength: {sigma}")

    # '''
    # want to plot the error integral against the coupling strength.
    # E_x = (1/T)*(sum from 0 to T of np.abs(x-x_prime))
    # E_y = (1/T)*(sum from 0 to T of np.abs(y-y_prime))
    # E_z = (1/T)*(sum from 0 to T of np.abs(z-z_prime))
    # '''

    # The OG system dynamics (uncoupled)
    # y, tpoints = init(dynamics)
    # display_dynamics(y, tpoints, f"OG System Dynamics for sigma = {sigma}", "blue", False)

    # component-wise dynamics and errors
    show_final_dynamics(sigma)

    # # x-coupled dynamics stochastic:
    # y, tpoints = init_system(dynamics_x_coupling_stochastic, sigma = sigma)
    # X = error_dynamics(y)
    # display_dynamics(y, tpoints, f"x-coupling stohastic dynamics for sigma = {sigma}", "lime", is_errors = False)
    # display_dynamics(X, tpoints, f"x-coupling stohastic error dynamics for sigma = {sigma}", "purple", is_errors = True)

    # # errors vs coupling strength
    # final_x_errors = []
    # final_y_errors = []
    # final_z_errors = []
    #
    # sigmas = []
    #
    # for j in range(1000):
    #
    #     sigmas.append(sigma)
    #
    #     x_errors = []
    #     y_errors = []
    #     z_errors = []
    #
    #     y, tpoints = init_system(dynamics_x_coupling, sigma = sigma)
    #     ex = error_dynamics(y)
    #     for i in ex:
    #         x_errors.append(i[0])
    #
    #     final_x_errors.append(sum(x_errors)/10000)
    #
    #     # y-coupled dynamics:
    #     y, tpoints = init_system(dynamics_y_coupling, sigma = sigma)
    #     ey = error_dynamics(y)
    #     for i in ey:
    #         y_errors.append(i[1])
    #
    #     final_y_errors.append(sum(y_errors)/10000)
    #
    #     # z-coupled dynamics:
    #     y, tpoints = init_system(dynamics_z_coupling, sigma = sigma)
    #     ez = error_dynamics(y)
    #     for i in ez:
    #         z_errors.append(i[2])
    #
    #     final_z_errors.append(sum(z_errors)/10000)
    #
    #     sigma += 0.01
    #
    # # # print(f"final_x_errors: {final_x_errors} ")
    # # # print(f"final_y_errors: {final_y_errors} ")
    # # # print(f"final_z_errors: {final_z_errors} ")
    # # # print(f"sigmas: {sigmas}")
    # #
    # # x-errors
    # plt.plot(sigmas, final_x_errors)
    # plt.title("coupling strength step size: 0.01")
    # plt.xlabel("coupling strength: sigma", fontsize = 12)
    # plt.ylabel("x-component errors", fontsize = 12)
    # plt.show()
    #
    # # y-errors
    # plt.plot(sigmas, final_y_errors)
    # plt.title("coupling strength step size: 0.01")
    # plt.xlabel("coupling strength: sigma", fontsize = 12)
    # plt.ylabel("y-component errors", fontsize = 12)
    # plt.show()
    #
    # # y-errors
    # plt.plot(sigmas, final_z_errors)
    # plt.title("coupling strength step size: 0.01")
    # plt.xlabel("coupling strength: sigma", fontsize = 12)
    # plt.ylabel("z-component errors", fontsize = 12)
    # plt.show()

    # sigma = 5.6
    #
    # # Pecorra-carrol drive response systems:
    # # x-driving dynamics:
    # y, tpoints = init_system(dynamics_x_driver, sigma = sigma)
    # X = error_dynamics(y)
    # display_dynamics(y, tpoints, f"System Dynamics: x-driving coupling strength = {sigma}", "blue", is_errors = False)
    # display_dynamics(X, tpoints, f"Error Dynamics: x-driving coupling strength = {sigma}", "red", is_errors = True)
    #
    # # y-driving dynamics:
    # y, tpoints = init_system(dynamics_y_driver, sigma = sigma)
    # X = error_dynamics(y)
    # display_dynamics(y, tpoints, f"System Dynamics: y-driving coupling strength = {sigma}", "lime", is_errors = False)
    # display_dynamics(X, tpoints, f"Error Dynamics: y-driving coupling strength = {sigma}", "red", is_errors = True)
    #
    # # z-driving dynamics:
    # y, tpoints = init_system(dynamics_z_driver, sigma = sigma)
    # X = error_dynamics(y)
    # display_dynamics(y, tpoints, f"System Dynamics: z-driving coupling strength = {sigma}", "orange", is_errors = False)
    # display_dynamics(X, tpoints, f"Error Dynamics: z-driving coupling strength = {sigma}", "red", is_errors = True)
