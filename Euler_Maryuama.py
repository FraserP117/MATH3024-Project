import numpy as np
import matplotlib.pyplot as plt

# forward Euler time step
Dt = 0.01

# dynamic noice variance
D = 0.000003

# def init():
#     global x, dynamics, t, timesteps
#     x = 0.1
#     dynamics = [x]
#     t = 0.0
#     timesteps = [t]

def init():
    global x, dynamics, t, timesteps, alpha_c, beta_c, a_c, b_c, sigma
    # assign values to all parameters for the dynamics
    sigma = 5.0
    alpha_c = 10.0; beta_c = 14.87; a_c = -1.27; b_c = -0.68; # sigma = sigma
    x = np.random.uniform(0.1,0.5,6) # make this 6D
    dynamics = [x]
    t = 0.0
    timesteps = [t]

def observe():
    global x, dynamics, t, timesteps, alpha_c, beta_c, a_c, b_c, sigma
    dynamics.append(x)
    timesteps.append(t)

def update():
    global x, dynamics, t, timesteps, alpha_c, beta_c, a_c, b_c, sigma
    xi = np.random.normal((1, 2)) # 2 dynamic noise samples from N(0, 1)

    # define the nonlinear function for each system
    fx = b_c * x[0] + 0.5 * (a_c-b_c) * (np.abs(x[0] + 1.0) - np.abs(x[0] - 1.0))
    fx_prime = b_c * x[3] + 0.5 * (a_c-b_c) * (np.abs(x[3] + 1.0) - np.abs(x[3] - 1.0))

    # The 6D vector x-couppled dynamics
    x[0] = alpha_c * (x[1] - x[0] - fx) + sigma * (x[3] - x[0] + np.sqrt(D)*(xi[0] - xi[1])) # x-coupling
    x[1] = x[0] - x[1] + x[2]
    x[2] = -beta_c * x[1]
    x[3] = alpha_c * (x[4] - x[3] - fx_prime) + sigma * (x[0]- x[3] + np.sqrt(D)*(xi[0] - xi[1])) # x-coupling
    x[4] = x[3] - x[4] + x[5]
    x[5] = -beta_c * x[4]

    t = t + Dt

init()

while t < 100:
    update()
    observe()

fig = plt.figure(figsize = (12,10))
plt.plot(timesteps, dynamics)
plt.xlabel(r'$t$', fontsize = 24)
plt.ylabel(r'$x(t)$', fontsize = 24)
plt.title(r'$\Delta t = 0.015$', fontsize = 24)
plt.show()
