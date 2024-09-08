import numpy as np
import matplotlib.pyplot as plt

# 1. Solving the damped wave equation using the explicit FDTD method
def explicit_FDTD(T=40, N=100, dt=0.01, c=1.0, kappa=0.001, rho=0.01):
    dx = 1/(N - 1)
    M = int(T/dt)  # Number of time steps

    # Initialize the wavefield
    u = np.zeros((M, N))
    x = np.linspace(0, 1, N)
    
    # Initial conditions: a small perturbation
    u[0, :] = np.sin(np.pi*x)  # Initial displacement
    u[1, :] = u[0, :]  # Assuming the string was at rest initially

    # Time-stepping loop (explicit)
    for n in range(1, M - 1):
        for i in range(1, N - 1):
            u[n+1, i] = (2 - 2*kappa*dt/rho)*u[n, i] \
                        + (kappa*dt/rho - 1)*u[n-1, i] \
                        + (c**2*dt**2/dx**2)*(u[n, i+1] - 2*u[n, i] + u[n, i-1])
    return u, x

# 2. Solving the damped wave equation using the implicit method
def implicit_method(T=40, N=100, dt=0.01, c=1.0, kappa=0.001, rho=0.01):
    dx = 1/(N - 1)
    M = int(T/dt)  # Number of time steps

    # Initialize the wavefield
    u = np.zeros((M, N))
    u_old = np.zeros(N)
    x = np.linspace(0, 1, N)

    # Initial conditions: a small perturbation
    u[0, :] = np.sin(np.pi*x)  # Initial displacement
    u[1, :] = u[0, :]  # Assuming the string was at rest initially

    # Coefficients for the matrix A
    alpha = (1/dt**2) + (kappa/(2*rho*dt)) + (2*c**2/dx**2)
    beta = -c**2/dx**2

    # Time-stepping loop using the implicit method
    for n in range(1, M - 1):
        A = np.diag(alpha*np.ones(N-2)) + np.diag(beta*np.ones(N-3), 1) + np.diag(beta*np.ones(N-3), -1)
        b = (2*u[n, 1:-1]/dt**2) - (u_old[1:-1]/dt**2) + (kappa*u_old[1:-1]/(2*rho*dt))
        u_new = np.linalg.solve(A, b)
        u_old[1:-1] = u[n, 1:-1].copy()
        u[n+1, 1:-1] = u_new

    return u, x

# 3. Solving the telegraph equation using the implicit method
def telegraph_equation(n, T=20, N=100, dt=0.01, c=1.0):
    dx = 1/(N - 1)  # Moved inside the function
    M = int(T/dt)  # Number of time steps

    # Initialize the wavefield
    u = np.zeros((M, N))
    x = np.linspace(0, 1, N)
    u[0, :] = np.sin(n*np.pi*x)  # Initial displacement
    u[1, :] = u[0, :]  # Assuming the string was at rest initially

    # Boundary conditions
    u[:, 0] = 0
    u[:, -1] = 0

    # Coefficients for the matrix A
    alpha = (1/dt**2) + (1/(2*dt)) + (2/dx**2) + 2
    beta = -1/dx**2

    # Time-stepping loop using the implicit method
    for step in range(1, M - 1):
        A = np.diag(alpha*np.ones(N-2)) + np.diag(beta*np.ones(N-3), 1) + np.diag(beta*np.ones(N-3), -1)
        b = (2*u[step, 1:-1]/dt**2) - (u[step-1, 1:-1]/dt**2) + (u[step-1, 1:-1]/(2*dt))
        u_new = np.linalg.solve(A, b)
        u[step+1, 1:-1] = u_new

    return u, x

# Plotting function for comparing explicit and implicit methods
def plot_comparison(explicit_u, implicit_u, x, time_idx, method='explicit'):
    plt.plot(x, explicit_u[time_idx, :], label='Explicit Method', color='r')
    plt.plot(x, implicit_u[time_idx, :], label='Implicit Method', linestyle='--', color='b')
    plt.xlabel('Position along the string')
    plt.ylabel('Displacement')
    plt.title(f'Comparison at t = {time_idx*0.01:.2f} s')
    plt.legend()
    plt.grid(True)
    plt.show()

# Plot results for different values of n for the telegraph equation
def plot_telegraph_modes():
    ns = [1, 2, 3, 4]
    time_values = [0.0, 1.0, 2.99, 8.0, 16.0, 19.99]  # Time points
    time_steps = [int(t/0.01) for t in time_values]

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.ravel()

    for i, n in enumerate(ns):
        u, x = telegraph_equation(n)
        for t in time_steps:
            axs[i].plot(x, u[t, :], label=f't={t*0.01:.2f}s')
        axs[i].set_title(f'Mode n={n}')
        axs[i].set_xlabel('Position along the string')
        axs[i].set_ylabel('Displacement')
        axs[i].legend()
        axs[i].grid(True)

    plt.suptitle("Telegraph Equation for Different Modes (n)")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# Part 1: Explicit Method (FDTD)
explicit_u, x = explicit_FDTD()

# Part 2: Implicit Method
implicit_u, x_implicit = implicit_method()

# Plot comparison at a specific time step
time_idx = 500  # Time index for comparison (500*dt = 5 seconds)
plot_comparison(explicit_u, implicit_u, x, time_idx)

# Part 3: Telegraph Equation for different modes
plot_telegraph_modes()