import numpy as np
import matplotlib.pyplot as plt

# Central Differences Method
def central_diff_advection(N=100, T=2, dt=0.01, c=1.0, dx=0.1):
    M = int(T/dt)
    x = np.linspace(0, 10, N)
    u = np.zeros((M, N))
    
    u[0, :] = np.exp(-10*(x - 1)**2)

    for j in range(0, M - 1):
        for i in range(1, N - 1):
            u[j+1, i] = u[j, i] - c*dt/(2*dx)*(u[j, i+1] - u[j, i-1])
    
    return u, x, M

# Upwind Scheme
def upwind_advection(N=100, T=2, dt=0.01, c=1.0, dx=0.1):
    M = int(T/dt)
    x = np.linspace(0, 10, N)
    u = np.zeros((M, N))
    
    u[0, :] = np.exp(-10*(x - 1)**2)

    for j in range(0, M - 1):
        for i in range(1, N):
            u[j+1, i] = u[j, i] - c*dt/dx*(u[j, i] - u[j, i-1])

        u[j+1, 0] = 0
        u[j+1, -1] = 0
    
    return u, x, M

# Downwind Scheme
def downwind_advection(N=100, T=2, dt=0.01, c=1.0, dx=0.1):
    M = int(T/dt)
    x = np.linspace(0, 10, N)
    u = np.zeros((M, N))
    
    u[0, :] = np.exp(-10*(x - 1)**2)

    for j in range(0, M - 1):
        for i in range(0, N - 1):
            u[j+1, i] = u[j, i] - c*dt/dx*(u[j, i+1] - u[j, i])

        u[j+1, 0] = 0
        u[j+1, -1] = 0
    
    return u, x, M

# Plotting function for all methods
def plot_advection(u, x, M, dt_val, method_name):
    fig, ax = plt.subplots()
    time_points = [0, int(M // 2), M - 1]
    
    for t in time_points:
        ax.plot(x, u[t, :], label=f"t={round(t*dt_val, 3)}s")
    
    ax.set_title(f"Advection Equation with {method_name}, dt={round(dt_val/dx, 3)}dx")
    ax.set_xlabel("x")
    ax.set_ylabel("Amplitude of the wave")
    ax.legend()
    plt.show()

# Parameters
N = 100
T = 2
dx = 0.1
dt = [1.5*dx, 0.01*dx]
c = 1.0

# Run each method and plot results
for dt_val in dt:
    u, x, M = central_diff_advection(N=N, T=T, dt=dt_val, c=c, dx=dx)
    plot_advection(u, x, M, dt_val, "Central_Differences")
    
    u, x, M = upwind_advection(N=N, T=T, dt=dt_val, c=c, dx=dx)
    plot_advection(u, x, M, dt_val, "Upwind_Scheme")

    u, x, M = downwind_advection(N=N, T=T, dt=dt_val, c=c, dx=dx)
    plot_advection(u, x, M, dt_val, "Downwind_Scheme")