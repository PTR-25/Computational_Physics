import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Section 1: Explicit Method for Heat Diffusion

# Parameters for Section 1
L = 1.0  # Length of the rod
N = 100  # Number of spatial points
D = 1e-2  # Diffusion coefficient
T0 = 100  # Initial temperature
T_boundary = 0  # Boundary temperature
dx = L/N  # Spatial step size
dt = dx**2/(2*D)  # Time step size satisfying the stability condition
t_max = 100  # Maximum time for simulation

# Initialize temperature array for Section 1
T = np.ones(N)*T0
T[0] = T_boundary
T[-1] = T_boundary

# Function to update temperature using the explicit method
def explicit_diffusion(T, D, dx, dt, t_max):
    T_time_evolution = [T.copy()]
    time_steps = int(t_max/dt)
    for _ in range(time_steps):
        T_new = T.copy()
        T_new[1:-1] = T[1:-1] + D*dt/dx**2*(T[2:] - 2*T[1:-1] + T[:-2])
        T = T_new.copy()
        T_time_evolution.append(T.copy())
    return np.array(T_time_evolution)

# Solve using the explicit method
T_time_evolution = explicit_diffusion(T, D, dx, dt, t_max)

# Create the meshgrid for plotting
x = np.linspace(0, L, N)
t = np.linspace(0, t_max, T_time_evolution.shape[0])

# Plot the temperature distribution over time using pcolormesh
fig, ax = plt.subplots()
cbar = ax.pcolormesh(x, t, T_time_evolution, cmap='hot', shading='auto')
ax.set_xlabel('Position along the rod (m)')
ax.set_ylabel('Time (s)')
ax.set_title('Temperature distribution over time using explicit method')
plt.colorbar(cbar, label="Temperature")
plt.show()

# Section 2: Crank-Nicholson Method for Heat Diffusion

# Parameters for Section 2
t_max = 100  # Maximum time for simulation
T_left_boundary = 50  # Temperature at the left boundary
T_right_boundary = 0  # Temperature at the right boundary

# Initialize temperature array for Section 2
T = np.ones(N)*T0
T[0] = T_left_boundary
T[-1] = T_right_boundary

# Create the A and B matrices for the Crank-Nicholson method
alpha = D*dt/(2*dx**2)

A = np.diag((1 + 2*alpha)*np.ones(N)) + np.diag(-alpha*np.ones(N-1), 1) + np.diag(-alpha*np.ones(N-1), -1)
B = np.diag((1 - 2*alpha)*np.ones(N)) + np.diag(alpha*np.ones(N-1), 1) + np.diag(alpha*np.ones(N-1), -1)

# Adjust the boundaries in matrices A and B
A[0, 0], A[0, 1], A[-1, -1], A[-1, -2] = 1, 0, 1, 0
B[0, 0], B[0, 1], B[-1, -1], B[-1, -2] = 1, 0, 1, 0

# Solve using the Crank-Nicholson method
T_time_evolution = [T.copy()]
time_steps = int(t_max/dt)
for _ in range(time_steps):
    T = np.linalg.solve(A, np.dot(B, T))
    T_time_evolution.append(T.copy())

# Convert list to numpy array for visualization
T_time_evolution = np.array(T_time_evolution)

# Create the meshgrid for plotting
x = np.linspace(0, L, N)
t = np.linspace(0, t_max, T_time_evolution.shape[0])

# Plot the temperature distribution over time using pcolormesh
fig, ax = plt.subplots()
cbar = ax.pcolormesh(x, t, T_time_evolution, cmap='hot', shading='auto')
ax.set_xlabel('Position along the rod (m)')
ax.set_ylabel('Time (s)')
ax.set_title('Temperature distribution over time using Crank-Nicholson method')
plt.colorbar(cbar, label="Temperature")
plt.show()

# Section 3: Crank-Nicholson Method with Non-Uniform Diffusion Coefficient

# Parameters for Section 3
D2 = 1e-3  # Diffusion coefficient in the region 0.4 <= x <= 0.6

# Initialize temperature array for Section 3
T = np.ones(N)*T0
T[0] = T_left_boundary
T[-1] = T_right_boundary

# Create the A and B matrices for the Crank-Nicholson method
alpha = D*dt/(2*dx**2)

A = np.diag((1 + 2*alpha)*np.ones(N)) + np.diag(-alpha*np.ones(N-1), 1) + np.diag(-alpha*np.ones(N-1), -1)
B = np.diag((1 - 2*alpha)*np.ones(N)) + np.diag(alpha*np.ones(N-1), 1) + np.diag(alpha*np.ones(N-1), -1)

# Modify A and B matrices for non-uniform diffusion coefficient in the region 0.4 <= x <= 0.6
region_start = int(0.4*N)
region_end = int(0.6*N)
alpha2 = D2*dt/(2*dx**2)

for i in range(region_start, region_end):
    A[i, i] = 1 + 2*alpha2
    if i > 0:
        A[i, i-1] = -alpha2
    if i < N - 1:
        A[i, i+1] = -alpha2
    B[i, i] = 1 - 2*alpha2
    if i > 0:
        B[i, i-1] = alpha2
    if i < N - 1:
        B[i, i+1] = alpha2

# Adjust the boundaries in matrices A and B
A[0, 0], A[0, 1], A[-1, -1], A[-1, -2] = 1, 0, 1, 0
B[0, 0], B[0, 1], B[-1, -1], B[-1, -2] = 1, 0, 1, 0

# Solve using the Crank-Nicholson method
T_time_evolution = [T.copy()]
time_steps = int(t_max/dt)
for _ in range(time_steps):
    T = np.linalg.solve(A, np.dot(B, T))
    T_time_evolution.append(T.copy())

# Convert list to numpy array for visualization
T_time_evolution = np.array(T_time_evolution)

# Create the meshgrid for plotting
x = np.linspace(0, L, N)
t = np.linspace(0, t_max, T_time_evolution.shape[0])

# Plot the temperature distribution over time using pcolormesh
fig, ax = plt.subplots()
cbar = ax.pcolormesh(x, t, T_time_evolution, cmap='hot', shading='auto')
ax.set_xlabel('Position along the rod (m)')
ax.set_ylabel('Time (s)')
ax.set_title('Temperature distribution over time with non-uniform diffusion coefficient')
plt.colorbar(cbar, label="Temperature")
plt.show()
