import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.integrate import quad

# Set scenario to either 'infinite_well' or 'harmonic_oscillator'
scenario = 'infinite_well'  # Default scenario

# Define parameters
L = 1.0  # Length of the spatial domain
N = 50  # Number of spatial points
dx = L/N  # Spatial step size
x = np.linspace(0, L, N)  # Spatial grid
t_max = 2.0 if scenario == 'infinite_well' else 10.0  # Maximum time
dt = (dx**2)/(4*0.5)  # Time step size ensuring stability
time_steps = int(t_max/dt)  # Number of time steps

# Define initial wave function
def psi(x):
    if scenario == 'infinite_well':
        return 1/np.sqrt(2)*np.sqrt(2/L)*(np.sin(np.pi*x/L) + np.sin(2*np.pi*x/L))
    elif scenario == 'harmonic_oscillator':
        return np.exp(-((x - L/2)**2)/(2*(0.1**2)))

# Normalize the wave function
def Psi_n(x):
    return psi(x)*np.conjugate(psi(x))

psi_n = quad(Psi_n, 0, L)[0]

# Initialize wave function
T = np.zeros((time_steps, N), dtype="complex_")
T[0] = psi(x)/np.sqrt(psi_n)  # Normalized initial wave function

# Potential: infinite well or harmonic oscillator
V = np.zeros(N) if scenario == 'infinite_well' else 0.5*(x - L/2)**2

# Crank-Nicholson matrices
alpha = dt/(2.0*dx**2)
A = np.diag((1 + 2j*alpha)*np.ones(N)) + np.diag(-1j*alpha*np.ones(N-1), 1) + np.diag(-1j*alpha*np.ones(N-1), -1)
B = np.diag((1 - 2j*alpha)*np.ones(N)) + np.diag(1j*alpha*np.ones(N-1), 1) + np.diag(1j*alpha*np.ones(N-1), -1)

# Adjust matrices for infinite well boundary conditions
if scenario == 'infinite_well':
    A[0, 0], A[0, 1], A[-1, -1], A[-1, -2] = 1, 0, 1, 0
    B[0, 0], B[0, 1], B[-1, -1], B[-1, -2] = 1, 0, 1, 0

# Time evolution of the wave function
for i in range(1, time_steps):
    T[i] = np.linalg.solve(A, np.dot(B, T[i-1] - dt*V*T[i-1]))
    if scenario == 'infinite_well':
        T[i][0], T[i][-1] = 0, 0  # Ensure wave function remains zero at the boundaries

# Visualization of the probability density at different time points
time_indices = [0, time_steps//2, time_steps - 1]  # Initial, midpoint, and final time steps
time_labels = [f't={i*dt:.2f}' for i in time_indices]

fig, ax = plt.subplots(figsize=(10, 6))
for i, label in zip(time_indices, time_labels):
    ax.plot(x, np.abs(T[i])**2, label=label)

ax.set_xlim(0, L)
ax.set_ylim(0, 1.5*np.max(np.abs(T)**2))
ax.set_xlabel('Position along the potential')
ax.set_ylabel(r'$|\psi(x, t)|^2$')
ax.set_title(f'Evolution of probability density at different times ({scenario})')
ax.legend()
plt.show()