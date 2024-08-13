import numpy as np
import matplotlib.pyplot as plt

# Define the Euler Forward method
def euler_forward(n_steps, dt, omega, alpha, x0, v0):
    x = np.zeros(n_steps)
    v = np.zeros(n_steps)
    x[0] = x0
    v[0] = v0
    for i in range(n_steps - 1):
        x[i+1] = x[i] + dt*v[i]
        v[i+1] = v[i] + dt*(-omega**2*x[i] - alpha*v[i])
    return x, v

# Define the Euler Backward method
def euler_backward(n_steps, dt, omega, alpha, x0, v0):
    x = np.zeros(n_steps)
    v = np.zeros(n_steps)
    x[0] = x0
    v[0] = v0
    for i in range(n_steps - 1):
        A = np.array([[1, -dt],
                      [dt*omega**2, 1 + dt*alpha]])
        b = np.array([x[i], v[i]])
        x[i+1], v[i+1] = np.linalg.solve(A, b)
    return x, v

# Define the Crank-Nicholson method
def crank_nicholson(n_steps, dt, omega, alpha, x0, v0):
    x = np.zeros(n_steps)
    v = np.zeros(n_steps)
    x[0] = x0
    v[0] = v0
    for i in range(n_steps - 1):
        A = np.array([[1, -0.5*dt],
                      [0.5*dt*omega**2, 1 + 0.5*dt*alpha]])
        b = np.array([x[i] + 0.5*dt*v[i],
                      v[i] + 0.5*dt*(-omega**2*x[i] - alpha*v[i])])
        x[i+1], v[i+1] = np.linalg.solve(A, b)
    return x, v

# Exact solution for the damped harmonic oscillator
def exact_solution(t, omega, alpha):
    omega_d = np.sqrt(omega**2 - (alpha/2)**2)
    return np.exp(-alpha*t/2)*np.cos(omega_d*t)

# Function to calculate absolute error between numerical and exact solutions
def calculate_error(numerical_solution, exact_solution):
    return np.abs(numerical_solution - exact_solution)

# Function to study error vs. time step
def study_error_vs_timestep(omega, alpha, x0, v0, t_max, timesteps):
    max_errors_f = []
    max_errors_b = []
    max_errors_cn = []
    
    for dt in timesteps:
        n_steps = int(t_max/dt)
        t = np.linspace(0, t_max, n_steps)
        x_exact = exact_solution(t, omega, alpha)
        
        x_f, _ = euler_forward(n_steps, dt, omega, alpha, x0, v0)
        x_b, _ = euler_backward(n_steps, dt, omega, alpha, x0, v0)
        x_cn, _ = crank_nicholson(n_steps, dt, omega, alpha, x0, v0)
        
        max_errors_f.append(np.max(calculate_error(x_f, x_exact)))
        max_errors_b.append(np.max(calculate_error(x_b, x_exact)))
        max_errors_cn.append(np.max(calculate_error(x_cn, x_exact)))
    
    return max_errors_f, max_errors_b, max_errors_cn

# Main execution
if __name__ == "__main__":
    # Parameters
    omega = 1.0  # Natural frequency
    alpha = 0.1  # Damping coefficient
    t_max = 10  # Maximum time
    x0, v0 = 1.0, 0.0  # Initial conditions

    # Solve using each method with a fixed time step
    dt = 0.01  # Time step size
    n_steps = int(t_max/dt)  # Number of time steps
    t = np.linspace(0, t_max, n_steps)

    x_euler_f, v_euler_f = euler_forward(n_steps, dt, omega, alpha, x0, v0)
    x_euler_b, v_euler_b = euler_backward(n_steps, dt, omega, alpha, x0, v0)
    x_crank_n, v_crank_n = crank_nicholson(n_steps, dt, omega, alpha, x0, v0)

    # Exact solution
    x_exact = exact_solution(t, omega, alpha)

    # Plot comparison of the numerical methods with the exact solution
    plt.figure(figsize=(10, 6))
    plt.plot(t, x_exact, label='Exact Solution', color='black', linestyle='--')
    plt.plot(t, x_euler_f, label='Position (Euler Forward)')
    plt.plot(t, x_euler_b, label='Position (Euler Backward)')
    plt.plot(t, x_crank_n, label='Position (Crank-Nicholson)')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.title('Comparison of Numerical Methods with Exact Solution')
    plt.legend()
    plt.show()

    # Plot errors
    error_euler_f = calculate_error(x_euler_f, x_exact)
    error_euler_b = calculate_error(x_euler_b, x_exact)
    error_crank_n = calculate_error(x_crank_n, x_exact)

    plt.figure(figsize=(10, 6))
    plt.plot(t, error_euler_f, label='Error (Euler Forward)')
    plt.plot(t, error_euler_b, label='Error (Euler Backward)')
    plt.plot(t, error_crank_n, label='Error (Crank-Nicholson)')
    plt.xlabel('Time')
    plt.ylabel('Absolute Error')
    plt.title('Error Analysis of Numerical Methods')
    plt.legend()
    plt.show()

    # Study the error vs. time step
    timesteps = np.logspace(-3, -1, 10)
    max_errors_f, max_errors_b, max_errors_cn = study_error_vs_timestep(omega, alpha, x0, v0, t_max, timesteps)

    # Plot error vs. time step
    plt.figure(figsize=(10, 6))
    plt.loglog(timesteps, max_errors_f, label='Euler Forward', marker='o')
    plt.loglog(timesteps, max_errors_b, label='Euler Backward', marker='o')
    plt.loglog(timesteps, max_errors_cn, label='Crank-Nicholson', marker='o')
    plt.xlabel('Time Step \(\Delta t\)')
    plt.ylabel('Maximum Absolute Error')
    plt.title('Error vs. Time Step for Numerical Methods')
    plt.legend()
    plt.show()
