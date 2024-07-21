import numpy as np
import time
import matplotlib.pyplot as plt

def TridiagonalSolver(d, o, u, r):
    n = len(d)
    # Create copies of the input arrays to avoid modifying them
    d = np.copy(d)
    r = np.copy(r)
    o = np.copy(o)
    u = np.copy(u)
    
    # Forward elimination
    for i in range(1, n):
        w = u[i-1] / d[i-1]
        d[i] -= w * o[i-1]
        r[i] -= w * r[i-1]
    
    # Backward substitution
    x = np.zeros(n)
    x[-1] = r[-1] / d[-1]
    for i in range(n-2, -1, -1):
        x[i] = (r[i] - o[i] * x[i+1]) / d[i]
    
    return x

# Function to generate random tridiagonal systems
def generate_tridiagonal_system(n):
    d = np.random.rand(n)
    o = np.random.rand(n-1)
    u = np.random.rand(n-1)
    r = np.random.rand(n)
    return d, o, u, r

# Function to measure time and solve using different methods
def measure_time_and_solve():
    sizes = [10, 100, 500, 1000, 2000, 10000]
    thomas_times = []
    inv_times = []
    solve_times = []
    
    for size in sizes:
        d, o, u, r = generate_tridiagonal_system(size)
        
        # Measure time for TridiagonalSolver
        start = time.time()
        TridiagonalSolver(d, o, u, r)
        thomas_times.append(time.time() - start)
        
        # Construct the full matrix for numpy solutions
        A = np.diag(d) + np.diag(o, k=1) + np.diag(u, k=-1)
        
        # Measure time for numpy.linalg.inv
        start = time.time()
        np.linalg.inv(A)
        inv_times.append(time.time() - start)
        
        # Measure time for numpy.linalg.solve
        start = time.time()
        np.linalg.solve(A, r)
        solve_times.append(time.time() - start)
    
    return sizes, thomas_times, inv_times, solve_times

# Plotting the results
def plot_results(sizes, thomas_times, inv_times, solve_times):
    plt.plot(sizes, thomas_times, label='Thomas Algorithm')
    plt.plot(sizes, inv_times, label='Numpy Inv')
    plt.plot(sizes, solve_times, label='Numpy Solve')
    plt.xlabel('Matrix Size')
    plt.ylabel('Calculation Time (s)')
    plt.legend()
    plt.title('Calculation Time vs Matrix Size')
    plt.show()

# Main execution
sizes, thomas_times, inv_times, solve_times = measure_time_and_solve()
plot_results(sizes, thomas_times, inv_times, solve_times)
