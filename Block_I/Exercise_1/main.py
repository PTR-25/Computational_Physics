import numpy as np
import timeit
import matplotlib.pyplot as plt
import gc
import cProfile
import pstats
import io

def TridiagonalSolver(d, o, u, r):
    n = len(d)
    d = np.copy(d)
    r = np.copy(r)
    o = np.copy(o)
    u = np.copy(u)
    
    # Forward elimination
    for i in range(1, n):
        w = u[i-1]/d[i-1]
        d[i] -= w*o[i-1]
        r[i] -= w*r[i-1]
    
    # Backward substitution
    x = np.zeros(n)
    x[-1] = r[-1]/d[-1]
    for i in range(n-2, -1, -1):
        x[i] = (r[i]-o[i]*x[i+1])/d[i]
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
    sizes = np.arange(10, 2001, 1)
    thomas_times = []
    inv_times = []
    solve_times = []
    
    repetitions = 10
    warmup_reps = 10

    for size in sizes:
        d, o, u, r = generate_tridiagonal_system(size)
        A = np.diag(d) + np.diag(o, k=1) + np.diag(u, k=-1)
        
        # Disable garbage collection
        gc.disable()
        
        # Perform warm-up runs
        for _ in range(warmup_reps):
            TridiagonalSolver(d, o, u, r)
            np.linalg.inv(A)
            np.linalg.solve(A, r)
        
        # Measure time for TridiagonalSolver
        thomas_timer = timeit.Timer(lambda: TridiagonalSolver(d, o, u, r))
        thomas_time = min(thomas_timer.repeat(repeat=5, number=repetitions)) / repetitions
        thomas_times.append(thomas_time)

        # Measure time for numpy.linalg.inv
        inv_timer = timeit.Timer(lambda: np.linalg.inv(A))
        inv_time = min(inv_timer.repeat(repeat=5, number=repetitions)) / repetitions
        inv_times.append(inv_time)
        
        # Measure time for numpy.linalg.solve
        solve_timer = timeit.Timer(lambda: np.linalg.solve(A, r))
        solve_time = min(solve_timer.repeat(repeat=5, number=repetitions)) / repetitions
        solve_times.append(solve_time)
        
        # Enable garbage collection
        gc.enable()
    
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
if __name__ == "__main__":
    sizes, thomas_times, inv_times, solve_times = measure_time_and_solve()
    plot_results(sizes, thomas_times, inv_times, solve_times)