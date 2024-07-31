import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time

def solve_sparse_matrix(N):
    V_boundary = 100
    alpha = 1
    
    # Define the diagonal matrix that will be repeated
    main_diag = 2*(1+alpha)*np.ones(N)
    off_diag = -1*np.ones(N-1)
    B_matrix = sparse.diags([main_diag, off_diag, off_diag], [0, -1, 1], shape=(N-2, N-2)).toarray()
    B_matrix = sparse.block_diag(([1], B_matrix, [1])).toarray()
    
    # Create identity matrices for tensor product to produce the block matrix
    Identity_full = sparse.eye(N).toarray()
    Identity_partial = sparse.diags([np.zeros(N), -np.ones(N), -np.ones(N)], [0, -1, 1], shape=(N, N)).toarray()
    Identity_adjusted = np.identity(N-2)
    Identity_adjusted = sparse.block_diag(([0], Identity_adjusted, [0])).toarray()
    
    # Perform the tensor products
    A_matrix = sparse.kron(Identity_full, B_matrix).toarray()+sparse.kron(Identity_partial, Identity_adjusted).toarray()
    
    # Include boundary conditions by modifying the resulting matrix
    A_matrix = sparse.block_diag((np.identity(N), A_matrix)).toarray()
    A_matrix[N:2*N, :N] = -1*Identity_adjusted
        
    # Vector of independent terms for the first row
    b_vector = np.zeros((N+1)*N)
    
    # Boundary conditions for the first row
    b_vector[0:N] = V_boundary
    
    # Solve the system of equations
    solution = np.linalg.solve(A_matrix, b_vector)
    solution = solution[:-N].reshape(N, N)
    return solution

def jacobi_method(N):
    V_boundary = 100
    tolerance = 1e-3
    potential_matrix = np.random.rand(N, N)
    
    # Boundary conditions
    potential_matrix[0, :] = V_boundary
    potential_matrix[1:, 0] = 0
    potential_matrix[N-1, 1:] = 0
    potential_matrix[N-1, :] = 0
    
    # Loop until the error is within the desired threshold
    while True:
        max_error = 0
        for i in range(1, N-1):
            for j in range(1, N-1):
                old_value = potential_matrix[i, j]
                potential_matrix[i, j] = 0.25*(potential_matrix[i+1, j]+potential_matrix[i-1, j]+potential_matrix[i, j+1]+potential_matrix[i, j-1])
                max_error = max(max_error, abs(potential_matrix[i, j]-old_value))
        if max_error < tolerance:
            break
    
    return potential_matrix

def solve_sparse_matrix_neumann(N):
    V_boundary = 100
    alpha = 1
    
    main_diag = 2*(1+alpha)*np.ones(N)
    off_diag = -1*np.ones(N-1)
    B_matrix = sparse.diags([main_diag, off_diag, off_diag], [0, -1, 1], shape=(N, N)).toarray()
    B_matrix[0, 0], B_matrix[N-1, N-1] = 1, 1
    
    Identity_full = sparse.eye(N).toarray()
    Identity_partial = sparse.diags([np.zeros(N), -np.ones(N-1), -np.ones(N-1)], [0, -1, 1], shape=(N, N)).toarray()
    Identity_adjusted = np.identity(N-2)
    Identity_adjusted = sparse.block_diag(([0], Identity_adjusted, [0])).toarray()
    
    A_matrix = sparse.kron(Identity_full, B_matrix).toarray()+sparse.kron(Identity_partial, Identity_adjusted).toarray()
    A_matrix = sparse.block_diag((np.identity(N), A_matrix)).toarray()
    A_matrix[N:2*N, :N] = -1*Identity_adjusted 
        
    b_vector = np.zeros((N+1)*N)
    b_vector[0:N] = V_boundary
    
    solution = np.linalg.solve(A_matrix, b_vector)
    solution = solution[:-N].reshape(N, N)
    return solution

N = 50

# Solving with sparse matrix method
solution_sparse = solve_sparse_matrix(N)

# Determine the range of potential values
vmin, vmax = solution_sparse.min(), solution_sparse.max()

plt.figure()
plt.contourf(solution_sparse, alpha=1, cmap=plt.cm.Oranges, vmin=vmin, vmax=vmax)
C = plt.contour(solution_sparse, colors='black', vmin=vmin, vmax=vmax)
plt.clabel(C, inline=1)
plt.title('Solution with Sparse Matrix Method')

plt.figure()
plt.imshow(solution_sparse, cmap=plt.cm.inferno, origin='lower', vmin=vmin, vmax=vmax)
plt.title('Solution with Sparse Matrix Method (imshow)')
plt.colorbar(label='Potential (V)')
plt.show()

# Solving with Jacobi method
solution_jacobi = jacobi_method(N)

# Determine the range of potential values for Jacobi solution
vmin, vmax = solution_jacobi.min(), solution_jacobi.max()

plt.figure()
plt.contourf(solution_jacobi, alpha=1, cmap=plt.cm.Oranges, vmin=vmin, vmax=vmax)
C = plt.contour(solution_jacobi, colors='black', vmin=vmin, vmax=vmax)
plt.clabel(C, inline=1)
plt.title('Solution with Jacobi Method')

plt.figure()
plt.imshow(solution_jacobi, cmap=plt.cm.inferno, origin='lower', vmin=vmin, vmax=vmax)
plt.title('Solution with Jacobi Method (imshow)')
plt.colorbar(label='Potential (V)')
plt.show()

# Solving with Neumann boundary conditions
solution_neumann = solve_sparse_matrix_neumann(N)

# Determine the range of potential values for Neumann solution
vmin, vmax = solution_neumann.min(), solution_neumann.max()

plt.figure()
plt.contourf(solution_neumann, alpha=1, cmap=plt.cm.Oranges, vmin=vmin, vmax=vmax)
C = plt.contour(solution_neumann, colors='black', vmin=vmin, vmax=vmax)
plt.clabel(C, inline=1)
plt.title('Solution with Neumann Boundary Conditions')

plt.figure()
plt.imshow(solution_neumann, cmap=plt.cm.inferno, origin='lower', vmin=vmin, vmax=vmax)
plt.title('Solution with Neumann Boundary Conditions (imshow)')
plt.colorbar(label='Potential (V)')
plt.show()

# Performance comparison
matrix_sizes = np.arange(3, 30, 1)
times_sparse = []
times_jacobi = []
for N in matrix_sizes:
    start_time_sparse = time.time()
    solve_sparse_matrix(N)
    end_time_sparse = time.time()
    
    start_time_jacobi = time.time()
    jacobi_method(N)
    end_time_jacobi = time.time()
    
    times_sparse.append(end_time_sparse-start_time_sparse)
    times_jacobi.append(end_time_jacobi-start_time_jacobi)
    
fig, ax = plt.subplots()
ax.scatter(matrix_sizes, times_jacobi, label='Jacobi Method', color='b')
ax.scatter(matrix_sizes, times_sparse, label='Sparse Matrices', color='g')
ax.tick_params(axis="y", which="minor")
ax.legend()
ax.set_xlabel('Matrix Order')
ax.set_ylabel('Time (s)')
ax.set_title('Comparison of Numerical Efficiency: Sparse vs. Jacobi')
plt.show()
