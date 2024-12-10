import numpy as np
from math import floor, ceil

def simplex(c, A, b):
    """
    Solves a linear programming problem using the simplex algorithm.
    Handles negative coefficients in the RHS and ensures feasibility.
    """
    # Convert inputs to NumPy arrays
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)
    
    m, n = A.shape

    # Create the tableau
    tableau = np.hstack([A, np.eye(m), b.reshape(-1, 1)])
    tableau = np.vstack([tableau, np.hstack([c, np.zeros(m + 1)])])

    # Ensure feasibility by handling negative RHS values
    for i in range(m):
        if tableau[i, -1] < 0:
            tableau[i, :] *= -1  # Flip the row to make the RHS non-negative

    # Pivot until we find the optimal solution
    while True:
        # Check if the current solution is optimal
        if np.all(tableau[-1, :-1] >= 0):
            break

        # Find the pivot column (most negative coefficient in the objective row)
        pivot_col = np.argmin(tableau[-1, :-1])

        # Handle unbounded solutions
        if np.all(tableau[:-1, pivot_col] <= 0):
            raise ValueError("Linear program is unbounded.")

        # Find the pivot row using the minimum ratio test
        ratios = np.divide(
            tableau[:-1, -1], tableau[:-1, pivot_col],
            out=np.full_like(tableau[:-1, -1], np.inf),  # Replace invalid ratios with infinity
            where=tableau[:-1, pivot_col] > 0  # Only divide where pivot column is positive
        )
        pivot_row = np.argmin(ratios)

        # Perform the pivot operation
        pivot_element = tableau[pivot_row, pivot_col]
        tableau[pivot_row] /= pivot_element
        for i in range(m + 1):
            if i != pivot_row:
                tableau[i] -= tableau[i, pivot_col] * tableau[pivot_row]

    # Extract the solution
    solution = np.zeros(n)
    for i in range(n):
        if np.sum(tableau[:, i] == 1) == 1 and np.sum(tableau[:, i] == 0) == m:
            solution[i] = tableau[np.where(tableau[:, i] == 1)[0][0], -1]

    # Optimal value is the last element of the last row
    return solution, tableau[-1, -1]



# Other methods omitted for brevity, including branch_and_bound and gomory_cuts.

def add_dummy_constraints(c, A, b):
    """
    Ensures the number of constraints matches the number of variables
    by adding dummy constraints with all coefficients and RHS set to 0.
    """
    num_variables = len(c)
    num_constraints = len(A)

    if num_constraints < num_variables:
        for _ in range(num_variables - num_constraints):
            A = np.vstack([A, np.zeros(num_variables)])
            b = np.append(b, 0)

    return A, b

def branch(problem, best_solution, best_value, depth=0, max_depth=100):
    c, A, b = problem
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)

    # Stop if max depth is reached
    if depth > max_depth:
        return best_solution, best_value

    # Solve the relaxed problem
    try:
        solution, value = simplex(c, A, b)
    except Exception as e:
        return best_solution, best_value  # Return current best if simplex fails

    # Check if the solution is better than the best so far
    if value > best_value:
        fractional_indices = [i for i, x in enumerate(solution) if not np.isclose(x, round(x), atol=1e-6)]

        if not fractional_indices:
            # All values are integers
            return solution, value

        # Branch on the first fractional variable
        i = fractional_indices[0]
        lower_bound = floor(solution[i])
        upper_bound = ceil(solution[i])

        # Add new constraints for the branches
        new_row_lower = np.zeros(len(c))
        new_row_lower[i] = 1

        new_row_upper = np.zeros(len(c))
        new_row_upper[i] = -1

        # Branch 1: x_i <= floor(solution[i])
        A_lower = np.vstack([A, new_row_lower])
        b_lower = np.append(b, lower_bound)
        best_solution, best_value = branch((c, A_lower, b_lower), best_solution, best_value, depth + 1, max_depth)

        # Branch 2: x_i >= ceil(solution[i])
        A_upper = np.vstack([A, new_row_upper])
        b_upper = np.append(b, -upper_bound)
        best_solution, best_value = branch((c, A_upper, b_upper), best_solution, best_value, depth + 1, max_depth)

    return best_solution, best_value

def branch_and_bound(c, A, b):
    # Ensure the number of constraints equals the number of variables
    A, b = add_dummy_constraints(c, A, b)

    # Call branch-and-bound with adjusted problem
    best_solution, best_value = branch((c, np.array(A), np.array(b)), None, float('-inf'))
    return best_solution, best_value



