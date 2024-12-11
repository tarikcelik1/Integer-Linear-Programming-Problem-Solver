import numpy as np
from math import floor, ceil

def primal_to_dual(A, b, c):
    """
    Given a primal problem of the form:
        min c^T x
        subject to A x >= b
                   x >= 0

    Returns the A, b, c for the dual problem:
        max b^T y
        subject to A^T y <= c
                   y >= 0
    """
    
    # Ensure A, b, c are numpy arrays
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    c = np.array(c, dtype=float)

    # The dual A is the transpose of the primal A
    dual_A = A.T
    
    # The dual b is the primal c
    dual_b = c
    
    # The dual c is the primal b
    dual_c = b

    return dual_A, dual_b, dual_c

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


def gomory_cuts(c, A, b, max_iterations=100, tol=1e-6):
    """
    Implements the Gomory cutting-plane method for Integer Linear Programming.
    
    Parameters:
        c (list or np.ndarray): The objective function coefficients (for minimization).
        A (list of lists or np.ndarray): The constraint matrix.
        b (list or np.ndarray): The right-hand side vector of constraints.
        max_iterations (int): Maximum number of Gomory cuts to add before giving up.
        tol (float): Tolerance for floating-point comparisons to integrality.
    
    Returns:
        (solution, value): The best integral solution found and its objective value.
                          If no integral solution is found within the iteration limit,
                          returns the best fractional solution found.
    """
    c = np.array(c, dtype=float)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    
    # Try solving the initial LP relaxation
    try:
        solution, value = simplex(c, A, b)
    except Exception as e:
        # If LP is infeasible or unbounded, return immediately
        return None, None
    
    best_solution = solution.copy()
    best_value = value
    
    def is_integral(sol, tolerance=tol):
        return all(abs(x - round(x)) <= tolerance for x in sol)
    
    # Check if the current solution is integral
    if is_integral(solution, tol):
        # Already integral, no cuts needed
        return solution, value
    
    # Iterate and add Gomory cuts
    for iteration in range(max_iterations):
        # Solve LP and get the tableau from simplex. 
        # We'll modify the simplex function slightly or replicate part of its logic 
        # to extract the final tableau. For simplicity, we re-run simplex and reconstruct 
        # a tableau. If your implementation differs, adapt accordingly.
        
        # Rebuild the tableau as in the simplex function:
        m, n = A.shape
        tableau = np.hstack([A, np.eye(m), b.reshape(-1, 1)])
        tableau = np.vstack([tableau, np.hstack([c, np.zeros(m + 1)])])
        
        # Make tableau feasible by handling negative RHS
        for i in range(m):
            if tableau[i, -1] < 0:
                tableau[i, :] *= -1
        
        # Run the simplex pivoting until optimality (we know it's feasible since simplex solved earlier)
        while True:
            if np.all(tableau[-1, :-1] >= -tol):
                # Optimal reached
                break
            pivot_col = np.argmin(tableau[-1, :-1])
            if np.all(tableau[:-1, pivot_col] <= tol):
                # Unbounded
                return best_solution, best_value
            ratios = np.divide(
                tableau[:-1, -1], tableau[:-1, pivot_col],
                out=np.full_like(tableau[:-1, -1], np.inf),
                where=tableau[:-1, pivot_col] > tol
            )
            pivot_row = np.argmin(ratios)
            pivot_element = tableau[pivot_row, pivot_col]
            tableau[pivot_row] /= pivot_element
            for i in range(m + 1):
                if i != pivot_row:
                    tableau[i] -= tableau[i, pivot_col] * tableau[pivot_row]
        
        # Extract current solution
        current_sol = np.zeros(n)
        for i in range(n):
            if np.sum(tableau[:, i] == 1) == 1 and np.sum(tableau[:, i] == 0) == m:
                current_sol[i] = tableau[np.where(tableau[:, i] == 1)[0][0], -1]
        
        current_val = tableau[-1, -1]
        
        # If integral, update best and return
        if is_integral(current_sol, tol):
            if current_val < best_value:
                best_solution = current_sol.copy()
                best_value = current_val
            return best_solution, best_value
        
        # Generate a Gomory cut from one of the fractional rows
        # We'll look at the slack/basic variable rows in the tableau (not the last row).
        # A Gomory cut is generated from a row with a fractional RHS.
        
        # Find a row to generate Gomory cut from
        cut_generated = False
        for row_idx in range(m):
            rhs = tableau[row_idx, -1]
            # Check if RHS is fractional
            frac_part = rhs - floor(rhs)
            if frac_part > tol and frac_part < 1 - tol:
                # Construct the Gomory cut:
                # Cut form: floor(x) >= sum of fractional parts of coefficients.
                # Gomory cut: 
                # sum_{j in N} ((coeff_j) - floor(coeff_j)) * x_j >= fractional_part_of_rhs
                # or equivalently:
                # sum_{j in N} floor(coeff_j)*x_j + x_j*(fractional_part_of_coeff_j) <= floor(rhs)
                # The standard form for a Gomory cut derived from a tableau row (for all non-basic vars):
                # Let row be: x_B = rhs - sum_j coeff_j*x_j (for j in non-basics)
                # The cut: sum_j (frac(coeff_j)) * x_j >= frac(rhs)
                # Move terms around to produce a standard form inequality: 
                # sum_j (-frac(coeff_j))* x_j <= -frac(rhs)
                #
                # We'll build a cut of form: A_new * x >= b_new 
                # where A_new and b_new come from the fractional parts of the coefficients.
                
                row_coeff = tableau[row_idx, :n]  # coefficients associated with original variables
                fractional_coeffs = row_coeff - np.floor(row_coeff)
                fractional_rhs = frac_part
                
                # Gomory cut (in >= form):
                # sum_j fractional_coeffs[j] * x_j >= fractional_rhs
                # Convert to standard form A_new x >= b_new:
                # Just use fractional_coeffs and fractional_rhs directly.
                
                # However, tableau row is in a form with slack variables too. We only need original variables (x), 
                # since slack variables are beyond index n, they form part of a basic solution.
                # Gomory cut involves nonbasic variables. 
                # Since we formed the tableau as [A | I_m | b], original variables are the first n columns.
                
                # Construct the cut:
                A_new = fractional_coeffs
                b_new = fractional_rhs
                
                # Add this cut to A, b:
                A = np.vstack([A, A_new])
                b = np.append(b, b_new)
                print(f"Gomory cut generated from row {row_idx}: {A_new} x >= {b_new}")
                cut_generated = True
                break
        
        if not cut_generated:
            # No suitable fractional row found to generate a Gomory cut.
            # This could mean we are stuck. Return the best known solution.
            return best_solution, best_value
        
    # If we reach here, we hit max_iterations without finding an integral solution
    return best_solution, best_value


