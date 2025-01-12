import numpy as np
from math import floor, ceil
import time

def simplex(c, A, b):
    c = np.array(c, dtype=float)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)

    m, n = A.shape
    tableau = np.hstack([A, np.eye(m), b.reshape(-1, 1)])
    tableau = np.vstack([tableau, np.hstack([-c, np.zeros(m), 0])])
    
    while True:
        if np.all(tableau[-1, :-1] >= 0):
            break

        pivot_col = np.argmin(tableau[-1, :-1])
        if np.all(tableau[:-1, pivot_col] <= 0): #if all elements in the pivot column are non-positive
            print("Unbounded problem.")
            return None, None

        ratios = [
            tableau[i, -1] / tableau[i, pivot_col]
            if tableau[i, pivot_col] > 0 else float("inf")
            for i in range(m)
        ]
        pivot_row = np.argmin(ratios)
        pivot_element = tableau[pivot_row, pivot_col]

        tableau[pivot_row] /= pivot_element
        for i in range(len(tableau)):
            if i != pivot_row:
                tableau[i] -= tableau[i, pivot_col] * tableau[pivot_row]

    solution = np.zeros(n)
    for i in range(n):
        col = tableau[:, i]
        if np.count_nonzero(col[:-1]) == 1 and np.isclose(np.sum(col[:-1]), 1):
            row = np.where(col[:-1] == 1)[0][0]
            solution[i] = tableau[row, -1]

    return solution, tableau[-1, -1]


def normalize_constraints(A, b, senses):
    """
    Converts constraints into `<=` form for the simplex method.
    Handles `>=` and `=` constraints.
    """
    normalized_A = []
    normalized_b = []

    for row, rhs, sense in zip(A, b, senses):
        if sense == ">=":
            normalized_A.append([-1 * coef for coef in row])
            normalized_b.append(-rhs)
        elif sense == "=":
            normalized_A.append(row)
            normalized_b.append(rhs)
            normalized_A.append([-1 * coef for coef in row])
            normalized_b.append(-rhs)
        else:
            normalized_A.append(row)
            normalized_b.append(rhs)

    return np.array(normalized_A, dtype=float), np.array(normalized_b, dtype=float)

def branch_and_bound(c, A, b, senses, max_depth=50):
    A, b = normalize_constraints(A, b, senses)
    best_solution, best_value = branch((c, A, b), None, float('-inf'), 0, max_depth)
    return best_solution, best_value


def custom_round(x, decimals=10):
    return np.round(x, decimals=decimals).astype(int)


def branch(problem, best_solution, best_value, depth=0, max_depth=50, visited_states=None):
    if visited_states is None:
        visited_states = set()

    if depth > max_depth or (depth > 10 and best_value >= 0):  #early pruning if no improvement; note: 'value' is not defined here yet
        return best_solution, best_value

    c, A, b = problem
    solution, value = simplex(c, A, b) 

    if solution is None:
        return best_solution, best_value

    # Generate state key with custom rounding
    state_key = (
        tuple(custom_round(solution)), 
        tuple(tuple(custom_round(row)) for row in A), 
        tuple(custom_round(b))
    )
    if state_key in visited_states:
        return best_solution, best_value
    visited_states.add(state_key)

    if value <= best_value:
        return best_solution, best_value

    fractional_indices = [i for i, x in enumerate(solution) if not np.isclose(x, round(x), atol=1e-8)]
    if not fractional_indices:
        return solution, value

    branch_var = max(fractional_indices, key=lambda i: abs(solution[i] - round(solution[i])))
    # print(f"Branching on x_{branch_var} = {solution[branch_var]}")

    # Branch on the variable
    lower_bound = floor(solution[branch_var])
    upper_bound = ceil(solution[branch_var])
    # Branch 1: x_branch_var <= lower_bound
    new_row_lower = np.zeros(A.shape[1])
    new_row_lower[branch_var] = 1
    A_lower = np.vstack([A, new_row_lower])
    b_lower = np.append(b, lower_bound)

    # Branch 2: x_branch_var >= upper_bound
    new_row_upper = np.zeros(A.shape[1])
    new_row_upper[branch_var] = -1
    A_upper = np.vstack([A, new_row_upper])
    b_upper = np.append(b, -upper_bound)

    #explore branches recursively, using the same set of visited states to prevent cycles
    best_solution, best_value = branch((c, A_lower, b_lower), best_solution, best_value, depth + 1, max_depth, visited_states)
    best_solution, best_value = branch((c, A_upper, b_upper), best_solution, best_value, depth + 1, max_depth, visited_states)

    return best_solution, best_value


