import numpy as np
from math import floor, ceil
from scipy.optimize import linprog

def simplex(c, A, b, constraint_sense):
    # Convert inputs to NumPy arrays
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)
    constraint_sense = np.array(constraint_sense)

    m, n = A.shape

    # Adjust A and b for >= constraints by multiplying by -1
    for i in range(m):
        if constraint_sense[i] == '>=':
            A[i] *= -1
            b[i] *= -1
        elif constraint_sense[i] == '=':
            # Equality constraints can be left as is or handled separately
            pass

    # Create the tableau
    tableau = np.hstack([A, np.eye(m), b.reshape(-1, 1)])
    tableau = np.vstack([tableau, np.hstack([c, np.zeros(m + 1)])])

    # Pivot until we find the optimal solution
    while True:
        # Check if the current solution is optimal
        if np.all(tableau[-1, :-1] >= 0):
            break

        # Find the pivot column
        pivot_col = np.argmin(tableau[-1, :-1])

        # Check for unboundedness
        if np.all(tableau[:-1, pivot_col] <= 0):
            raise Exception("Linear program is unbounded.")

        # Find the pivot row
        ratios = np.divide(
            tableau[:-1, -1], tableau[:-1, pivot_col],
            out=np.full_like(tableau[:-1, -1], np.inf),
            where=tableau[:-1, pivot_col] > 0
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
        col = tableau[:, i]
        if np.count_nonzero(col[:-1]) == 1 and np.isclose(col[-1], 0):
            row = np.where(np.isclose(col[:-1], 1))[0][0]
            solution[i] = tableau[row, -1]

    return solution, tableau[-1, -1]

def branch_and_bound(c, A, b, constraint_sense):
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)
    constraint_sense = np.array(constraint_sense)

    best_solution = None
    best_value = float('-inf')
    nodes = []

    # Initial problem
    nodes.append((A.copy(), b.copy(), constraint_sense.copy()))

    while nodes:
        current_A, current_b, current_constraint_sense = nodes.pop()

        # Solve LP relaxation
        try:
            solution, value = simplex(c, current_A, current_b, current_constraint_sense)
        except Exception as e:
            continue  # Infeasible or unbounded, skip this node

        # Check if solution is better than the best so far
        if value > best_value:
            # Check if solution is integer
            fractional_indices = [i for i, x in enumerate(solution) if not np.isclose(x, round(x), atol=1e-6)]

            if not fractional_indices:
                # All variables are integer
                best_solution = solution
                best_value = value
            else:
                # Branch on the first fractional variable
                i = fractional_indices[0]
                xi = solution[i]
                lower_bound = floor(xi)
                upper_bound = ceil(xi)

                # Create new constraints for branching
                # Left node: x_i <= floor(xi)
                A_left = np.vstack([current_A, np.eye(1, len(c), i)])
                b_left = np.append(current_b, lower_bound)
                constraint_sense_left = np.append(current_constraint_sense, '<=')

                # Right node: x_i >= ceil(xi)
                A_right = np.vstack([current_A, -np.eye(1, len(c), i)])
                b_right = np.append(current_b, -upper_bound)
                constraint_sense_right = np.append(current_constraint_sense, '<=')

                # Add new nodes to the stack
                nodes.append((A_left, b_left, constraint_sense_left))
                nodes.append((A_right, b_right, constraint_sense_right))

    if best_solution is not None:
        return {"Optimal Value": -best_value, "Solution": best_solution}
    else:
        return {"Optimal Value": None, "Solution": None, "Status": "No feasible integer solution found"}


def gomory_cuts(c, A, b, constraint_sense):
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    
    for i in range(len(constraint_sense)):
        if constraint_sense[i] == '<=':
            A_ub.append(A[i])
            b_ub.append(b[i])
        elif constraint_sense[i] == '>=':
            A_ub.append([-a for a in A[i]])
            b_ub.append(-b[i])
        elif constraint_sense[i] == '=':
            A_eq.append(A[i])
            b_eq.append(b[i])

    while True:
        res = linprog(
            c, 
            A_ub=A_ub if A_ub else None, 
            b_ub=b_ub if A_ub else None,
            A_eq=A_eq if A_eq else None, 
            b_eq=b_eq if A_eq else None, 
            method='highs'
        )
        
        if res.success and np.all(np.mod(res.x, 1) == 0):
            print(f"Integer solution found: {res.x}, Optimal value: {res.fun}")
            return res
        elif not res.success:
            print("No solution found.")
            return res
        
        fractional = np.where(np.mod(res.x, 1) != 0)[0]
        if fractional.size == 0:
            break
        
        i = fractional[0]
        cut = res.x[i] - np.floor(res.x[i])
        A_ub.append(np.floor(A[i]))
        b_ub.append(np.floor(b[i]))
    
    print("Gomory cuts method did not converge.")
    return res


def get_user_input():
    num_vars = int(input("Enter the number of variables: "))
    num_constraints = int(input("Enter the number of constraints: "))

    print("\nEnter the coefficients of the objective function:")
    c = [float(input(f"Coefficient for x_{i+1}: ")) for i in range(num_vars)]
    c = [-coef for coef in c]  # Convert to maximization problem

    A = []
    b = []
    constraint_sense = []

    for i in range(num_constraints):
        print(f"\nConstraint {i+1}:")
        row = [float(input(f"Coefficient for x_{j+1}: ")) for j in range(num_vars)]
        sense = input("Enter the constraint sense (<=, >=, =): ").strip()
        rhs = float(input("Enter the right-hand side value: "))

        A.append(row)
        b.append(rhs)
        constraint_sense.append(sense)

    return c, A, b, constraint_sense

def test_methods():
    c, A, b, constraint_sense = get_user_input()

    simplex_solution, simplex_value = simplex(c, A, b, constraint_sense)
    print(f"\nSimplex solution (LP Relaxation): {simplex_solution}, Optimal value: {-simplex_value}")

    print("\nBranch and Bound:")
    result = branch_and_bound(c, A, b, constraint_sense)
    if result["Solution"] is not None:
        print(f"Optimal Value: {result['Optimal Value']}, Solution: {result['Solution']}")
    else:
        print("No feasible integer solution found.")

    print("\nGomory Cuts:")
    gomory_cuts(c, A, b, constraint_sense)
    result2 = gomory_cuts(c, A, b, constraint_sense)
    if result2.success:
        print(f"\nGomory Cuts solution: {result2.x}, Optimal value: {result2.fun}")
    else:
        print("\nGomory Cuts could not find a solution.")

    # Compare with SciPy's implementation
    # Adjust A and b for constraint senses
    A_ub = []
    b_ub = []
    A_eq = []
    b_eq = []

    for i in range(len(constraint_sense)):
        if constraint_sense[i] == '<=':
            A_ub.append(A[i])
            b_ub.append(b[i])
        elif constraint_sense[i] == '>=':
            A_ub.append([-a for a in A[i]])
            b_ub.append(-b[i])
        elif constraint_sense[i] == '=':
            A_eq.append(A[i])
            b_eq.append(b[i])

    res = linprog(c, A_ub=A_ub if A_ub else None, b_ub=b_ub if b_ub else None,
                  A_eq=A_eq if A_eq else None, b_eq=b_eq if b_eq else None, method='highs')

    if res.success:
        print(f"\nSciPy solution: {res.x}, Optimal value: {res.fun}")
    else:
        print("\nSciPy could not find a solution.")

if __name__ == "__main__":
    test_methods()
