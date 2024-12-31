import time
from branch_and_bound_max import branch_and_bound

# Initialize arrays for multiple problems
c_array = [
    [16, 22,12,8],          # Objective coefficients for Problem 1
    [8000,11000,6000,4000], # Objective coefficients for Problem 2
]

A_array = [
    [
        [5, 7,4,3],       # Constraints for Problem 1
    ],
    [
        [5000, 0,0,3000],      # Constraints for Problem 2
        [0, 7000,0,3000],
        [1, 1,0,0],
        [1, 1,1,1],
    ]
    
]

b_array = [
    [14],         # Right-hand side for Problem 1
    [14000,1400,1,2],         # Right-hand side for Problem 2

]

constraint_types_array = [
    ["<="],    # Constraint types for Problem 1
    ["<=", "<=", "<=", "<="]  # Constraint types for Problem 2
]

start_time_array = []   
end_time_array = []

for i in range(len(c_array)):
    start_time_array.append(time.time())
    best_sol, best_val = branch_and_bound(c_array[i], A_array[i], b_array[i], constraint_types_array[i])
    end_time_array.append(time.time())

    print(f"Problem {i + 1}:")
    print("Best solution found:", best_sol)
    print("Best objective value:", best_val)

# Calculate averagae time taken for all problems
avg_time = sum(end_time_array[i] - start_time_array[i] for i in range(len(c_array))) / len(c_array)
print("Average time taken (seconds):", avg_time)
