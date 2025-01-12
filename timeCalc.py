import time
from branch_and_bound import branch_and_bound
from gomory import solve_with_gomory
import matplotlib.pyplot as plt
import numpy as np


def generate_problem(n, m):
    """Generates a single problem and returns execution times for both methods."""
    A = np.random.randint(5, 30, size=(n, m))
    b = np.random.randint(10, 100, size=n)
    c = np.random.randint(1, 20, size=m)
    constraint_types = ["<="] * n

    # Measure Branch and Bound execution time
    start_time_bb = time.time()
    branch_and_bound(c, A, b, constraint_types)
    end_time_bb = time.time()



    return end_time_bb - start_time_bb


def test_average_time(n, num_tests=50):
    """Tests average execution time for n variables/constraints over num_tests iterations."""
    total_time_bb = 0
    total_time_gom = 0

    for _ in range(num_tests):
        time_bb= generate_problem(n, n)
        total_time_bb += time_bb

    avg_time_bb = (total_time_bb / num_tests) * 1e9  # Convert seconds to nanoseconds

    return avg_time_bb


def main():
    max_problem_size = 10  # Maximum number of variables/constraints to test
    num_tests_per_problem = 20 # Number of tests per problem size

    avg_times_bb = []
    avg_times_gom = [7177989,
                    91698964,
                    78182817,
                    163012573,
                    513077378,
                    486239195,
                    540265242,
                    767060689,
                    533826973,
                    339983128]

    for n in range(1, max_problem_size + 1):
        avg_time_bb = test_average_time(n, num_tests_per_problem)
        print(avg_time_bb)
        avg_times_bb.append(avg_time_bb)

    # Plotting results
    problem_sizes = list(range(1, max_problem_size + 1))
    plt.plot(problem_sizes, avg_times_bb, label="Branch and Bound")
    plt.plot(problem_sizes, avg_times_gom, label="Gomory")
    plt.xlabel('Number of variables/constraints (n x n)')
    plt.ylabel('Average time taken (ns)')
    plt.title('Branch and Bound vs Gomory (Average Execution Time)')
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
