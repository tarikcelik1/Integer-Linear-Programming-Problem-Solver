import numpy as np

class GomoryCut:
    def __init__(self, A, b, c,constraint_type):
        self.A = np.array(A, dtype=float)
        self.b = np.array(b, dtype=float).reshape(-1, 1)
        self.constraint_type = constraint_type
        self.c = np.array(-np.array(c, dtype=float)).reshape(-1, 1)
        self.precision = 1e-6
        self.pdigits = 6
        for i in range(len(self.constraint_type)):
            if self.constraint_type[i] == ">=":
                self.A[i] = -self.A[i]
                self.b[i] = -self.b[i]


    def primal_simplex(self):
        while np.min(self.tableau[0, 1:]) < -self.precision:
            j = 1
            while j < self.cols and self.tableau[0, j] >= -self.precision:
                j += 1
            if j == self.cols:
                #no entering column found
                return False

            l = 0
            min_ratio = 0
            for i in range(1, self.rows):
                if self.tableau[i, j] > self.precision:
                    ratio = self.tableau[i, 0] / self.tableau[i, j]
                    if (l == 0 or ratio < min_ratio or
                        (self.basis[l - 1] > self.basis[i - 1] and ratio == min_ratio)):
                        l = i
                        min_ratio = ratio

            if l == 0:
                return False

            self.tableau[l] /= self.tableau[l, j]
            for i in range(self.rows):
                if i != l:
                    self.tableau[i] -= self.tableau[l] * self.tableau[i, j]
            self.basis[l - 1] = j
        return True

    def dual_simplex(self):
        while np.min(self.tableau[1:, 0]) < -self.precision:
            l = 1
            while l < self.rows and self.tableau[l, 0] >= -self.precision:
                l += 1
            if l == self.rows:
                #no leaving row found
                return False

            j = 0
            min_ratio = 0
            for i in range(1, self.cols):
                if self.tableau[l, i] < -self.precision:
                    ratio = -self.tableau[0, i] / self.tableau[l, i]
                    if j == 0 or ratio < min_ratio:
                        j = i
                        min_ratio = ratio

            if j == 0:
                return False

            self.tableau[l] /= self.tableau[l, j]
            for i in range(self.rows):
                if i != l:
                    self.tableau[i] -= self.tableau[l] * self.tableau[i, j]
            self.basis[l - 1] = j

        return True

    def solve(self):
        self.m = self.A.shape[0]
        self.n = self.A.shape[1]

        self.rows = self.m + 1
        self.cols = self.m + self.n + 1

        #construct initial tableau
        self.tableau = np.zeros((self.rows, self.cols))
        self.tableau[0, 0] = 0
        self.tableau[0, 1:] = np.concatenate((self.c.T, np.zeros((1, self.m))), axis=1)
        self.tableau[1:, 0] = self.b.flatten()
        self.tableau[1:, 1 : self.n + 1] = self.A
        self.tableau[1:, self.n + 1 :] = np.eye(self.m)

        self.basis = np.arange(self.n + 1, self.n + self.m + 1)

        while True:
            t = self.primal_simplex()
            if not t:
                return None, None

            #check for integrality
            k = 0
            for i in range(1, self.rows):
                val = self.tableau[i, 0]
                frac = val - np.floor(val)
                # Check if fraction close to integer
                if (frac > self.precision and 1 - frac > self.precision):
                    k = i
                    break

            if k == 0:
                #all integer
                solution = np.zeros(self.n)
                for i in range(1, self.rows):
                    if self.basis[i - 1] <= self.n:
                        solution[self.basis[i - 1] - 1] = self.tableau[i, 0]
                optimal_value = -self.tableau[0, 0]
                return solution, -optimal_value

            #add Gomory cut
            self.tableau = np.concatenate((self.tableau, np.zeros(self.rows).reshape(-1, 1)), axis=1)
            self.cols += 1
            self.tableau = np.concatenate((self.tableau, np.zeros(self.cols).reshape(1, -1)), axis=0)
            self.rows += 1

            for i in range(0, self.cols - 1):
                val = self.tableau[k, i]
                frac = val - np.floor(val)
                if frac > self.precision and 1 - frac > self.precision:
                    self.tableau[self.rows - 1, i] = -frac
            self.tableau[self.rows - 1, self.cols - 1] = 1

            self.basis = np.concatenate((self.basis, np.array([self.cols - 1]).reshape(1)))

            t = self.dual_simplex()
            if not t:
                return None, None

        return None, None
    
def solve_with_gomory(A, b, c, constraint_types):
    gomory = GomoryCut(A, b, c, constraint_types)
    return gomory.solve()



