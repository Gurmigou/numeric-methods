import numpy as np

np.set_printoptions(precision=5)

A = np.array([[11, 7, 3, 7], [7, 10, -1, 4], [3, -1, 16, -7], [7, 4, -7, 15]])
b = np.array([[2], [2], [-3], [5]])
n = len(A)

def find_S_D(n):
    D = np.zeros([n, n])
    S = np.zeros([n, n])

    for i in range(0, n):
        for j in range(i, n):
            if i == j:
                # find D[i][i]
                sum = 0
                for p in range(0, i):
                    sum += pow(S[p][i], 2) * D[p][p]
                D[i][i] = np.sign(A[i][i] - sum)

                # find S[i][i]
                S[i][i] = np.sqrt(abs(A[i][i] - sum))
            else:
                sum = 0
                for p in range(0, i):
                    sum += S[p][i] * D[p][p] * S[p][j]
                S[i][j] = (A[i][j] - sum) / (D[i][i] * S[i][i])

    return S, D


def S_T_mul_D(S, D):
    S_T = S.transpose()
    return np.matmul(S_T, D)


# A * x = B: x is result
def find_vector(A, B):
    inv = np.linalg.inv(A)
    res = np.matmul(inv, B)
    for i in range(0, len(res)):
        res[i][0] = round(res[i][0], 7)
    return res


A_T = A.transpose()
print("A transpose:")
print(A_T)

# Check if matrix A and its transpose are equal
print("\nCheck A and A transpose equality:")
print(np.equal(A, A_T))

S, D = find_S_D(n)
print("\nMatrix S:")
print(S)

print("\nMatrix D:")
print(D)

S_T_mul_D = S_T_mul_D(S, D)
print("\nMatrix S transpose * matrix D")
print(S_T_mul_D)

Y_vector = find_vector(S_T_mul_D, b)
print("\nY vector:")
print(Y_vector)

X_vector = find_vector(S, Y_vector)
print("\nX vector:")
print(X_vector)

# ----------------------------------------------- Seidel`s method -----------------------------------------------


def check_seidel_method_converges(n):
    A_T = A.transpose()
    print("\nCheck A and A transpose equality:")
    eq = np.equal(A, A_T)

    # check all elements of matrix eq are True
    for i in range(0, n):
        for j in range(0, n):
            if not eq[i][j]:
                print("Seidel's method doesn't converge")
                return False

    def create_matrix(rows, columns):
        matrix = np.zeros([rows, columns])
        for i in range(0, rows):
            for j in range(0, columns):
                matrix[i][j] = A[i][j]
        return matrix

    matrix = A
    while True:
        if not np.linalg.det(matrix) > 0:
            print("Seidel's method doesn't converge")
            return False

        n -= 1
        if n < 2:
            print("Seidel's method converges")
            return True

        matrix = create_matrix(n, n)


converges = check_seidel_method_converges(n)
if not converges:
    exit(1)


def get_constants():
    constants = []
    for i in range(0, n):
        constants.append(A[i][i])
    return constants


def seidels():
    x_0 = [0] * n
    x_1 = [0] * n
    constants = get_constants()
    epsilon = 0.00001

    def calc_x_i(i):
        sum = 0
        for p in range(0, n):
            if p == i:
                continue
            if p > i:
                sum += A[i][p] * x_0[p]
            else:
                sum += A[i][p] * x_1[p]

        return (b[i][0] - sum) / constants[i]

    def calc_norm():
        elements = []
        for i in range(0, n):
            elements.append(x_1[i] - x_0[i])
        return max(elements)

    loop = converges
    while loop:
        for i in range(0, n):
            x_1[i] = calc_x_i(i)

        norm = calc_norm()
        loop = norm > epsilon
        for i in range(0, n):
            x_0[i] = x_1[i]

    for i in range(0, n):
        x_1[i] = round(x_1[i], 7)

    return x_1


x_seidels_solution = seidels()
print(x_seidels_solution)