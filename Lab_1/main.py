import math

epsilon = pow(10, -4)


def func(x):
    return math.sinh(x) - 12 * math.tanh(x) - 0.311


def func_derivative(x):
    return math.cosh(x) - 12 * pow((1 / math.cosh(x)), 2)


def func_derivative_2(x):
    return math.sinh(x) + (24 * math.sinh(x) / pow(math.cosh(x), 3))


def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def dichotomy_method(a, b):
    prev_root = b

    print("0. x = " + str(a) + " f(x) = " + str(func(a)))

    index = 1
    while True:
        cur_root = (a + b) / 2
        cur_root_func = func(cur_root)

        print(str(index) + ". x = " + str(cur_root) + " f(x) = " + str(func(cur_root)))
        index += 1

        if cur_root_func == 0:
            return cur_root_func

        if abs(cur_root - prev_root) < epsilon:
            return cur_root

        # root is left
        if sign(func(a)) == sign(cur_root_func):
            a = cur_root
        else:
            b = cur_root

        prev_root = cur_root


def dichotomy_method_a_priori_estimate(a, b):
    return math.floor(math.log2((b - a) / epsilon)) + 1


def modified_newton_method(x_0):
    func_derivative_x_0 = func_derivative(x_0)
    x_prev = x_0

    print("0. x = " + str(x_0) + " f(x) = " + str(func(x_0)))

    index = 1
    while True:
        x_cur = x_prev - func(x_prev) / func_derivative_x_0

        print(str(index) + ". x = " + str(x_cur) + " f(x) = " + str(func(x_cur)))
        index += 1

        if abs(x_cur - x_prev) < epsilon:
            return x_cur

        x_prev = x_cur


def modified_newton_method_a_priori_estimate(x_abs, q):
    return math.floor(math.log2(math.log(x_abs / epsilon) / math.log(1 / q)) + 1) + 1


def check_derivative_2_not_change_sign(a, b):
    init_sign = sign(func_derivative_2(a))
    while a < b:
        if sign(func_derivative_2(a)) != init_sign:
            return False
        a += 0.001
    return True


def find_m(a, b):
    min_a = abs(func_derivative(a))
    min_b = abs(func_derivative(b))
    return min(min_a, min_b)


def find_M(a, b):
    max_a = abs(func_derivative_2(a))
    max_b = abs(func_derivative_2(b))
    return max(max_a, max_b)


def find_q(m, M, x_abs):
    return (M * x_abs) / (2 * m)


if __name__ == '__main__':
    a = -3.5
    b = -3

    # Dichotomy method
    print("------------ Dichotomy method ------------")
    print("1. Check if f(a) * f(b) < 0: " + str(func(a)) + " * " + str(func(b)) + " = " + str(
        func(a) * func(b)) + ": " + str(func(a) * func(b) < 0))
    print("Dichotomy method result: x = " + str(dichotomy_method(a, b)))
    print("A priori estimate: " + str(dichotomy_method_a_priori_estimate(a, b)))

    # Modified Newton method
    print("\n------------ Modified Newton method ------------")
    print("1. f(a) * f(b) = f(-4) * f(-2) < 0: " + str(func(a)) + " * " + str(func(b)) + " = " + str(
        func(a) * func(b)) + " : " + str(func(a) * func(b) < 0))
    print("2. Check if f''(x) does not change sign in [a; b]: " + str(check_derivative_2_not_change_sign(a, b)))

    x_0_const = -3.25
    print("3. f(x_0) * f''(x_0) > 0: " + str(func(x_0_const)) + " * " + str(func_derivative_2(x_0_const)) + " = " + str(
        func(x_0_const) * func_derivative_2(x_0_const)) + " : " + str(
        func(x_0_const) * func_derivative_2(x_0_const) > 0))

    m = find_m(a, b)
    M = find_M(a, b)

    # |x_0 - x_*|
    # x_0 = -3
    x_abs = 0.25
    q = find_q(m, M, x_abs)
    print("4. Check if q < 1: (" + str(M) + " * " + str(x_abs) + ") / 2 * " + str(m) + " < 1: " + str(q) + " < 1: " + str(q < 1))
    print("Modified Newton method result: x = " + str(modified_newton_method(x_0_const)))
    print("A priori estimate: " + str(modified_newton_method_a_priori_estimate(x_abs, q)))
