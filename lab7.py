import math
import matplotlib.pyplot as plt
import numpy as np
 
def f(x):
    f = np.log10(3 * x - 1) + np.exp(2 * x - 1)
    return f
 
def dfdx(x):
    dfdx = 3 / (np.log(10) * (3 * x - 1)) - 2 * np.exp(2 * x - 1)
    return dfdx
 
def d2fdx2(x):
    d2fdx2 = - 9 / (np.log(10) * (3 * x - 1) ** 2) - 4 * np.exp(2 * x - 1)
    return d2fdx2
 
def chord_tangent_method(left, right, epsilon):
    x_c = x_t = 0
    a_prev = a_next = left
    b_prev = b_next = right
    print(f"f({left})f''({left}) > 0 - {(f(left) * d2fdx2(left)) > 0}")
    print(f"f({right})f''({right}) > 0 - {(f(right) * d2fdx2(right)) > 0}")
    if ((f(left) * d2fdx2(left)) > 0):
        print(f'Неподвижный конец отрезка для метода хорд - {left}')
        while (abs(b_next - a_next) >= epsilon):
            a_next = a_prev - f(a_prev) / dfdx(a_prev)
            b_next = b_prev - (f(b_prev) * (b_prev - a_prev)) / (f(b_prev) - f(a_prev))
            a_prev = a_next
            b_prev = b_next
 
        x_c = b_next 
        x_t = a_next
 
    if ((f(right) * d2fdx2(right)) > 0):
        print(f'Неподвижный конец отрезка для метода хорд - {right}')
        while (abs(a_next - b_next) >= epsilon):
            a_next = a_prev - (f(a_prev) * (b_prev - a_prev)) / (f(b_prev) - f(a_prev))
            b_next = b_prev - f(b_prev) / dfdx(b_prev)
            a_prev = a_next
            b_prev = b_next
 
        x_c = a_next 
        x_t = b_next
 
    return x_c, x_t
 
epsilon = 0.0001
 
x_a = np.linspace(0, 2, 200)
plt.plot(x_a, np.log10(3 * x_a - 1), label='y=log10(3 * x - 1)')
plt.plot(x_a, np.exp(2 * x_a - 1), label='y=exp(2 * x - 1)')
 
plt.legend()
 
print(f'Корни уравнения log10(3 * x - 1) + exp(2 * x - 1) = 0 с точностью epsilon = {epsilon}:')
 
print(f'Первый корень:')
x_c, x_t = chord_tangent_method(0.25, 0.75, epsilon)
print(f'Приближение по методу касательных: {x_t}')
print(f'Приближение по методу хорд: {x_c}')
 
print(f'Второй корень:')
x_c, x_t = chord_tangent_method(1.25, 1.5, epsilon)
print(f'Приближение по методу касательных: {x_t}')
print(f'Приближение по методу хорд: {x_c}')
 
plt.show()