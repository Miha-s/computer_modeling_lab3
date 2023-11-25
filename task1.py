#!/usr/bin/python3

import sympy as sp


x = sp.Symbol('x')
y = sp.Symbol('y')
f = 0.25*y+y*(x**2)
initial = -1 
a = 0
b = 1
h = 0.1

def compute_runge_kuta_4(func, initial, a, n, h):
    xn = a
    yn = initial
    result_table = [[xn, yn]]
    current_n = 1
    print("k1, k2, k3, k4, xn, yn")

    while current_n < n:
        temp_x = xn
        temp_y = yn
        k1 = func.subs(x, temp_x).subs(y, temp_y) * h
        temp_x = xn + h/2
        temp_y = yn + k1/2
        k2 = func.subs(x, temp_x).subs(y, temp_y) * h
        temp_x = xn + h/2
        temp_y = yn + k2/2
        k3 = func.subs(x, temp_x).subs(y, temp_y) * h
        temp_x = xn + h
        temp_y = yn + k3
        k4 = func.subs(x, temp_x).subs(y, temp_y) * h
        yn = yn + (1/6)*(k1 + 2*k2+2*k3+k4)
        xn = xn + h
        form = "{:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(k1, k2, k3, k4, xn, yn)
        print(form)
        result_table.append([xn, yn])
        current_n += 1
    
    return result_table

def extrapol_adams_4(y0, yk, yk_1, yk_2, yk_3, h):
    return y0 + (h/24)*(55*yk - 59*yk_1 + 37*yk_2 - 9*yk_3)

def adams(value_table, h, func, n):
    result_table = []
    for value in value_table:
        result_table.append([value[0], value[1], func.subs(x, value[0]).subs(y, value[1])])
    
    i = 0
    while i < n:
        x_val = result_table[-1][0]
        x_val += h
        y_val = extrapol_adams_4(result_table[-1][1], result_table[-1][2], result_table[-2][2], result_table[-3][2], result_table[-4][2], h)
        value = [x_val, y_val, func.subs(x, x_val).subs(y, y_val)]
        result_table.append(value)
        i += 1

    return result_table


res = compute_runge_kuta_4(f, initial, a, 4, h)
res = adams(res, h, f, 7)


for value in res:
    print("{:.3f}, {:.3f}".format(value[0], value[1]))