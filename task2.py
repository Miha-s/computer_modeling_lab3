#!/usr/bin/python3

import sympy as sp


x = sp.Symbol('x')
u = sp.Symbol('u')
w = sp.Symbol('w')
v = sp.Symbol('v')

w_der = -x*u + 2*x
v_der = w
u_der = v

u0 = 0
v0 = 1
w0 = 0

h = 0.1
a = 0
b = 1

def euler_koshi():
    xn = a

    un = u0
    vn = v0
    wn = w0
    un_t = un + h * u_der.subs(v, vn)
    vn_t = vn + h * v_der.subs(w, wn)
    wn_t = wn + h * w_der.subs(x, xn).subs(u, un)

    helper_table = [[xn, un, vn, wn, un_t, vn_t, wn_t]]

    xn = 0
    un = 1
    vn = 2
    wn = 3
    un_t = 4
    vn_t = 5
    wn_t = 6

    while helper_table[-1][xn] + h <= b:
        new_row = [xn, un, vn, wn, un_t, vn_t, wn_t]
        last = helper_table[-1]
        new_row[xn] = last[xn] + h
        new_row[un] = last[un] + h * (u_der.subs(v, last[vn]) + u_der.subs(v, last[un_t]))/2
        new_row[vn] = last[vn] + h * (v_der.subs(w, last[wn]) + v_der.subs(w, last[vn_t]))/2
        new_row[wn] = last[wn] + h * (w_der.subs(x, last[xn]).subs(u, last[un]) + w_der.subs(x, new_row[xn]).subs(u, last[un_t]))
        new_row[un_t] = new_row[un] + h * u_der.subs(v, new_row[vn])
        new_row[vn_t] = new_row[vn] + h * v_der.subs(w, new_row[wn])
        new_row[wn_t] = new_row[wn] + h * w_der.subs(x, new_row[xn]).subs(u, new_row[xn])

        helper_table.append(new_row)
    
    result_table = []
    for el in helper_table:
        result_table.append([el[xn], el[wn]])
    return result_table

res = euler_koshi()

for value in res:
    print("{:.3f}, {:.3f}".format(value[0], value[1]))