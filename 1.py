import sympy as sp

# Define symbols

# weno9
a, b, c, d, e = sp.symbols('a b c d e')

# weno17
# a, b, c, d, e, f, g, h, i = sp.symbols('a b c d e f g h i')


x = sp.Symbol('x')

u = a + b*x + c*x**2 + d*x**3 + e*x**4

# u = a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6 + h*x**7 + i*x**8


# Define equations for cell averages over intervals [0,1], [1,2], ..., [4,5]
eqs = []
for j in range(9):
    avg = sp.integrate(u, (x, j, j+1))
    eqs.append(avg)

# Define the unknowns (cell averages u0 to u4)
# u_vals = sp.symbols('U1 U2 U3 U4 U5')
u_vals = sp.symbols('u0 u1 u2 u3 u4')
# u_vals = sp.symbols('u0 u1 u2 u3 u4 u5 u6 u7 u8')
# Form equations: avg = u_j
equations = [eqs[i] - u_vals[i] for i in range(5)]
# equations = [eqs[i] - u_vals[i] for i in range(9)]

# Solve for a, b, c, d, e
solution = sp.solve(equations, (a, b, c, d, e), simplify=True)
# solution = sp.solve(equations, (a, b, c, d, e, f, g, h, i), simplify=True)



# xv = 4.0 - sp.sqrt(3.0/5.0)/2 + 0.5
xv = 5

# Use the previously solved values for a, b, c, d, e
Q = solution[a] + solution[b]*xv + solution[c]*xv*xv + solution[d]*xv*xv*xv + solution[e]*xv*xv*xv*xv
Qx = solution[b] + 2*xv*solution[c] + 3*xv*xv*solution[d] + 4*xv*xv*xv*solution[e]
Qxx = 2*solution[c] + 6*xv*solution[d] + 12*xv*xv*solution[e]
Qxxx = 6*solution[d] + 24*xv*solution[e]
Qxxxx = 24*solution[e]

# Q = solution[a] + solution[b]*xv + solution[c]*xv**2 + solution[d]*xv**3 + solution[e]*xv**4 + solution[f]*xv**5 + solution[g]*xv**6 + solution[h]*xv**7 + solution[i]*xv**8
# Qx = solution[b] + 2*xv*solution[c] + 3*xv**2*solution[d] + 4*xv**3*solution[e] + 5*xv**4*solution[f] + 6*xv**5*solution[g] + 7*xv**6*solution[h] + 8*xv**7*solution[i]
# Qxx = 2*solution[c] + 6*xv*solution[d] + 12*xv**2*solution[e] + 20*xv**3*solution[f] + 30*xv**4*solution[g] + 42*xv**5*solution[h] + 56*xv**6*solution[i]
# Qxxx = 6*solution[d] + 24*xv*solution[e] + 60*xv**2*solution[f] + 120*xv**3*solution[g] + 210*xv**4*solution[h] + 336*xv**5*solution[i]
# Qxxxx = 24*solution[e] + 120*xv*solution[f] + 360*xv**2*solution[g] + 840*xv**3*solution[h] + 1680*xv**4*solution[i]

# Simplify the expressions
Q_simplified = sp.simplify(Q)
Qx_simplified = sp.simplify(Qx)
Qxx_simplified = sp.simplify(Qxx)
Qxxx_simplified = sp.simplify(Qxxx)
Qxxxx_simplified = sp.simplify(Qxxxx)

Q_simplified, Qx_simplified, Qxx_simplified, Qxxx_simplified, Qxxxx_simplified
print(f"Q: {Q_simplified};")
print(f"Qx: {Qx_simplified};")
print(f"Qxx: {Qxx_simplified};")
print(f"Qxxx: {Qxxx_simplified};")
print(f"Qxxxx: {Qxxxx_simplified};")