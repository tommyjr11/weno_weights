import sympy as sp

# Define symbols

# weno9
# a, b, c, d, e = sp.symbols('a b c d e')

# weno11
a, b, c, d, e, f = sp.symbols('a b c d e f')


x = sp.Symbol('x')

# u = a + b*x + c*x**2 + d*x**3 + e*x**4

u = a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5


# Define equations for cell averages over intervals [0,1], [1,2], ..., [4,5]
eqs = []
for j in range(6):
    avg = sp.integrate(u, (x, j, j+1))
    eqs.append(avg)


# u_vals = sp.symbols('u0 u1 u2 u3 u4')
u_vals = sp.symbols('u0 u1 u2 u3 u4 u5')

# Form equations: avg = u_j
# equations = [eqs[i] - u_vals[i] for i in range(5)]
equations = [eqs[i] - u_vals[i] for i in range(6)]


# Solve for a, b, c, d, e
# solution = sp.solve(equations, (a, b, c, d, e), simplify=True)
solution = sp.solve(equations, (a, b, c, d, e,f), simplify=True)




# xv = 4.0 - sp.sqrt(3.0/5.0)/2 + 0.5
# xv = 0+sp.Rational(1, 2) - sp.sqrt(sp.Rational(3, 5)) / 2
xv = 0+sp.Rational(1, 2)

# Use the previously solved values for a, b, c, d, e
# Q = solution[a] + solution[b]*xv + solution[c]*xv*xv + solution[d]*xv*xv*xv + solution[e]*xv*xv*xv*xv
# Qx = solution[b] + 2*xv*solution[c] + 3*xv*xv*solution[d] + 4*xv*xv*xv*solution[e]
# Qxx = 2*solution[c] + 6*xv*solution[d] + 12*xv*xv*solution[e]
# Qxxx = 6*solution[d] + 24*xv*solution[e]
# Qxxxx = 24*solution[e]

Q = solution[a] + solution[b]*xv + solution[c]*xv*xv + solution[d]*xv*xv*xv + solution[e]*xv*xv*xv*xv + solution[f]*xv*xv*xv*xv*xv
Qx = solution[b] + 2*xv*solution[c] + 3*xv*xv*solution[d] + 4*xv*xv*xv*solution[e] + 5*xv*xv*xv*xv*solution[f]
Qxx = 2*solution[c] + 6*xv*solution[d] + 12*xv*xv*solution[e] + 20*xv*xv*xv*solution[f]
Qxxx = 6*solution[d] + 24*xv*solution[e] + 60*xv*xv*solution[f]
Qxxxx = 24*solution[e] + 120*xv*solution[f]
Qxxxxx = 120*solution[f]

# Simplify the expressions
Q_simplified = sp.simplify(Q)
Qx_simplified = sp.simplify(Qx)
Qxx_simplified = sp.simplify(Qxx)
Qxxx_simplified = sp.simplify(Qxxx)
Qxxxx_simplified = sp.simplify(Qxxxx)
Qxxxxx_simplified = sp.simplify(Qxxxxx)


# Q_simplified, Qx_simplified, Qxx_simplified, Qxxx_simplified, Qxxxx_simplified
Q_simplified, Qx_simplified, Qxx_simplified, Qxxx_simplified, Qxxxx_simplified, Qxxxxx_simplified
print(f"Q: {Q_simplified};")
print(f"Qx: {Qx_simplified};")
print(f"Qxx: {Qxx_simplified};")
print(f"Qxxx: {Qxxx_simplified};")
print(f"Qxxxx: {Qxxxx_simplified};")
print(f"Qxxxxx: {Qxxxxx_simplified};")