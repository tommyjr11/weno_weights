import sympy as sp
import math
# Define symbolic cell averages for cells -4 to 4 (9 cells)
U = {i: sp.symbols(f"U{i}") for i in range(-4, 5)}

# Symbolic spatial variable
x = sp.Symbol('x')

def reconstruct_stencil(cells, eval_x):
    """Return the reconstruction value at eval_x for a given 5-cell stencil."""
    a, b, c, d, e = sp.symbols('a b c d e')
    P = a + b*x + c*x**2 + d*x**3 + e*x**4
    
    eqs = []
    for j in cells:
        avg = sp.integrate(P, (x, j, j + 1))
        eqs.append(avg - U[j])
    
    sol = sp.solve(eqs, (a, b, c, d, e), rational=True)
    P_sub = P.subs(sol)
    return sp.simplify(P_sub.subs(x, eval_x))

def reconstruct_central(eval_x):
    """Return the 9-cell, degree-8 central reconstruction value at eval_x."""
    coeffs = sp.symbols('A0:9')  # A0, A1, ..., A8
    P = sum(coeffs[k] * x**k for k in range(9))
    
    eqs = []
    for j in range(-4, 5):
        avg = sp.integrate(P, (x, j, j + 1))
        eqs.append(avg - U[j])
    
    sol = sp.solve(eqs, coeffs, rational=True)
    P_sub = P.subs(sol)
    return sp.simplify(P_sub.subs(x, eval_x))

# Evaluation point: interface between cell 0 and 1 (x = 0)
eval_x = sp.Rational(0.5 + sp.sqrt(3/5)/2)

# Define the five 5-cell stencils S0 … S4 (left-biased for x = 0⁺)
stencils = [list(range(-4 + s, 1 + s)) for s in range(5)]

# Compute reconstruction expressions for each stencil at eval_x
Q_stencils = [reconstruct_stencil(st, eval_x) for st in stencils]

# Compute the optimal (degree-8) central reconstruction
Q_central = reconstruct_central(eval_x)

# Extract linear coefficients of each expression with respect to U[-4] … U[4]
coeff_indices = list(range(-4, 5))

def coeff_vector(expr):
    return [sp.simplify(sp.diff(expr, U[i])) for i in coeff_indices]

C_central = coeff_vector(Q_central)
C_stencils = [coeff_vector(q) for q in Q_stencils]

# Unknown optimal weights w0 … w4
w = sp.symbols('w0:5')

# Build and solve the linear system  Σ w_j * C_stencil_j = C_central
eqs_weights = []
for k in range(len(coeff_indices)):
    eq = sum(w[j] * C_stencils[j][k] for j in range(5)) - C_central[k]
    eqs_weights.append(eq)

solution_weights = sp.solve(eqs_weights, w, rational=True)

# Display the optimal linear weights for WENO-9 (right boundary, left-biased)
print(solution_weights)
