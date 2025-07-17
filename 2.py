import sympy as sp
import math
S = sp.S
def compute_weno_weights(r=5, eval_x=sp.Rational(1, 2)):
    
    
    U = {i: sp.symbols(f"U{i}") for i in range(-(r-1), r)}
    x = sp.Symbol('x')

    
    def reconstruct_stencil(cells):
        coeffs = sp.symbols(f'a0:{r}')  
        P = sum(coeffs[k] * x ** k for k in range(r))
        eqs = [
            sp.integrate(P, (x, j, j + 1)) - U[j]
            for j in cells
        ]
        sol = sp.solve(eqs, coeffs, rational=True)
        return sp.expand(P.subs(sol).subs(x, eval_x))

    stencils = [list(range(-(r-1)+s, s+1)) for s in range(r)]
    Q_stencils = [reconstruct_stencil(st) for st in stencils]

    
    m = 2 * r - 1            
    deg = 2 * r - 2          
    coeffs_central = sp.symbols(f'A0:{deg+1}')
    P_central = sum(coeffs_central[k] * x ** k for k in range(deg + 1))
    eqs_central = [
        sp.integrate(P_central, (x, j, j + 1)) - U[j]
        for j in range(-(r-1), r)
    ]
    sol_central = sp.solve(eqs_central, coeffs_central, rational=True)
    Q_central = sp.expand(P_central.subs(sol_central).subs(x, eval_x))

    
    coeff_indices = list(range(-(r-1), r))
    def coeff_vector(expr):
        return [sp.diff(expr, U[i]) for i in coeff_indices]

    C_central = coeff_vector(Q_central)
    C_stencils = [coeff_vector(q) for q in Q_stencils]

    
    w = sp.symbols(f'w0:{r}')
    eqs_weights = [
        sum(w[j] * C_stencils[j][k] for j in range(r)) - C_central[k]
        for k in range(len(coeff_indices))
    ]
    sol_weights = sp.solve(eqs_weights, w, rational=True)
    return sol_weights

# ---------------- 示例 -----------------
if __name__ == "__main__":
    
    # In cases you do not trust, you can try these
    
    # # r = 5 ⇒ WENO-9, using the same Gauss points as in the original script
    # eval_x = sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2
    # weno9_weights = compute_weno_weights(r=5, eval_x=eval_x)
    # print("WENO-9 (r=5) weights:", weno9_weights)

    # r = 3 ⇒ WENO-5, interface right side
    # eval_x = 1
    # weno5_weights = compute_weno_weights(r=3, eval_x=eval_x)
    # print("WENO-5 (r=3) weights:", weno5_weights)
    
    # r = 3 ⇒ WENO-11, interface right side
    # eval_x = 1
    # weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    # print("WENO-11 (r=6) weights1 sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2:", weno11_weights)
    
    
    # ——————————————————————————————————————weno11 3 gauss point——————————————————————————————
    eval_x = sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights1 sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2:", weno11_weights)
    
    eval_x = sp.Rational(1, 2) - sp.sqrt(sp.Rational(3, 5)) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights2 sp.Rational(1, 2) - sp.sqrt(sp.Rational(3, 5)) / 2:", weno11_weights)
    
    eval_x = 0.0
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights3 0:", weno11_weights)
    
    
    # ——————————————————————————————————————weno11 4 gauss point——————————————————————————————
    sqrt6_5 = sp.sqrt(S(6)/5)
    
    eval_x = -sp.sqrt((3 + 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights1 -sp.sqrt((3 + 2*sqrt6_5)/7) / 2:", weno11_weights)
    
   
    eval_x = -sp.sqrt((3 - 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights2 -sp.sqrt((3 - 2*sqrt6_5)/7) / 2:", weno11_weights)
    
    
    eval_x = sp.sqrt((3 - 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights3 sp.sqrt((3 - 2*sqrt6_5)/7) / 2:", weno11_weights)
    
    eval_x = sp.sqrt((3 + 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=6, eval_x=eval_x)
    print("WENO-11 (r=6) weights4 sp.sqrt((3 + 2*sqrt6_5)/7) / 2:", weno11_weights)
    

    # ——————————————————————————————————————weno13 3 gauss point——————————————————————————————

    eval_x = sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights1 sp.Rational(1, 2) + sp.sqrt(sp.Rational(3, 5)) / 2:", weno11_weights)
    
    eval_x = sp.Rational(1, 2) - sp.sqrt(sp.Rational(3, 5)) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights2 sp.Rational(1, 2) - sp.sqrt(sp.Rational(3, 5)) / 2:", weno11_weights)
    
    eval_x = 0.0
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights3 0:", weno11_weights)
    
    # ——————————————————————————————————————weno13 4 gauss point——————————————————————————————
    sqrt6_5 = sp.sqrt(S(6)/5)
    
    eval_x = -sp.sqrt((3 + 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights1 -sp.sqrt((3 + 2*sqrt6_5)/7) / 2:", weno11_weights)
    
   
    eval_x = -sp.sqrt((3 - 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights2 -sp.sqrt((3 - 2*sqrt6_5)/7) / 2:", weno11_weights)
    
    
    eval_x = sp.sqrt((3 - 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights3 sp.sqrt((3 - 2*sqrt6_5)/7) / 2:", weno11_weights)
    
    eval_x = sp.sqrt((3 + 2*sqrt6_5)/7) / 2
    weno11_weights = compute_weno_weights(r=7, eval_x=eval_x)
    print("WENO-13 (r=7) weights4 sp.sqrt((3 + 2*sqrt6_5)/7) / 2:", weno11_weights)