# aip.py
from sympy import *
from implemented_structures import Poly_with_degrees, Poly_system, Solution, Solutions_set
from config import current_values, round_parameter

# Глобальные переменные
all_roots = None
initial_system = None

def solve_system(system: Poly_system) -> list:
    roots = set()
    for p in system.polys:
        cur_roots = set(r.round(round_parameter) for r in real_roots(Poly(p.poly, system.x_m)))
        if len(roots) == 0:
            roots = cur_roots
        else:
            roots &= cur_roots
    return [round(float(r), round_parameter) for r in roots]

def R(p: Poly_with_degrees, q: Poly_with_degrees, n: int = 1) -> Poly_with_degrees:
    n1 = p.n
    n2 = q.n

    if (n == 0) or (p.x_m_index != q.x_m_index) or (p.x_m_index == -1) or (n1 < n2):
        return p
    
    p_n1 = p.lead_coeff()
    q_n2 = q.lead_coeff()

    # g = Poly_with_degrees(gcd(p_n1.poly, q_n2.poly))
    # p_n1 /= g
    # q_n2 /= g

    res = (q_n2 * p) - p.x_m ** (n1 - n2) * (p_n1 * q)
    
    return R(res, q, n-1)

def alpha(system: Poly_system) -> Poly_system:
    print("\nALPHA")
    res = system.copy()
    p1 = res.p1
    dp1_dx_m = p1.diff(p1.x_m)
    res.insert(dp1_dx_m)
    res.insert(p1.n * p1 - p1.x_m * dp1_dx_m)
    res.erase(p1)
    return res

def beta(system: Poly_system) -> Poly_system:
    print("\nBETA")
    res = system.copy()
    res.erase(res.p1)
    return res

def gamma(system: Poly_system) -> Poly_system:
    print("\nGAMMA")
    res = Poly_system([])
    for p in system.polys:
        res.insert(R(p, system.p1, p.n - system.n1 + 1))
    res.insert(system.p1)
    return res

def delta(system: Poly_system) -> Poly_system:
    print("\nDELTA")
    res = system.copy()
    p1 = res.p1
    p1_lead_coeff = p1.lead_coeff()
    res.insert(p1 - p1.x_m ** p1.n * p1_lead_coeff)
    res.insert(p1_lead_coeff)
    res.erase(p1)
    return res

def aip(system: Poly_system):
    global all_roots, initial_system
    
    if system.l != 0:
        print(system) 
    else:
        s = Solution(current_values)
        if s != Solution(dict()):
            print(f"Корни: {s}")
        else:
            print("Нет корней")
            return

    if system.l == 0:
        if initial_system.subs(current_values).l != 0:
            print("Ложные корни")
        else:
            s = Solution(current_values)
            all_roots.add(s)
        return

    x_m, roots = system.eliminate_varibles()

    for r in roots:
        print()
        current_values[x_m] = r
        print(f"Замена {x_m} = {r}:")
        aip(system.subs(x_m, r))

    if x_m == None and system.l != 0:
        if system.l == 1:
            aip(alpha(system))
        else:
            if system.k == 1:
                aip(alpha(system))
                aip(beta(system))
            else:
                aip(gamma(system))
                aip(delta(system))