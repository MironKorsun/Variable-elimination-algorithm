from sympy import *
from implemented_structures import Poly_with_degrees, Poly_system, Solution, Solutions_set
from config import current_values, round_parameter

all_roots = None
initial_system = None

def R(p: Poly_with_degrees, q: Poly_with_degrees, n: int = 1) -> Poly_with_degrees:
    n1 = p.n
    n2 = q.n

    if (n == 0) or (p.x_m_index != q.x_m_index) or (p.x_m_index == -1) or (n1 < n2):
        return p
    
    p_n1 = p.lead_coeff()
    q_n2 = q.lead_coeff()

    g = Poly_with_degrees(gcd(p_n1.poly, q_n2.poly))
    p_n1 /= g
    q_n2 /= g

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

def ves(system: Poly_system) -> bool:
    global all_roots, initial_system
    
    if system.l == 0:
        s = Solution(current_values)
        if s == Solution(dict()):
            print("Нет корней")
            return False
        print(f"Корни: {s}")
        subbed_system = initial_system.subs(s.solution)
        if subbed_system.l != 0:
            if len(s.solution.items()) == len(current_values.items()):
                print("Ложные корни")
                return False
            print("\nПодстановка этих корней в исходную систему:")
            return ves(subbed_system)
        all_roots.add(s)
        return
        
    print(system) 
    x_m, roots = system.eliminate_varibles()

    root_found = False

    for r in roots:
        print()
        current_values[x_m] = r
        print(f"Замена {x_m} = {round(r, round_parameter)}:")
        if ves(system.subs(x_m, r)) == True:
            root_found = True
        current_values[x_m] = None

    if root_found:
        return True

    if x_m == None and system.l != 0:
        if system.l == 1:
            if ves(alpha(system)) == True:
                return True
        else:
            if system.k == 1:
                if ves(alpha(system)) == True:
                    return True
                if ves(beta(system)) == True:
                    return True
            else:
                if ves(gamma(system)) == True:
                    return True
                if ves(delta(system)) == True:
                    return True
        return False