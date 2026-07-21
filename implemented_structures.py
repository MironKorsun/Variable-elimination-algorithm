from sympy import *
from sortedcontainers import SortedSet
from config import variables, round_parameter, close_zero

class Solution:
    def __init__(self, dic):
        self.solution = {var: val for var, val in dic.items() if val is not None}
    
    def __repr__(self):
        return '{' + ', '.join([f"{var}: {round(val, round_parameter)}" for var, val in self.solution.items()]) + '}'
    
    __str__ = __repr__

    def __eq__(self, other):
        if not isinstance(other, Solution):
            return False
        return tuple(self.solution.values()) == tuple(other.solution.values())

    def __hash__(self):
        return tuple(self.solution.values()).__hash__()

class Solutions_set:
    empty_solution = Solution(dict())

    def __init__(self, solutions):
        self.set = SortedSet(solutions, key=self._sorting_key)
        if self.empty_solution in self.set:
            self.set.remove(self.empty_solution)

    def add(self, other):
        if (other in self.set) or (other == self.empty_solution):
            return
        self.set.add(other)

    def _sorting_key(self, s):
        return tuple(s.solution.values())
    
    def __iter__(self):
        return iter(self.set)
    
    def __str__(self):
        if len(self.set) == 0:
            return 'Пустое множество {}'
        return '\n'.join(str(s) for s in self.set)

class Poly_with_degrees: 
    def __init__(self, expr):
        try:
            if abs(float(expr)) < close_zero:
                expr = 0
        except:
            pass
        self.poly = Poly(simplify(sympify(str(expr))), variables)
        self.degrees = tuple(self.poly.degree(v) for v in variables)
        self.x_m_index = -1
        self.x_m = None
        for i in range(len(variables)):
            if isinstance(self.degrees[i], int) and self.degrees[i] != 0:
                self.x_m_index = i
                self.x_m = variables[i]
                break
        self.n = self.degrees[self.x_m_index] if self.x_m_index != -1 else 0

    def lead_coeff(self):
        if self.x_m == None:
            return self
        return Poly_with_degrees(self.poly.as_expr().coeff(self.x_m, self.n))
    
    def copy(self):
        return Poly_with_degrees(self.poly)
    
    def subs(self, *args):
        if len(args) == 2:
            args = {args[0]: args[1]}
        else: 
            args = args[0]
        p = Poly(self.poly.as_expr(), variables)
        for x_i in args:
            p = p.subs(x_i, args[x_i])
        return Poly_with_degrees(p)
        
    def diff(self, x_m):
        return Poly_with_degrees(diff(self.poly, x_m))

    def __eq__(self, other):
        if not isinstance(other, Poly_with_degrees):
            return False
        return self.poly == other.poly
    
    def __hash__(self):
        return self.poly.__hash__()

    def __str__(self):
        return str(Poly(N(self.poly.as_expr(), round_parameter), variables[self.x_m_index])).split(', ')[0].replace('Poly(', '').replace('**', '^') + " = 0"
    
    __repr__ = __str__

    def __mul__(self, other):
        if not isinstance(other, Poly_with_degrees):
            other = Poly_with_degrees(other)
        return Poly_with_degrees(self.poly.as_expr() * other.poly.as_expr())
    
    __rmul__ = __mul__
    
    def __sub__(self, other):
        if not isinstance(other, Poly_with_degrees):
            other = Poly_with_degrees(other)
        return Poly_with_degrees(self.poly - other.poly)
    
    def __itruediv__(self, other):
        if not isinstance(other, Poly_with_degrees):
            other = Poly_with_degrees(other)
        self = Poly_with_degrees(self.poly/other.poly)
        return self

class Poly_system:
    zero_poly = Poly_with_degrees("0")
    def __init__(self, polys):
        self.polys = SortedSet(polys, key=self._sorting_key)
        if self.zero_poly in self.polys:
            self.polys.remove(self.zero_poly)
        if len(self.polys) > 0:
            self._update_data()
        else:
            self._clear_data()

    def _sorting_key(self, p):
        return (p.x_m_index, p.degrees)
    
    def _update_data(self):
        self.l = len(self.polys)
        self.p1 = self.polys[0]
        self.x_m_index = self.p1.x_m_index
        self.x_m = self.p1.x_m
        self.n1 = self.p1.n
        self.k = 0
        for item in self.polys:
            if item.x_m_index == self.x_m_index:
                self.k += 1
            else:
                break

    def _clear_data(self):
        self.l = 0
        self.p1 = None
        self.x_m_index = -1
        self.x_m = None
        self.n1 = 0
        self.k = 0

    def copy(self):
        return Poly_system(self.polys)

    def insert(self, p):
        if isinstance(p, Poly_with_degrees):
            p = p.poly.as_expr()
        elif isinstance(p, Poly):
            p = p.as_expr()
        else:
            p = sympify(p)
        # mx_coeff = max([abs(part.as_coeff_Mul()[0]) for part in p.as_ordered_terms()])
        # if mx_coeff != 0:
        #     p /= mx_coeff
        # p = Add(*[part if abs(part.as_coeff_Mul()[0]) > 0.001 else 0 for part in p.as_ordered_terms()])
        p = Poly_with_degrees(p)
        if p in self.polys or p == self.zero_poly:
            return
        self.polys.add(p)
        self._update_data()

    def erase(self, p: Poly_with_degrees):
        if not (p in self.polys):
            return
        self.polys.remove(p)
        if len(self.polys) > 0:
            self._update_data()
        else:
            self._clear_data()

    def has_nonzero_constants(self):
        for p in self.polys:
            if p.x_m_index > -1:
                return False
            if p.poly.as_expr() != 0:
                return True
        return False

    def solve_system(self) -> list:

        roots = set()
        for p in self.polys:
            cur_roots = real_roots(Poly(p.poly, self.x_m))
            cur_roots = set(round(float(r), 20) for r in cur_roots)
            if len(roots) == 0:
                roots = cur_roots
            else:
                roots &= cur_roots
        return roots

    def eliminate_varibles(self):
        if self.has_nonzero_constants():
            print("Нет решений")
            self.polys.clear()
            self._clear_data()
            return None, []
        system = Poly_system([])
        i = 0
        main_x_m = None
        while i < self.l:
            p = self.polys[i]
            if p.degrees.count(0) == len(p.degrees) - 1:
                if main_x_m == None:
                    main_x_m = p.x_m
                if p.x_m == main_x_m: 
                    system.insert(p)
                    self.erase(p)
                else:
                    i += 1
            else:
                i += 1
        roots = system.solve_system()
        if main_x_m != None and len(roots) == 0:
            print("Нет решений")
        return main_x_m, roots
    
    def subs(self, *args):
        res = Poly_system([])
        for p in self.polys:
            res.insert(p.subs(*args))
        return res

    def __repr__(self):
        return "\n".join([str(p) for p in self.polys]) #  + " '" + str(p.x_m_index) + "' " + str(p.degrees)
    
    def show_data(self):
        print(f"L = {self.l}")
        print(f"x_m = {self.x_m}")
        print(f"k = {self.k}")
        print(f"p1 = {self.p1}")
        print(f"n1 = {self.n1}")