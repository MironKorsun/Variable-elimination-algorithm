"""Алгоритм исключения переменных (АИП) В.Н. Нефедова.

Каждый шаг заменяет систему p_i = 0 несколькими альтернативными системами с меньшей
степенью по главной переменной x_m (первая переменная, входящая в систему):

    alpha - l = 1:            ∂p1/∂x_m = 0  и  n*p1 - x_m*∂p1/∂x_m = 0
    beta  - k = 1, l >= 2:    многочлен p1 отбрасывается (x_m исключена)
    gamma - k >= 2:           остальные многочлены понижаются по степени остатками R^x_m(p_i, p1)
    delta - k >= 2:           старший коэффициент p1 равен нулю: lc(p1) = 0 и p1 - x_m^n * lc(p1) = 0

Как только в системе появляется многочлен от одной переменной, берутся его вещественные
корни (общие корни всех таких многочленов), корень подставляется в систему и процесс
повторяется. Решения, найденные не для всех переменных, достраиваются подстановкой в
исходную систему; лишние корни отбрасываются.

Все ветви обязательно перебираются (объединение решений альтернативных задач даёт все
решения исходной), поэтому досрочного выхода после первых найденных корней нет.
"""
from __future__ import annotations

from itertools import combinations, product

import config
from implemented_structures import (
    ONE,
    Context,
    Poly_system,
    Poly_with_degrees,
    Solution,
    Solutions_set,
    max_abs,
)
from roots import Root


def R(p: Poly_with_degrees, q: Poly_with_degrees, n: int = 1) -> Poly_with_degrees:
    """n-кратное понижение степени p по x_m с помощью q (псевдоостаток, Remark 1: общий НОД
    старших коэффициентов сокращается).  R(p, q) = lc(q)*p - x^(n1-n2) * lc(p) * q."""
    ctx = p.ctx
    numeric = p.numeric or q.numeric
    while True:
        if n == 0 or p.x_m_index != q.x_m_index or p.x_m_index == -1 or p.n < q.n:
            return p
        g = ctx.gens[p.x_m_index]
        pn, qn = p.lead_coeff_elem(), q.lead_coeff_elem()
        if not numeric:
            d = pn.gcd(qn)
            pn, qn = pn.exquo(d), qn.exquo(d)
        res = qn * p.poly - (g ** (p.n - q.n)) * pn * q.poly
        ref = max_abs(qn) * max_abs(p.poly) + max_abs(pn) * max_abs(q.poly)
        p = Poly_with_degrees(ctx, res, numeric, ref)
        n -= 1


def alpha(system: Poly_system) -> Poly_system:
    res = system.copy()
    p1 = res.p1
    g = system.ctx.gens[p1.x_m_index]
    d = p1.poly.diff(g)
    scale = max_abs(p1.poly) * p1.n
    res.erase(p1)
    res.insert(d, ref=scale)
    res.insert(p1.poly * p1.n - g * d, ref=2 * scale)
    return res


def beta(system: Poly_system) -> Poly_system:
    res = system.copy()
    res.erase(res.p1)
    return res


def gamma(system: Poly_system) -> Poly_system:
    res = Poly_system(system.ctx, [], system.numeric)
    for p in system.polys:
        res.insert(R(p, system.p1, p.n - system.n1 + 1))
    res.insert(system.p1)
    return res


def delta(system: Poly_system) -> Poly_system:
    res = system.copy()
    p1 = res.p1
    g = system.ctx.gens[p1.x_m_index]
    lc = p1.lead_coeff_elem()
    res.erase(p1)
    res.insert(p1.poly - (g ** p1.n) * lc, ref=max_abs(p1.poly))
    res.insert(lc, ref=max_abs(lc))
    return res


def _f(v) -> float:
    return int(v.numerator) / int(v.denominator)


class VES:
    """Решатель системы p_i = 0: все вещественные решения (изолированные) в виде Solutions_set."""

    def __init__(self, ctx: Context, system: Poly_system, verbose: bool = True, out=print):
        self.ctx = ctx
        self.initial = system
        self.verbose = verbose
        self.out = out
        self.all_roots = Solutions_set()
        self.assigned: dict[int, Root] = {}
        self._expanding: set = set()

    def _say(self, text=""):
        if self.verbose:
            self.out(text)

    def run(self, preassigned: dict | None = None) -> Solutions_set:
        """preassigned: значения переменных, уже зафиксированных снаружи (точные, например границы)."""
        self.assigned = {i: Root(v, True) for i, v in (preassigned or {}).items()}
        self._ves(self.initial)
        return self.all_roots

    # ---- рекурсия ------------------------------------------------------
    def _ves(self, system: Poly_system):
        if system.l == 0:
            self._leaf()
            return

        self._say(system)
        if system.has_nonzero_constants():
            self._say("No solutions")
            return

        main, roots, rest = system.split_univariate()
        if main is not None:
            if not roots:
                self._say("No solutions")
                return
            for r in roots:
                self._say()
                self.assigned[main] = r
                self._say(f"Substituting {self.ctx.names[main]} = {round(_f(r.value), config.ROUND_PARAMETER)}:")
                self._ves(rest.subs(main, r.value, r.exact))
                del self.assigned[main]
            return

        if system.l == 1:
            self._branch("ALPHA", alpha(system))
        elif system.k == 1:
            self._branch("ALPHA", alpha(system))
            self._branch("BETA", beta(system))
        else:
            self._branch("GAMMA", gamma(system))
            self._branch("DELTA", delta(system))

    def _branch(self, name: str, system: Poly_system):
        self._say(f"\n{name}")
        self._ves(system)

    def _leaf(self):
        s = Solution(self.ctx, {i: r.value for i, r in self.assigned.items()})
        if len(s) == 0:
            self._say("No roots")
            return
        self._say(f"Roots: {s}")
        sub = self.initial
        for i, r in self.assigned.items():
            sub = sub.subs(i, r.value, r.exact)
        if sub.l != 0:
            if len(s) == self.ctx.m:
                self._say("Extraneous roots")
                return
            key = tuple(s.float_items())
            if key in self._expanding:
                return
            self._say("\nSubstituting these roots into the original system:")
            self._expanding.add(key)
            self._ves(sub)
            self._expanding.discard(key)
            return
        self.all_roots.add(s)


def solve_system(ctx: Context, polynomials, verbose: bool = True, out=print) -> Solutions_set:
    """polynomials: строки или элементы кольца ctx.ring."""
    system = Poly_system(ctx, [ctx.parse(p) if isinstance(p, str) else p for p in polynomials])
    return VES(ctx, system, verbose, out).run()


# ---- глобальный минимум на параллелепипеде ---------------------------------
def _evaluate(elem, point: dict):
    total = 0 * ONE
    for mon, c in elem.terms():
        term = c
        for i, e in enumerate(mon):
            if e:
                term = term * point[i] ** e
        total = total + term
    return total


class MinimumResult:
    def __init__(self):
        self.candidates: list = []  # (Solution с полным набором переменных, значение, свободные переменные)
        self.value = None
        self.points: list[Solution] = []


def global_minimum(ctx: Context, p, bounds: dict, verbose: bool = True, out=print) -> tuple[MinimumResult, Solutions_set]:
    """Глобальный минимум многочлена p на параллелепипеде bounds = {индекс: (нижняя, верхняя)}.

    Критические точки ищутся внутри параллелепипеда, на всех гранях, рёбрах и в вершинах:
    для каждого подмножества переменных, зафиксированных на границах, решается система
    градиента по остальным переменным.
    """
    m = ctx.m
    all_roots = Solutions_set()

    def run(subst: dict):
        ps = p
        for i, v in subst.items():
            ps = ps.subs(ctx.gens[i], v)
        system = Poly_system(ctx, [ps.diff(g) for g in ctx.gens])
        for s in VES(ctx, system, verbose, out).run(preassigned=subst):
            all_roots.add(s)

    if verbose:
        out("\nWithout considering bounds:")
    run({})

    for k in range(1, m + 1):
        for keys in combinations(range(m), k):
            for combo in product(*[bounds[i] for i in keys]):
                subst = dict(zip(keys, combo))
                if verbose:
                    shown = ", ".join(f"{ctx.names[i]}: {_f(v)}" for i, v in subst.items())
                    out(f"\nBound {{{shown}}}:")
                run(subst)

    result = MinimumResult()
    if verbose:
        out("\nCritical points in the domain:")
    for s in all_roots:
        point, free = {}, []
        for i in range(m):
            if i in s.values:
                point[i] = s.values[i]
            else:
                point[i] = bounds[i][0]
                free.append(ctx.names[i])
        inside = all(
            _f(bounds[i][0]) - 1e-12 * (1 + abs(_f(bounds[i][0]))) <= _f(point[i])
            <= _f(bounds[i][1]) + 1e-12 * (1 + abs(_f(bounds[i][1])))
            for i in range(m)
        )
        if not inside:
            continue
        value = _f(_evaluate(p, point))
        full = Solution(ctx, point)
        result.candidates.append((full, value, free))
        if verbose:
            note = f" (free variables set to lower bound: {', '.join(free)})" if free else ""
            out(f"Point #{len(result.candidates)}: {full} Value: {round(value, config.ROUND_PARAMETER)}{note}")

    for full, value, _ in result.candidates:
        if result.value is None or value < result.value - config.SOLUTION_TOL * (1 + abs(result.value)):
            result.value, result.points = value, [full]
        elif abs(value - result.value) <= config.SOLUTION_TOL * (1 + abs(result.value)):
            result.points.append(full)
    return result, all_roots
