"""Структуры данных алгоритма исключения переменных (АИП).

Многочлены - элементы кольца QQ[x1..xm] из sympy.polys.rings (разреженные словари
мономов), поэтому вся арифметика точная и в десятки раз быстрее прежней связки
Poly(simplify(sympify(str(expr)))).

Два режима хранения:
  * точный (numeric=False)  - примитивный целочисленный многочлен с положительным
                              старшим коэффициентом (пропорциональные многочлены совпадают);
  * численный (numeric=True) - после подстановки приближённого (иррационального) корня:
                              коэффициенты нормированы на максимальный, округлены до
                              шага 2^-PRECISION_BITS, «шумовые нули» отброшены.
"""
from __future__ import annotations

from sympy import QQ, Symbol, ring as sympy_ring, sympify

import config
from roots import common_real_roots

ZERO = QQ(0)
ONE = QQ(1)


def _round_fixed(c):
    """Округление рационального числа до ближайшего кратного 2^-PRECISION_BITS."""
    bits = config.PRECISION_BITS
    num, den = int(c.numerator), int(c.denominator)
    scaled = (2 * num * (1 << bits) + den) // (2 * den)
    return QQ(scaled, 1 << bits)


def max_abs(elem):
    return max((abs(c) for _, c in elem.terms()), default=ZERO)


class Context:
    """Кольцо многочленов QQ[x1..xm]; первая переменная исключается первой."""

    def __init__(self, names):
        names = [str(n) for n in names]
        if len(set(names)) != len(names):
            raise ValueError("Имена переменных должны быть различны")
        self.names = names
        self.m = len(names)
        self.ring, *gens = sympy_ring(names, QQ)
        self.gens = tuple(gens)

    def parse(self, text) -> "object":
        """Разбор многочлена из строки (^ и ** - степень, десятичные дроби - точные рациональные)."""
        symbols_ = {n: Symbol(n) for n in self.names}
        expr = sympify(str(text), locals=symbols_, rational=True)
        try:
            return self.ring.from_expr(expr)
        except ValueError as e:
            raise ValueError(f"'{text}' не является многочленом от переменных {self.names}") from e

    @staticmethod
    def parse_number(text):
        """Число из строки (целое, десятичная дробь или a/b) -> точное рациональное QQ."""
        value = sympify(str(text), rational=True)
        if not value.is_Rational:
            raise ValueError(f"'{text}' не является рациональным числом")
        return QQ(int(value.p), int(value.q))

    def poly(self, expr, numeric=False, ref=None) -> "Poly_with_degrees":
        if isinstance(expr, (str, int)) or not hasattr(expr, "terms"):
            expr = self.parse(expr)
        return Poly_with_degrees(self, expr, numeric, ref)


def _coeff_in(ctx, elem, i, k):
    """Коэффициент при x_i^k (многочлен от остальных переменных)."""
    d = {}
    for mon, c in elem.terms():
        if mon[i] == k:
            d[mon[:i] + (0,) + mon[i + 1:]] = c
    return ctx.ring.from_dict(d)


class Poly_with_degrees:
    """Многочлен + степени по каждой переменной + главная переменная x_m (первая вошедшая)."""

    __slots__ = ("ctx", "poly", "numeric", "degrees", "x_m_index", "x_m", "n", "_lc_degrees")

    def __init__(self, ctx: Context, elem, numeric: bool = False, ref=None):
        self.ctx = ctx
        self.numeric = numeric
        self._lc_degrees = None
        if numeric:
            elem = self._normalize_numeric(ctx, elem, ONE if ref is None else ref)
        else:
            elem = self._normalize_exact(elem)
        self.poly = elem
        if not elem:
            self.degrees = (0,) * ctx.m
        else:
            self.degrees = tuple(int(d) for d in elem.degrees())
        self.x_m_index = -1
        self.x_m = None
        for i, d in enumerate(self.degrees):
            if d != 0:
                self.x_m_index = i
                self.x_m = ctx.names[i]
                break
        self.n = self.degrees[self.x_m_index] if self.x_m_index != -1 else 0

    # ---- нормировка ----------------------------------------------------
    @staticmethod
    def _normalize_exact(elem):
        if not elem:
            return elem
        _, prim = elem.primitive()
        if prim.LC < 0:
            prim = -prim
        return prim

    @staticmethod
    def _normalize_numeric(ctx, elem, ref):
        terms = list(elem.terms())
        if not terms:
            return elem
        thr = config.CHOP * ref
        kept = [(m, c) for m, c in terms if abs(c) > thr]
        if not kept:
            return ctx.ring.zero
        mx = max(abs(c) for _, c in kept)
        thr_rel = config.CHOP * mx
        out = {}
        for m, c in kept:
            if abs(c) > thr_rel:
                out[m] = _round_fixed(c / mx)
        return ctx.ring.from_dict(out)

    # ---- свойства ------------------------------------------------------
    @property
    def is_zero(self):
        return not self.poly

    def lead_coeff_elem(self):
        """Старший коэффициент по главной переменной (многочлен от остальных переменных)."""
        if self.x_m_index == -1:
            return self.poly
        return _coeff_in(self.ctx, self.poly, self.x_m_index, self.n)

    def lead_coeff(self) -> "Poly_with_degrees":
        return Poly_with_degrees(self.ctx, self.lead_coeff_elem(), self.numeric, max_abs(self.lead_coeff_elem()))

    def lc_degrees(self):
        if self._lc_degrees is None:
            lc = self.lead_coeff_elem()
            self._lc_degrees = tuple(int(d) for d in lc.degrees()) if lc else (0,) * self.ctx.m
        return self._lc_degrees

    def sort_key(self):
        """Константы идут первыми; далее по главной переменной, её степени n1 (минимальная
        положительная степень - p1), степеням старшего коэффициента (Remark 3 статьи)."""
        if self.x_m_index == -1:
            return (-1,)
        return (self.x_m_index, self.n, self.lc_degrees(), self.degrees, len(self.poly))

    def copy(self):
        return Poly_with_degrees(self.ctx, self.poly, self.numeric, ONE)

    def diff(self, index: int) -> "Poly_with_degrees":
        d = self.poly.diff(self.ctx.gens[index])
        return Poly_with_degrees(self.ctx, d, self.numeric, max_abs(self.poly))

    def subs(self, index: int, value, exact: bool = True) -> "Poly_with_degrees":
        """Подстановка x_index = value. exact=False: value - приближение иррационального числа."""
        ctx = self.ctx
        if exact and not self.numeric:
            return Poly_with_degrees(ctx, self.poly.subs(ctx.gens[index], value), False)
        value = QQ(value)
        powers = {0: ONE}
        acc, mag = {}, {}
        avalue = abs(value)
        apowers = {0: ONE}
        for mon, c in self.poly.terms():
            e = mon[index]
            if e not in powers:
                top = max(powers)
                while top < e:
                    top += 1
                    powers[top] = powers[top - 1] * value
                    apowers[top] = apowers[top - 1] * avalue
            key = mon[:index] + (0,) + mon[index + 1:]
            acc[key] = acc.get(key, ZERO) + c * powers[e]
            mag[key] = mag.get(key, ZERO) + abs(c) * apowers[e]
        # коэффициент, который сократился до шума относительно суммы модулей слагаемых, - ноль
        kept = {k: v for k, v in acc.items() if abs(v) > config.CHOP * mag[k]}
        elem = ctx.ring.from_dict(kept)
        return Poly_with_degrees(ctx, elem, True, ZERO)

    def __eq__(self, other):
        return isinstance(other, Poly_with_degrees) and self.poly == other.poly

    def __hash__(self):
        return hash(self.poly)

    def __str__(self):
        if not self.poly:
            return "0 = 0"
        expr = self.poly.as_expr()
        if self.numeric:
            expr = expr.evalf(config.ROUND_PARAMETER + 1)
        return str(expr).replace("**", "^") + " = 0"

    __repr__ = __str__


class Solution:
    """Решение: значения части (или всех) переменных. Значения - рациональные числа QQ."""

    def __init__(self, ctx: Context, values: dict):
        self.ctx = ctx
        self.values = {i: v for i, v in sorted(values.items()) if v is not None}

    @property
    def solution(self):
        return {self.ctx.names[i]: v for i, v in self.values.items()}

    def float_items(self):
        return [(i, int(v.numerator) / int(v.denominator)) for i, v in self.values.items()]

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return "{" + ", ".join(f"{self.ctx.names[i]}: {round(f, config.ROUND_PARAMETER)}"
                               for i, f in self.float_items()) + "}"

    __str__ = __repr__

    def same_as(self, other: "Solution") -> bool:
        if self.values.keys() != other.values.keys():
            return False
        for (_, a), (_, b) in zip(self.float_items(), other.float_items()):
            if abs(a - b) > config.SOLUTION_TOL * (1 + abs(a)):
                return False
        return True

    def __eq__(self, other):
        return isinstance(other, Solution) and self.same_as(other)

    __hash__ = None


class Solutions_set:
    """Множество решений без повторов (повторы определяются с допуском SOLUTION_TOL)."""

    def __init__(self, solutions=()):
        self.items: list[Solution] = []
        for s in solutions:
            self.add(s)

    def add(self, s: Solution):
        if len(s) == 0 or any(s.same_as(t) for t in self.items):
            return
        self.items.append(s)

    def _key(self, s):
        return tuple(s.float_items())

    def __iter__(self):
        return iter(sorted(self.items, key=self._key))

    def __len__(self):
        return len(self.items)

    def __str__(self):
        if not self.items:
            return "Пустое множество {}"
        return "\n".join(str(s) for s in self)


class Poly_system:
    """Система уравнений p_i = 0; polys отсортированы так, что p1 = polys[0] выбирается по статье."""

    def __init__(self, ctx: Context, polys=(), numeric: bool = False):
        self.ctx = ctx
        self.numeric = numeric
        self.polys: list[Poly_with_degrees] = []
        self._seen = set()
        for p in polys:
            self.insert(p)
        self._update_data()

    # ---- служебное -----------------------------------------------------
    def _update_data(self):
        self.polys.sort(key=lambda p: p.sort_key())
        self.l = len(self.polys)
        if self.l == 0:
            self.p1 = None
            self.x_m_index, self.x_m, self.n1, self.k = -1, None, 0, 0
            return
        self.p1 = self.polys[0]
        self.x_m_index = self.p1.x_m_index
        self.x_m = self.p1.x_m
        self.n1 = self.p1.n
        self.k = sum(1 for p in self.polys if p.x_m_index == self.x_m_index)

    def copy(self):
        return Poly_system(self.ctx, self.polys, self.numeric)

    def insert(self, p, ref=None):
        """Добавить многочлен (Poly_with_degrees или элемент кольца); нули и повторы игнорируются."""
        if isinstance(p, Poly_with_degrees):
            if self.numeric and not p.numeric:
                p = Poly_with_degrees(self.ctx, p.poly, True, max_abs(p.poly))
            elif p.numeric and not self.numeric:
                self.numeric = True
                self.polys = [Poly_with_degrees(self.ctx, q.poly, True, max_abs(q.poly)) for q in self.polys]
                self._seen = {q.poly for q in self.polys}
        else:
            p = Poly_with_degrees(self.ctx, p, self.numeric, ref)
        if p.is_zero or p.poly in self._seen:
            return
        self._seen.add(p.poly)
        self.polys.append(p)
        self._update_data()

    def erase(self, p):
        if p.poly not in self._seen:
            return
        self._seen.discard(p.poly)
        self.polys = [q for q in self.polys if q.poly != p.poly]
        self._update_data()

    def has_nonzero_constants(self):
        return any(p.x_m_index == -1 and not p.is_zero for p in self.polys)

    # ---- подстановка и корни ------------------------------------------
    def subs(self, index: int, value, exact: bool = True) -> "Poly_system":
        res = Poly_system(self.ctx, [], self.numeric or not exact)
        for p in self.polys:
            res.insert(p.subs(index, value, exact))
        return res

    def split_univariate(self):
        """Найти многочлены от одной переменной (первой по порядку) и их общие вещественные корни.

        Возвращает (индекс переменной или None, корни, остальные многочлены системы).
        """
        m = self.ctx.m
        main = None
        for p in self.polys:
            if p.x_m_index != -1 and p.degrees.count(0) == m - 1:
                main = p.x_m_index
                break
        if main is None:
            return None, [], self
        uni = [p for p in self.polys if p.x_m_index == main and p.degrees.count(0) == m - 1]
        rest = Poly_system(self.ctx, [p for p in self.polys if p not in uni], self.numeric)
        roots = common_real_roots([p.poly for p in uni], main, self.numeric)
        return main, roots, rest

    def __repr__(self):
        return "\n".join(str(p) for p in self.polys)

    def show_data(self):
        print(f"L = {self.l}")
        print(f"x_m = {self.x_m}")
        print(f"k = {self.k}")
        print(f"p1 = {self.p1}")
        print(f"n1 = {self.n1}")
