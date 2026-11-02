"""Вещественные корни одномерных многочленов.

Два режима:
  * точный  - коэффициенты рациональные и точные: общие корни = корни НОД
              (с отбрасыванием кратностей, Remark 2 статьи); рациональные корни
              возвращаются точно, иррациональные - с шириной интервала 10^-ROOT_EPS_DIGITS;
  * численный - коэффициенты приближённые (после подстановки иррационального корня):
              берутся корни наименьшего по степени многочлена, а также корни его
              производных (кратный корень после возмущения может «уйти» с вещественной
              оси, а корень производной остаётся рядом), и оставляются те кандидаты,
              для которых невязка всех многочленов пренебрежимо мала.
"""
from __future__ import annotations

from fractions import Fraction
from typing import NamedTuple

from sympy import Poly, QQ as SYM_QQ, Rational, symbols

import config

_T = symbols("t")
_EPS = Rational(1, 10**config.ROOT_EPS_DIGITS)


class Root(NamedTuple):
    value: object  # значение QQ (PythonMPQ / mpq)
    exact: bool    # True - точное рациональное число


def _q(frac: Fraction):
    return config.QQ(frac.numerator, frac.denominator)


def univariate_poly(elem, index: int) -> Poly:
    """Многочлен кольца (зависящий только от переменной с номером index) -> Poly от t над QQ."""
    coeffs = {}
    for mon, c in elem.terms():
        coeffs[mon[index]] = c
    deg = max(coeffs)
    lst = [Rational(int(coeffs[k].numerator), int(coeffs[k].denominator)) if k in coeffs else Rational(0)
           for k in range(deg, -1, -1)]
    return Poly(lst, _T, domain=SYM_QQ)


def _frac(r: Rational) -> Fraction:
    return Fraction(int(r.p), int(r.q))


def _horner(poly: Poly, x: Fraction) -> tuple[Fraction, Fraction]:
    """Значение многочлена в точке x и сумма модулей слагаемых (масштаб для оценки невязки)."""
    value = Fraction(0)
    mag = Fraction(0)
    ax = abs(x)
    for c in poly.all_coeffs():
        fc = _frac(c)
        value = value * x + fc
        mag = mag * ax + abs(fc)
    return value, mag


def _isolate(poly: Poly) -> list[tuple[Fraction, Fraction]]:
    """Изолирующие интервалы вещественных корней (многочлен свободен от квадратов)."""
    if poly.degree() <= 0:
        return []
    intervals = poly.intervals(eps=_EPS, sqf=True)
    out = []
    for item in intervals:
        a, b = item[0] if isinstance(item[0], tuple) else item
        out.append((_frac(Rational(a)), _frac(Rational(b))))
    return out


def _try_rational(poly: Poly, approx: Fraction) -> Fraction | None:
    """Если рядом с approx есть рациональный корень с небольшим знаменателем - вернуть его."""
    cand = approx.limit_denominator(10**9)
    value, _ = _horner(poly, cand)
    return cand if value == 0 else None


def exact_common_roots(polys: list[Poly]) -> list[Root]:
    g = polys[0]
    for q in polys[1:]:
        g = g.gcd(q)
    if g.degree() <= 0:
        return []
    g = g.sqf_part()
    roots = []
    for lo, hi in _isolate(g):
        if lo == hi:
            roots.append(Root(_q(lo), True))
            continue
        mid = (lo + hi) / 2
        rat = _try_rational(g, mid)
        if rat is not None:
            roots.append(Root(_q(rat), True))
        else:
            roots.append(Root(_q(mid), False))
    return roots


def numeric_common_roots(polys: list[Poly]) -> list[Root]:
    base = min((p for p in polys if p.degree() > 0), key=lambda p: p.degree(), default=None)
    if base is None:
        return []
    candidates: list[Fraction] = []
    current = base
    for _ in range(3):  # сам многочлен, производная, вторая производная
        if current.degree() <= 0:
            break
        sqf = current.sqf_part()
        for lo, hi in _isolate(sqf):
            candidates.append((lo + hi) / 2)
        current = current.diff(_T)
    candidates.sort()

    tol = Fraction(config.RESIDUAL_TOL.numerator, config.RESIDUAL_TOL.denominator)
    accepted: list[Fraction] = []
    for r in candidates:
        if accepted and abs(r - accepted[-1]) <= Fraction(1, 10**25) * (1 + abs(r)):
            continue
        ok = True
        for p in polys:
            value, mag = _horner(p, r)
            if abs(value) > tol * mag:
                ok = False
                break
        if ok:
            accepted.append(r)
    return [Root(_q(r), False) for r in accepted]


def common_real_roots(elems: list, index: int, numeric: bool) -> list[Root]:
    """Общие вещественные корни многочленов (элементов кольца), зависящих только от переменной index."""
    polys = [univariate_poly(e, index) for e in elems]
    return numeric_common_roots(polys) if numeric else exact_common_roots(polys)
