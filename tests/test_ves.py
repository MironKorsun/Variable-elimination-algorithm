"""Тесты алгоритма исключения переменных (запуск: pytest)."""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

import pytest
import sympy as sp

from VES import global_minimum, solve_system
from implemented_structures import Context, Poly_system, Poly_with_degrees, Solutions_set

ROOT = pathlib.Path(__file__).resolve().parents[1]


def floats(solutions):
    return sorted(tuple(round(f, 9) for _, f in s.float_items()) for s in solutions)


def residual_ok(ctx, polys, solutions, tol=1e-12):
    """Каждое найденное решение обращает все уравнения в ноль (проверка в точной арифметике)."""
    syms = [sp.Symbol(n) for n in ctx.names]
    exprs = [sp.sympify(p.replace("^", "**")) for p in polys]
    for s in solutions:
        assert len(s) == ctx.m
        point = {sym: sp.Rational(int(v.numerator), int(v.denominator)) for sym, (_, v) in zip(syms, s.values.items())}
        for e in exprs:
            assert abs(float(e.subs(point))) < tol, (e, s)


# ---- разбор и структуры ------------------------------------------------------
def test_parse_is_exact():
    ctx = Context(["x", "y"])
    p = ctx.parse("0.5*x^2 + y**2 - 1/3")
    assert str(p.as_expr()) == "x**2/2 + y**2 - 1/3"


def test_parse_rejects_unknown_symbol():
    with pytest.raises(ValueError):
        Context(["x"]).parse("x + q")


def test_proportional_polynomials_are_merged():
    ctx = Context(["x", "y"])
    system = Poly_system(ctx, [ctx.parse("2*x + 4*y"), ctx.parse("-x - 2*y"), ctx.parse("0")])
    assert system.l == 1


def test_main_variable_and_degrees():
    ctx = Context(["x", "y", "z"])
    p = Poly_with_degrees(ctx, ctx.parse("y^3 + z"))
    assert p.x_m == "y" and p.n == 3 and p.degrees == (0, 3, 1)


# ---- системы уравнений -------------------------------------------------------
def test_rational_roots():
    ctx = Context(["x", "y"])
    got = solve_system(ctx, ["x^2 - 3*x + 2", "y - x"], verbose=False)
    assert floats(got) == [(1.0, 1.0), (2.0, 2.0)]


def test_irrational_roots_circle_line():
    ctx = Context(["x", "y"])
    got = solve_system(ctx, ["x^2 + y^2 - 1", "x - y"], verbose=False)
    r = 2 ** -0.5
    assert floats(got) == [(round(-r, 9), round(-r, 9)), (round(r, 9), round(r, 9))]


def test_two_conics_four_points():
    ctx = Context(["x", "y"])
    polys = ["x^2 + y^2 - 4", "x*y - 1"]
    got = solve_system(ctx, polys, verbose=False)
    assert len(got) == 4
    residual_ok(ctx, polys, got, tol=1e-30)


def test_root_not_found_by_sympy_solve_poly_system():
    # x^3 = y, y^2 = x + 2  <=>  x^6 - x - 2 = 0: два вещественных корня, x = -1 и x = 1.21486...
    ctx = Context(["x", "y"])
    polys = ["x^3 - y", "y^2 - x - 2"]
    got = solve_system(ctx, polys, verbose=False)
    assert floats(got)[0] == (-1.0, -1.0)
    assert len(got) == 2 and abs(floats(got)[1][0] - 1.2148623) < 1e-6
    residual_ok(ctx, polys, got, tol=1e-30)


def test_multiple_root_isolated_point():
    ctx = Context(["x", "y"])
    got = solve_system(ctx, ["x^2 + y^2", "x - y"], verbose=False)
    assert floats(got) == [(0.0, 0.0)]


def test_tangency_after_irrational_substitution():
    # после подстановки x = ±sqrt(2) многочлен по y имеет двойной корень (численный режим)
    ctx = Context(["x", "y"])
    polys = ["x^2 - 2", "y^2 - 2*x*y + x^2"]
    got = solve_system(ctx, polys, verbose=False)
    r = 2 ** 0.5
    assert floats(got) == [(round(-r, 9), round(-r, 9)), (round(r, 9), round(r, 9))]


def test_inconsistent_system_has_no_roots():
    ctx = Context(["x", "y"])
    assert len(solve_system(ctx, ["x^2 + 1", "y - x"], verbose=False)) == 0
    assert len(solve_system(ctx, ["x*y - 1", "x*y - 2"], verbose=False)) == 0


def test_three_variables_triple_root():
    # x+y+z = 3, x^2+y^2+z^2 = 3, xyz = 1  <=>  (t - 1)^3 = 0: единственное решение (1, 1, 1)
    ctx = Context(["x", "y", "z"])
    got = solve_system(ctx, ["x + y + z - 3", "x^2 + y^2 + z^2 - 3", "x*y*z - 1"], verbose=False)
    assert floats(got) == [(1.0, 1.0, 1.0)]


def test_three_variables_no_real_solutions():
    # x+y+z = 1, x^2+y^2+z^2 = 1, xyz = 1/27: кубика t^3 - t^2 - 1/27 имеет один вещественный корень
    ctx = Context(["x", "y", "z"])
    got = solve_system(ctx, ["x^2 + y^2 + z^2 - 1", "x + y + z - 1", "x*y*z - 1/27"], verbose=False)
    assert len(got) == 0


# ---- глобальный минимум --------------------------------------------------------
def minimum(names, poly, lo, hi, **kw):
    ctx = Context(names)
    bounds = {i: (Context.parse_number(str(lo)), Context.parse_number(str(hi))) for i in range(ctx.m)}
    result, _ = global_minimum(ctx, ctx.parse(poly), bounds, verbose=False)
    return result


def test_min_paraboloid():
    r = minimum(["x", "y"], "x^2 + y^2", -2, 2)
    assert r.value == 0.0 and floats(r.points) == [(0.0, 0.0)]


def test_min_rosenbrock():
    r = minimum(["x", "y"], "(1 - x)^2 + 100*(y - x^2)^2", -2, 2)
    assert abs(r.value) < 1e-12 and floats(r.points) == [(1.0, 1.0)]


def test_min_on_boundary_without_interior_critical_points():
    r = minimum(["x", "y"], "x + 2*y", 0, 1)
    assert r.value == 0.0 and floats(r.points) == [(0.0, 0.0)]


def test_min_at_face_point():
    # y^2 - x: внутри критических точек нет, минимум на грани x = 1 при y = 0
    r = minimum(["x", "y"], "y^2 - x", -1, 1)
    assert r.value == -1.0 and floats(r.points) == [(1.0, 0.0)]


def test_min_with_free_variable_does_not_crash():
    # раньше программа падала на exit(): многочлен не зависит от y
    r = minimum(["x", "y"], "x^2", -1, 1)
    assert r.value == 0.0 and len(r.points) >= 1


def test_min_three_variables_two_global_minima():
    r = minimum(["x", "y", "z"],
                "(-2*x - 2*y + 3*z^2 + 1)^2 + 3*(-2*x - 3*y^2 + 3*z + 2)^2 + 2*(-x^2 + 3*y - z - 1)^2", -4, 4)
    assert abs(r.value) < 1e-12
    pts = floats(r.points)
    assert (1.0, 1.0, 1.0) in pts and len(pts) == 2


def test_min_five_variables_degree_seven_finds_all_faces():
    r = minimum(["a", "b", "c", "d", "e"], "a^7 + b^7 + c^7 + d^7 + e^7 + a*b*c*d*e", -2, 2)
    assert r.value == -672.0 and floats(r.points) == [(-2.0,) * 5]
    # На грани a = -2, e = 2 система 7b^6 = 4cd, 7c^6 = 4bd, 7d^6 = 4bc имеет решение b = c = d = -t,
    # t = (4/7)^(1/4) = 0.8693...: эта критическая точка раньше терялась (сравнение float-корней).
    t = (4 / 7) ** 0.25
    found = {tuple(round(f, 6) for _, f in full.float_items()) for full, _, _ in r.candidates}
    assert (-2.0, round(-t, 6), round(-t, 6), round(-t, 6), 2.0) in found


# ---- консольный интерфейс -----------------------------------------------------
def run_cli(stdin, *args):
    return subprocess.run([sys.executable, str(ROOT / "main.py"), *args], input=stdin, capture_output=True,
                          text=True, cwd=ROOT, timeout=300)


def test_cli_solve_quiet():
    out = run_cli("x y\n1\n2\nx^2 - 1\ny - x\n", "-q")
    assert out.returncode == 0, out.stderr
    assert "{x: -1.0, y: -1.0}" in out.stdout and "{x: 1.0, y: 1.0}" in out.stdout
    assert "ALPHA" not in out.stdout


def test_cli_min_matches_original_format():
    out = run_cli("x y\n2\nx^2 + y^2\n-2\n2\n-2\n2\n")
    assert out.returncode == 0, out.stderr
    assert "Global minimum = 0.0 and is achieved at:" in out.stdout
    assert "Point #1: {x: 0.0, y: 0.0}" in out.stdout


def example_io(name):
    """Ввод и ожидаемый итог из файла examples/<name>.txt (формат: вывод консольной сессии)."""
    text = (ROOT / "examples" / f"{name}.txt").read_text(encoding="utf-8")
    prompts = (r"Enter variables separated by spaces|Enter 1 or 2|Enter the number of polynomials|"
               r"Enter polynomial(?: #\d+)?|Lower bound for \w+|Upper bound for \w+")
    inputs = [m.group(1) for m in re.finditer(rf"^(?:{prompts}): (.*)$", text, flags=re.M)]
    tail = text[text.index("Global minimum"):].strip()
    return "\n".join(inputs) + "\n", tail


@pytest.mark.parametrize("name", ["00_example", "01_example", "02_example", "04_example"])
def test_examples_reproduce(name):
    stdin, expected_tail = example_io(name)
    out = run_cli(stdin, "-q")
    assert out.returncode == 0, out.stderr
    got_tail = out.stdout[out.stdout.index("Global minimum"):].strip()
    assert got_tail == expected_tail
