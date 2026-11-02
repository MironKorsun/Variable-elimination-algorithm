"""Консольный интерфейс: решение системы полиномиальных уравнений или глобальный минимум
многочлена на координатном параллелепипеде (алгоритм исключения переменных, VES).

Запуск:  python main.py [-q]      (-q / --quiet: не печатать ход алгоритма)
"""
from __future__ import annotations

import argparse

import config
from VES import global_minimum, solve_system
from implemented_structures import Context


def main(argv=None):
    parser = argparse.ArgumentParser(description="Variable Elimination Algorithm (VES)")
    parser.add_argument("-q", "--quiet", action="store_true", help="не печатать ход алгоритма, только ответ")
    args = parser.parse_args(argv)
    verbose = not args.quiet

    names = input("Enter variables separated by spaces: ").replace(",", " ").split()
    ctx = Context(names)

    print("\nDo you want to: 1) Solve the system or 2) Find the global minimum?")
    mode = input("Enter 1 or 2: ").strip()

    if mode == "1":
        count = int(input("\nEnter the number of polynomials: "))
        polys = [ctx.parse(input(f"Enter polynomial #{i + 1}: ")) for i in range(count)]
        if verbose:
            print("\nOriginal problem:")
        roots = solve_system(ctx, polys, verbose)
        print("\nAnswer:")
        print(roots)

    elif mode == "2":
        p = ctx.parse(input("\nEnter polynomial: "))
        print("\nEnter bounds:")
        bounds = {}
        for i, name in enumerate(ctx.names):
            lo = ctx.parse_number(input(f"Lower bound for {name}: "))
            hi = ctx.parse_number(input(f"Upper bound for {name}: "))
            if lo > hi:
                raise SystemExit(f"Lower bound for {name} is greater than the upper bound")
            bounds[i] = (lo, hi)

        result, _ = global_minimum(ctx, p, bounds, verbose)
        if result.value is None:
            print("\nNo critical points in the domain were found.")
            return
        print(f"\nGlobal minimum = {round(result.value, config.ROUND_PARAMETER)} and is achieved at:")
        for i, point in enumerate(result.points):
            print(f"Point #{i + 1}: {point}")

    else:
        raise SystemExit("Enter 1 or 2")


if __name__ == "__main__":
    main()
