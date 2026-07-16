# from AIP import *
import config
from implemented_structures import Poly_with_degrees, Poly_system, Solutions_set
import AIP

AIP.all_roots = Solutions_set([])
AIP.initial_system = Poly_system([])

print("\nХотите: 1) Решить систему или 2) Найти глобальный минимум?")
match input("Введите 1 или 2: "):
    case '1':
        for i in range(int(input("\nВведите количество многочленов: "))):
            AIP.initial_system.insert(input(f"Введите многочлен номер {i+1}: "))

        print("\nИсходная задача")
        AIP.aip(AIP.initial_system)
        print("\nОтвет:")
        print(AIP.all_roots)

    case '2':
        p = Poly_with_degrees(input("\nВведите многочлен: "))
        for v in config.variables:
            AIP.initial_system.insert(p.diff(v))

        print("\nВведите границы:")
        borders = [[0, 0] for _ in config.variables]
        for i in range(len(config.variables)):
            borders[i][0] = float(input(f"Левая граница по {config.variables[i]}: "))
            borders[i][1] = float(input(f"Правая граница по {config.variables[i]}: "))

        print("\nОбласть внутри")
        AIP.aip(AIP.initial_system)

        for i in range(len(config.variables)):
            for b in borders[i]:
                p_subs = p.subs(config.variables[i], b)
                AIP.initial_system = Poly_system([])
                for v in config.variables:
                    AIP.initial_system.insert(p_subs.diff(v))
                print(f"\nГрань {config.variables[i]} = {b}")
                config.current_values[config.variables[i]] = b
                AIP.aip(AIP.initial_system)

        print("\nЭкстремумы:")
        print(AIP.all_roots)

        print("\nОтвет:")
        global_min = None
        min_value = None
        for r in AIP.all_roots:
            value = float(p.subs(r.solution).poly.as_expr())
            if min_value == None or value < min_value:
                min_value = value
                global_min = r
        print(f"Значение: {min_value}")
        print(f"Точка: {global_min}")

       
