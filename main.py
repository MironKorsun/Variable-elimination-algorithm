# from AIP import *
import config
from implemented_structures import Poly_with_degrees, Poly_system, Solutions_set
import VES

VES.all_roots = Solutions_set([])
VES.initial_system = Poly_system([])

print("\nХотите: 1) Решить систему или 2) Найти глобальный минимум?")
match input("Введите 1 или 2: "):
    case '1':
        for i in range(int(input("\nВведите количество многочленов: "))):
            VES.initial_system.insert(input(f"Введите многочлен номер {i+1}: "))

        print("\nИсходная задача")
        VES.ves(VES.initial_system)
        print("\nОтвет:")
        print(VES.all_roots)

    case '2':
        p = Poly_with_degrees(input("\nВведите многочлен: "))
        for v in config.variables:
            VES.initial_system.insert(p.diff(v))

        print("\nВведите границы:")
        borders = dict()
        for var in config.variables:
            borders[var] = [
                float(input(f"Левая граница по {var}: ")),
                float(input(f"Правая граница по {var}: "))
            ]
        print("\nБез учета границ")
        VES.ves(VES.initial_system)

        for x_i in borders:
            for b in borders[x_i]:
                p_subs = p.subs(x_i, b)
                VES.initial_system = Poly_system([])
                for v in config.variables:
                    VES.initial_system.insert(p_subs.diff(v))
                print(f"\nГрань {x_i} = {b}")
                config.current_values[x_i] = b
                VES.ves(VES.initial_system)

        print("\nЭкстремумы на рассматриваемой области:")
        global_min = None
        min_value = None
        for r in VES.all_roots:
            if all(b[0] <= r.solution[var] and r.solution[var] <= b[1] for var, b in borders.items()):
                value = float(p.subs(r.solution).poly.as_expr())
                print("Точка:", r, "Значение:", round(value, config.round_parameter))
                if min_value == None or value < min_value:
                    min_value = value
                    global_min = [r]
                elif value == min_value:
                    global_min.append(r)
        print(f"\nГлобальный минимум = {min_value} и достигается в:")
        for i, point in enumerate(global_min):
            print(f"Точке {i+1}: {point}")
        

       
