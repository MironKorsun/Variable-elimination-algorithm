from itertools import combinations, product
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
                float(input(f"Нижняя граница по {var}: ")),
                float(input(f"Верхняя граница по {var}: "))
            ]
        print("\nБез учета границ")
        VES.ves(VES.initial_system)

        for i in range(len(config.variables)):
            
            for keys in combinations(borders.keys(), i+1):
                values = [borders[key] for key in keys]
                for combo in product(*values):
                    substitution = dict(zip(keys, combo))
                    for var, val in substitution.items():
                        config.current_values[var] = val

                    p_subs = p.subs(substitution)
                    VES.initial_system = Poly_system([])
                    for var in config.variables:
                        VES.initial_system.insert(p_subs.diff(var))
                    
                    print(f"\nГраница {substitution}:")
                    VES.ves(VES.initial_system)

                    for v in substitution:
                        config.current_values[v] = None  

        print("\nКритические точки на рассматриваемой области:")
        global_min = None
        min_value = None
        for i, r in enumerate(VES.all_roots):
            try:
                if all(b[0] <= r.solution[var] and r.solution[var] <= b[1] for var, b in borders.items()):
                    value = float(p.subs(r.solution).poly.as_expr())
                    print(f"Точка {i+1}: {r} Значение: {round(value, config.round_parameter)}")
                    if min_value == None or value < min_value:
                        min_value = value
                        global_min = [r]
                    elif value == min_value:
                        global_min.append(r)
            except:
                print(r)
                exit()
        print(f"\nГлобальный минимум = {min_value} и достигается в:")
        for i, point in enumerate(global_min):
            print(f"Точке {i+1}: {point}")