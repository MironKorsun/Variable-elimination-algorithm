from itertools import combinations, product
import config
from implemented_structures import Poly_with_degrees, Poly_system, Solutions_set
import VES

VES.all_roots = Solutions_set([])
VES.initial_system = Poly_system([])

print("\nDo you want to: 1) Solve the system or 2) Find the global minimum?")
match input("Enter 1 or 2: "):
    case '1':
        for i in range(int(input("\nEnter the number of polynomials: "))):
            VES.initial_system.insert(input(f"Enter polynomial #{i+1}: "))

        print("\nOriginal problem:")
        VES.ves(VES.initial_system)
        print("\nAnswer:")
        print(VES.all_roots)

    case '2':
        p = Poly_with_degrees(input("\nEnter polynomial: "))
        for v in config.variables:
            VES.initial_system.insert(p.diff(v))

        print("\nEnter bounds:")
        borders = dict()
        for var in config.variables:
            borders[var] = [
                float(input(f"Lower bound for {var}: ")),
                float(input(f"Upper bound for {var}: "))
            ]
        print("\nWithout considering bounds:")
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
                    
                    print(f"\nBound {substitution}:")
                    VES.ves(VES.initial_system)

                    for v in substitution:
                        config.current_values[v] = None  

        print("\nCritical points in the domain:")
        global_min = None
        min_value = None
        for i, r in enumerate(VES.all_roots):
            try:
                if all(b[0] <= r.solution[var] and r.solution[var] <= b[1] for var, b in borders.items()):
                    value = float(p.subs(r.solution).poly.as_expr())
                    print(f"Point #{i+1}: {r} Value: {round(value, config.round_parameter)}")
                    if min_value == None or value < min_value:
                        min_value = value
                        global_min = [r]
                    elif value == min_value:
                        global_min.append(r)
            except:
                print(r)
                exit()
        print(f"\nGlobal minimum = {min_value} and is achieved at:")
        for i, point in enumerate(global_min):
            print(f"Point #{i+1}: {point}")
