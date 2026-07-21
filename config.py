from sympy import symbols

syms = symbols(input("Введите переменные через пробел: "))
try:
    variables = list(syms)
except:
    variables = [syms]

current_values = dict.fromkeys(variables, None)
round_parameter = 5  # кол-во знаков при округлении для печаи в терминал
close_zero = 0.0001