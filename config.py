from sympy import symbols

syms = symbols(input("Enter variables separated by spaces: "))
try:
    variables = list(syms)
except:
    variables = [syms]

current_values = dict.fromkeys(variables, None)
# Number of decimal places for rounding when printing to terminal
round_parameter = 5
close_zero = 0.0001
