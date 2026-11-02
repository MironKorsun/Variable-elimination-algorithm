# Variable Elimination Algorithm - VES

![Static Badge](https://img.shields.io/badge/MironKorsun-VES-green)
![GitHub top language](https://img.shields.io/github/languages/top/MironKorsun/Variable-elimination-algorithm)
![GitHub](https://img.shields.io/github/license/MironKorsun/Variable-elimination-algorithm)
![GitHub Repo stars](https://img.shields.io/github/stars/MironKorsun/Variable-elimination-algorithm)
![GitHub issues](https://img.shields.io/github/issues/MironKorsun/Variable-elimination-algorithm)

## About
This project implements the **Variable Elimination Algorithm (VEA)** for:
- Solving systems of polynomial equations
- Finding the global minimum of a polynomial on a coordinate parallelepiped

Based on the paper by **V.N. Nefedov** ([arXiv:2607.13658](https://arxiv.org/abs/2607.13658)).

## Installation (Windows)
Follow these steps to set up the project on Windows:

1. **Clone the repository:**
```git clone https://github.com/MironKorsun/Variable-elimination-algorithm.git```

2. **Navigate to the project directory:**
```cd Variable-elimination-algorithm```

3. **Create a virtual environment:**
```python -m venv venv```

4. **Activate the virtual environment:**
```venv\Scripts\activate.bat```

5. **Install dependencies:**
```pip install -r requirements.txt```

6. **Run the main script:**
```python main.py```

   Add `-q` (`python main.py -q`) to print only the answer instead of the whole trace of the algorithm.

7. **Run the tests (optional):**
```pip install -r requirements-dev.txt```
```pytest```

## Documentation
Full user documentation is available [here](./docs/user-guide.md).
For examples, check out the [examples](./examples) directory.

## Support and contacts
If you encounter any difficulties or have questions regarding the use of the code, please create a 
[discussion](https://github.com/MironKorsun/Variable-elimination-algorithm/issues/new/choose) in this repository or send an email to <Miron.Korsun@hotmail.com>.

## Requirements
- **Python** 3.9 or higher (the minimum supported by `sympy` 1.14; tested on Python 3.13)
- Dependencies are listed in `requirements.txt` (only `sympy`)

## How it works / what is exact
Polynomials are stored in `QQ[x1..xm]` with exact rational coefficients, so the symbolic steps (the
reductions R, the GCD of leading coefficients from Remark 1, derivatives, the α/β/γ/δ transformations)
involve no rounding. Real roots are isolated rigorously (Sturm sequences) and multiple roots are removed by
`p / gcd(p, p')`. Only after an irrational root is substituted into the system does the arithmetic become
approximate (fixed-point rationals with a step of 2^-300); see [the user guide](./docs/user-guide.md).
All branches of the algorithm are explored — the paper guarantees that the union of the alternative
problems contains every solution, so there is no early exit after the first roots.

The method is practical for low-dimensional problems (the paper mentions up to 4–5 variables) and polynomials
of moderate degree.

> **Note:** If you discover that it works on a lower version or fails on a specific version, please let me know via [support](https://github.com/MironKorsun/Variable-elimination-algorithm#support-and-contacts).
