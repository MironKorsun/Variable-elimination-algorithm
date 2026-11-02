## Preface
This project addresses the problem of finding the **global minimum of a polynomial** on a coordinate parallelepiped, as well as solving systems of an arbitrary number of polynomials using the **Variable Elimination Algorithm (VEA)**.

### Core Idea
The algorithm reduces the degree of a selected variable through a sequence of **α-, β-, γ-, δ-transformations**, which depend on specific characteristics of the polynomial system. This process continues until the variable is completely eliminated.

### Theoretical Foundation
According to the paper by **V.N. Nefedov**, these transformations preserve all roots of the system, ensuring that no solutions are lost during the elimination process.

> **Reference:** Nefedov V.N. "The problem of finding the global minimum of a polynomial on a coordinate parallelepiped." arXiv:2607.13658 (2026). [https://arxiv.org/abs/2607.13658](https://arxiv.org/abs/2607.13658)

### User Workflow
The user independently defines which variables are used in their problem. It is important that no extra variables are introduced, and that no used variable is forgotten. Variables can consist of a single character or multiple characters.

After entering the variables, the user selects the mode of operation: find the global minimum of a single polynomial and the points at which it is attained, or solve a system of polynomials.

In the case of finding the global minimum, the problem reduces to solving a system of **m** polynomials — the partial derivatives with respect to each of the **m** variables. The user also enters the bounds of the search domain, which are combinatorially substituted into the original polynomial, after which it is differentiated.

### Algorithm Overview
The essence of solving the system lies in the recursive application of the VEA. As long as the system exists (i.e., it contains more than 0 polynomials), the algorithm either searches for roots in univariate polynomials, substitutes them into the system, and restarts the VEA, or applies the VEA to each of the **α-, β-, γ-, δ-transformations**.

All branches are always explored: by the paper, the union of the solutions of the alternative problems contains every solution of the original one, so stopping after the first roots would lose solutions. A branch ends when it has no real roots ("No solutions"), when a nonzero constant appears, or when a root is substituted and the remaining system is solved recursively. Candidate roots are finally checked against the original system, and extraneous roots are discarded.

Roots of the univariate polynomials are handled as follows:
- While the system is exact (no irrational root has been substituted yet), the common roots of all univariate polynomials are the roots of their GCD, with multiple roots removed (`p / gcd(p, p')`, Remark 2 of the paper). Rational roots are found exactly, irrational ones to within 10^-70.
- After an irrational root is substituted, the system is stored with fixed-point rational coefficients (step 2^-300), and noise coefficients are zeroed relative to the size of the terms that produced them. Common roots are then taken from the polynomial of the lowest degree and from its derivatives (a multiple root may detach from the real axis under perturbation, but a root of the derivative stays close) and are accepted only if every univariate polynomial vanishes there up to a tiny residual.

### Data Structures
To optimize the search for characteristics of the polynomial system, such as `p1` from Nefedov's paper or the number `k` for selecting the appropriate transformation, specialized classes were introduced:

- **`Context`** — the polynomial ring `QQ[x1..xm]` (from `sympy.polys.rings`) for the variables entered by the user. Variables are entered once, in `main.py`; importing `config.py` no longer asks for input. The first variable is eliminated first.
- **`Poly_with_degrees`** — a polynomial with the degree of each variable stored as a tuple. This allows for rapid identification of the main variable and its degree for elimination. Coefficients are exact rationals: input such as `0.5*x` is read as `1/2`, and no floating-point arithmetic or `simplify` is used. In exact mode every polynomial is stored as a primitive integer polynomial with a positive leading coefficient, so proportional polynomials are the same object.
- **`Poly_system`** — the system of equations, kept sorted by (main variable, its degree, degrees of the leading coefficient, ...), which gives `p1` as in Remarks 2–3 of the paper, together with `l`, `k`, `n1`. No third-party sorted container is needed any more.
- **`Solution`, `Solutions_set`** — solutions as exact rationals (approximate for irrational roots); duplicates are detected with a relative tolerance of 1e-9.

Additionally, a case was handled where a polynomial consists of a single constant — in such a case, the system has no solutions, and the branch of the recursive VEA tree is terminated.

### Practical Considerations
As noted in the paper, the method is most effective for **low-dimensional problems** and polynomials of **not very high degree**. While theoretically applicable to any polynomial optimization problem, the implementation may lead to a large number of alternative problems concerning the real roots of polynomials, which can be of very high degree.
