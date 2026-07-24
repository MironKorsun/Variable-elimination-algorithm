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

It is important to note that roots are found in only one branch of the VEA. Once roots are found, there is no need to continue searching in other branches — the program execution is interrupted, and the algorithm exits the recursive tree.

### Data Structures
To optimize the search for characteristics of the polynomial system, such as `p1` from Nefedov's paper or the number `k` for selecting the appropriate transformation, specialized classes were introduced:

- **`Poly_with_degrees`** — a polynomial with the degree of each variable stored as a tuple. This allows for rapid identification of the main variable and its degree for elimination.
- **`Poly_system`** — a storage system based on a sorted set from the `sortedcontainers` library. Sorting is performed by comparing the tuples of variable degrees of the polynomials in the system, providing fast access to the first polynomial (`p1`) and other system characteristics.

Additionally, a case was handled where a polynomial consists of a single constant — in such a case, the system has no solutions, and the branch of the recursive VEA tree is terminated.

### Practical Considerations
As noted in the paper, the method is most effective for **low-dimensional problems** and polynomials of **not very high degree**. While theoretically applicable to any polynomial optimization problem, the implementation may lead to a large number of alternative problems concerning the real roots of polynomials, which can be of very high degree.
