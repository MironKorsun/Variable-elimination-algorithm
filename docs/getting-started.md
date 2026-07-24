## Preface
This project addresses the problem of finding the **global minimum of a polynomial** on a coordinate parallelepiped, as well as solving systems of an arbitrary number of polynomials using the **Variable Elimination Algorithm (VEA)**.

### Core Idea
The algorithm reduces the degree of a selected variable through a sequence of **α-, β-, γ-, δ-transformations**, which depend on specific characteristics of the polynomial system. This process continues until the variable is completely eliminated.

### Theoretical Foundation
According to the paper by **V.N. Nefedov**, these transformations preserve all roots of the system, ensuring that no solutions are lost during the elimination process.

> **Reference:** Nefedov V.N. "The problem of finding the global minimum of a polynomial on a coordinate parallelepiped." arXiv:2607.13658 (2026). [https://arxiv.org/abs/2607.13658](https://arxiv.org/abs/2607.13658)

### Practical Considerations
As noted in the paper, the method is most effective for **low-dimensional problems** and polynomials of **not very high degree**. While theoretically applicable to any polynomial optimization problem, the implementation may lead to a large number of alternative problems concerning the real roots of polynomials, which can be of very high degree.
