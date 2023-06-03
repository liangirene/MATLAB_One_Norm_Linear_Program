# One-Norm-LP
The OneNormLP function solves a L^1 norm regression problem
- obj denotes the objective value 
- x denotes the x s.t |Ax-b| achieve the minimum value
- status contains several messages:
  - 1) degenerate
  - 2) optimal reached --> obj converges
The status of the problem is determined based on two conditions, as there are only four possible situations resulting from the duality theorem. The situation where the primal is feasible but lacks a lower bound is excluded due to the presence of the one_norm. Furthermore, since no restrictions are imposed on x, it is considered feasible. As a result, the only two remaining conditions to consider are whether the problem is degenerate
