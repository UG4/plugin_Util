# Tests

This directory contains tests for the C++ implementation of the solver utilities.

## Successfully tested

- Laplace problem with BiCGStab and geometric multigrid.
- Results compared with the original Lua implementation.
- MGStats disabled (`mgStats = nil`).
- Standard MGStats configuration (`mgStats = "standard"`).
- Custom MGStats configuration.
- MGStats parameters are correctly resolved and passed to UG4.
- Standard and custom MGStats tests reproduce the convergence behavior of the original Lua implementation.
- Generated `.vec` and `.vtu` results matched the corresponding original Lua results in the tested Laplace cases.

- Poisson problem using BiCGStab with geometric multigrid.
- Cooler problem using BiCGStab with geometric multigrid.
- Stationary Henry problem using Newton, BiCGStab, geometric multigrid, ILU smoothing, convergence checks, and standard line search.
- Poisson, Cooler, and Henry tests reproduce the solver behavior of the corresponding original Lua examples.
- Generated solution files for the tested Poisson, Cooler, and Henry cases matched the corresponding results of the original Lua implementation.
- The Henry test exercises the nested solver configuration
  `Newton -> BiCGStab -> GMG -> ILU`, including convergence checks and line search.
