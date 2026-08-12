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
