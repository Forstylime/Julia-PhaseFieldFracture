# PhaseFieldFracture

A Julia research codebase for phase-field fracture simulations, focused on a
computationally efficient monolithic arc-length method based on the
Efendiev-Mielke scheme.

The project skeleton mirrors the paper workflow:

- physics and constitutive models in `src/physics`
- finite-element assembly primitives in `src/fem`
- monolithic, staggered, and arc-length solvers in `src/solvers`
- simulation entry points in `scripts`
- verification tests in `test`
- post-processing notebooks in `notebooks`

Simulation outputs are written under `data/`, which is intentionally ignored by
Git.
