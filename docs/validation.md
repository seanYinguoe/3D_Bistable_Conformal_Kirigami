# Validation record

Checked on 8 October 2026 with MATLAB R2025a Update 1 on macOS. The manuscript is a preprint **under review**.

## Automated checks

4/4 tests passed:

- Triangle area under uniform scaling.
- Edge-stretch calculation under uniform expansion.
- Matching 3D/UV connectivity, finite coordinates and positive areas for all three selected target meshes.
- Required unit-library columns and the triangle-angle sum.

The test count is checked explicitly, preventing an empty suite from being reported as successful.

## Conformal-mesh example

`demo_conformal_mesh('quarter_dome')` read 997 vertices and 1,812 faces. The area-derived scale ranged approximately from 0.999 to 1.332. It saved coordinates, connectivity, scales and a figure. This checks the supplied map; it does not rerun BFF.

## Energy example

`demo_unit_energy` ran 81 strain steps for each of two tilting angles, with eight chain segments per ligament and the settings saved by the example.

| Tilting angle | Energy solves | All exit flags positive | Maximum constraint residual |
|---|---:|---|---:|
| 0 | 81 | Yes | 1.78e-15 |
| pi/40 | 81 | Yes | 2.22e-15 |

All returned energies were finite and nonnegative. The plotted curves and conformal-mesh illustration were visually inspected. These checks establish numerical completion of the examples, not global optimality or experimental calibration.

### Geometry and energy are checked separately

The original isotropic routine reconstructs a final rigid-link geometry after calculating elastic energies. At the end of the example's strain range, that additional geometry solve did not satisfy its equations. The new energy-only option skips that reconstruction; it does not change the energy equations or their solver settings. The example records energy-solver diagnostics and fails if a solve is unsuccessful or infeasible. The optional geometry solver now exposes its residual and warns on failure.

## Not rerun

- Full-surface reparameterisation, every unit assignment and assembled-structure validation.
- The complete library-generation sweeps.
- COMSOL finite-element studies and full experimental comparisons.
- Exact reproduction of every preprint figure.

The retained full drivers and models are research workflows, separate from the verified small examples.
