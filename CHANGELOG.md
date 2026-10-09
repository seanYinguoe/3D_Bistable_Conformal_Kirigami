# Changes

## Public release — 9 October 2026

- Rechecked all 10 fast MATLAB tests and the dependency closure of the three main programs.
- Removed COMSOL model files and sidecars from branch histories, preserving a local backup of the previous repository. The maintained source and input data are unchanged.
- Updated provenance references to the rewritten commits and strengthened model-file exclusions.

## Focused Project 2 cleanup — 8 October 2026

- Replaced legacy drivers and introductory demos with three entry points: `main_surface`, `main_unit_energy` and `main_unit_library`.
- Kept the 41 source functions reachable from those programs and organised them by surface mapping, geometry, design, mechanics, library processing, plotting and export.
- Removed unused exploratory programs, experiment comparisons, duplicate/raw mesh variants, saved COMSOL models and obsolete documentation from the current tree. Git history preserves earlier source versions; model files are excluded from the public branch histories.
- Retained three target/UV meshes, the selected historical unit library, citation metadata and the OBJ reader's licence.
- Moved run choices into `config/`; isolated generated runs under ignored `results/` folders.
- Added saved reparameterisation, assignment and deployed-unit diagnostics. A nearest library match is reported separately from meeting the strain tolerance.
- Added convergence checks to the anisotropic energy path. Testing exposed an infeasible second deployment step in the unscaled solver (exit flag 0, constraint residual 0.0472). Scaling the objective and its gradient by a positive stiffness constant resolved the checked cases; reported energy still uses the original expression. This is a numerical-conditioning change, not a new constitutive law, and older numerical trajectories may differ.
- Prevented invalid units from silently disappearing from an energy sum.
- Added compact usage instructions, input provenance and a validation record.

### Provenance

The original source baseline is Git commit `d978245`. The first organisational pass is preserved at `1e6bb81`, immediately before this cleanup. The full surface workflow comes from `Xiaoyuan_bistable_kirigami_unit/main_3D.m`; the library workflow comes from `Inverse_design/ternary_study.m` within that source tree. The unused files and old source manifest can be inspected in those commits.

The original research directory, including uncommitted quarter-dome fabrication changes, was not edited. This maintained version is not represented as the exact code snapshot used for every manuscript figure. The manuscript remains **under review**.
