# Running and extending the code

## Requirements

| Task | Requirements |
|---|---|
| Supplied conformal-mesh demo and data checks | MATLAB; checked with R2025a |
| Unit geometry / energy minimisation | Optimization Toolbox (`fsolve`, `fmincon`) |
| Optional library parameter sweeps | Parallel Computing Toolbox where `parfor` is used |
| Some exploratory filtering scripts | Statistics and Machine Learning Toolbox |
| Creating a new target's UV map | Boundary First Flattening or an equivalent, separately installed mapping tool |
| Optional FE comparison | COMSOL Multiphysics with suitable structural-mechanics capability; LiveLink for MATLAB for automation |
| Legacy image-based experiment comparison | Image Processing Toolbox |

The minimum supported software versions have not been established. No Octave compatibility or untested COMSOL compatibility is claimed.

## Example outputs

- `demo_conformal_mesh`: paired UV/3D coordinates, face connectivity, area-derived scale, MATLAB version and PNG.
- `demo_unit_energy`: geometry settings, strains, raw energies, solver exit flags and constraint residuals, MATLAB version, CSV and PNG. CSV columns are strain, energy for beta = 0, energy for beta = pi/40.

Outputs are written to unique folders under `results/`. Plot energies are normalized separately for readability; saved energies are unnormalized. The inherited unit model contains `Emod = 4.3e11` and `b = 1.0`; these constants have not been silently replaced with new material data. The example illustrates the model's response and is not an absolute-energy calibration against an experiment. Check consistent units and material parameters before scientific use.

## Full surface workflow

After `setup_project`, the original driver can be inspected or run with `run('legacy/main_3D.m')`. Its default is `double_dome`; the included alternatives are `quarter_dome` and `hemisphere`. It loads the supplied mesh and lookup table using repository-relative paths. It creates interactive figures and may take substantially longer than the introduction examples.

The driver uses grid edge length 16 and the committed ligament-thickness settings. Recent uncommitted fabrication trials in the local research folder were not imported. The numerical library comes from `ternary_refine_20260325_161516/anisotropy_filter_refine.mat`, exactly the file selected by that driver, rather than the many similarly named intermediate tables.

Preserve the third output of `assign_opt_beta` when adapting this workflow. Its fields include `flag`, `eps_pred`, `eta_pred`, `beta_error` and local sample counts. Boundary units are handled specially: the current implementation sets their beta to `pi/12`; some old comments describing beta = 0 were stale. Unit-library interpolation is not an independent simulation of the assembled structure.

The `legacy/main_2D.m` driver is an additional planar workflow. COMSOL models live in `models/comsol`; they are optional and were not part of the quick-start validation. Original experiment-analysis scripts and their small input files are retained together in `legacy/experiments`.

## Source organisation and limits

`src/conformal` contains the original `Conformal_mapping` functions, `src/units` the `Generate_hexagon` code, and `src/inverse_design` the `Inverse_design` code. Other original drivers remain under `legacy`. `source-manifest.csv` maps every baseline file to the new tree (or preserved Git history), and `source-commit.txt` records that baseline.

This cleanup retains numerical equations. Changes concern entry points, input/output locations, documentation, UTF-8 text where needed, and file-read errors. Generated trial exports, duplicate table snapshots, editor files and unrelated example meshes were omitted from the current tree; their historical versions remain in Git. No experimental research folder on the original machine was modified.

The code is not a complete immutable archive of all figures and experiments in the paper. Keep a tagged code version and the exact input table with future publication results.

## Third-party reader

`third_party/Objread` retains the original BSD-3-Clause licence. The cleanup changes path joining to `fullfile`, reports unreadable files clearly and closes files on errors. The numerical mesh contents are unchanged. BFF itself is not bundled; the included meshes contain precomputed mapping results.
