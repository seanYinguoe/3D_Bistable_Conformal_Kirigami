# Working with Project 2

Run `setup_project` from the repository root. It configures the current MATLAB session without changing the saved path. Use `cfg = surface_config`, `unit_config` or `library_config`, edit the returned struct, and pass it to the corresponding `main_*` function.

## 1. Design a surface pattern

```matlab
setup_project
cfg = surface_config;
result = main_surface(cfg);
```

The surface workflow follows the original `main_3D` driver:

1. Read a target OBJ containing both curved vertices and planar UV coordinates. Align UV indices to mesh connectivity.
2. Fit a regular triangular grid to the planar domain. Barycentrically map its vertices to the surface and attempt the original surface-constrained reparameterisation.
3. Uniformly rescale the target using the admissible library strain range. Convert each unit's three edge stretches into angle/strain descriptors.
4. Interpolate the library to select the closest strain-matching tilt `beta`. Choose a ligament width from the predicted energy-barrier measure.
5. Assemble compact units, optionally solve deployed unit geometries, and export the cut pattern.

The UV map is an **input**: this program does not run Boundary First Flattening. The final target may be uniformly larger than the original mesh; `result.rescaling` records this factor. OBJ coordinates do not specify physical units.

### Main settings

| Setting | Default | Meaning |
|---|---|---|
| `mesh` | `quarter_dome_flat.obj` | Other supplied choices: `double_dome`, `hemisphere` |
| `edge_length` | 16 | Reference triangular-grid edge length, in mesh units |
| `flank_length_ratio`, `flank_width_ratio` | 0.85, 0.05 | Flank geometry relative to grid edge length |
| `ligament_width_ratios` | `[0.015 0.025]` | Narrow and wide ligament choices |
| `library_eta_min` | 0.10 | Filter on the supplied library's bistability measure |
| `thickness_threshold` | 0.40 | Predicted eta at or above this uses the wider ligament |
| `assignment.strainTol` | 0.015 | Maximum strain mismatch counted as within tolerance |
| `build_deployed` | `true` | Also solve each deployed unit's geometry |
| `export_svg`, `fillet_radius` | `true`, 0.30 | Filleted SVG; radius uses mesh units |
| `make_plot` | `true` | Write a summary PNG |
| `evaluate_energy` | `false` | Optional, expensive independent-unit energy sum |

Keep length and material units consistent. SVG fillet radius is not an automatic cutter-kerf correction. The thin-ligament defaults are the committed research settings; thicker uncommitted fabrication trials were left in the original research directory.

### Files and diagnostics

- `surface_result.mat`: configuration, connectivity, flat/target vertices, stretch demands, selected units, geometry and diagnostics.
- `unit_assignments.csv`: rim flag, beta (radians), ligament width, predicted eta, assignment flag, absolute strain mismatch and tolerance pass/fail.
- `cut_pattern.svg`: generated when `export_svg` is enabled.
- `surface_overview.png`: generated when `make_plot` is enabled.

Inspect `result.reparameterization.used_barycentric_fallback` and `stop_reason`. Inspect `geometry_diagnostics.exitflag` and `max_constraint` for the deployed-unit solves.

The inherited assignment flag `eps_matched` means the **nearest available strain match was selected**; it does not guarantee the requested tolerance. Use `assignment_within_tolerance`. Rim units have beta = pi/12, a `rim_undeployed` flag and no predicted bistability value. This boundary prescription is not evidence of a mechanically equilibrated rim.

A generated cut pattern is a candidate for further assessment. The supplied library describes a particular unit geometry; changing ligament width or flank geometry can change its energy landscape. Re-evaluate those units rather than treating interpolated eta as a fresh mechanical simulation.

If `evaluate_energy = true`, `energy.points`, `energy.segments` and `energy.use_parallel` control the HBM sum. It sums independently prescribed unit paths, with no global sheet-equilibrium solve. Invalid unit results are not silently omitted from the total.

## 2. Examine unit energy and bistability

```matlab
setup_project
cfg = unit_config;
cfg.edge_stretches = [1.63 1.63 1.63]; % [lambda12 lambda23 lambda31]
result = main_unit_energy(cfg);
```

The default compares beta = 0 and pi/40 at edge length 15, flank length 12.75, flank width 0.75 and ligament width 0.225. Unequal edge stretches prescribe anisotropic deformation, provided they form a valid triangle.

The Hencky bar-chain model minimises bending and extension energy during incremental deployment, warm-starting each solve from the preceding configuration. The detector looks for a persistent energy maximum followed by a meaningful minimum, with the source code's smoothness and energy-barrier checks. It uses eight segments per ligament and a deployment increment of 1/200. A curve can stop early once classified; `no_transition` means none was detected on the sampled path.

Outputs are `unit_energy.mat`, one CSV per beta, and an optional `unit_energy.png`. The MAT file retains the detector status, extrema indices, sampled energy and solver diagnostics. The figure normalises each curve separately; MAT and CSV retain unnormalised energies.

The inherited model uses `E = 4.3e11` and `b = 1.0`. These have not been calibrated to a new material or mesh-unit convention. The anisotropic solver now divides its objective and gradient by a positive stiffness scale for numerical conditioning, then reports energy using the original expression. This preserves the objective's minimisers mathematically but can change numerical trajectories relative to older runs. Unsuccessful or infeasible solves stop with an explicit error instead of being classified as bistable.

## 3. Build a unit library

Start with a single-case check:

```matlab
setup_project
cfg = library_config;
cfg.alpha2 = pi/3;
cfg.alpha3 = pi/3;
cfg.beta = pi/40;
cfg.postprocess = false;
cfg.refine = false;
result = main_unit_library(cfg);
```

For a full sweep, use `cfg = library_config` without these overrides. It samples ordered target angles with `a1 = pi - a2 - a3` and `a1 >= a2 >= a3`, and 20 beta values from 0 to pi/20. Angles are in radians. The original sweep settings use flank-length ratio **0.80**, while the supplied historical table's refinement metadata uses **0.85**. Set `cfg.flank_length_ratio = 0.85` for a study intended to align with that geometry; this still does not promise bit-for-bit regeneration of the historical table.

The default sweep runs serially. Set `cfg.use_parallel = true` to use Parallel Computing Toolbox. It can be expensive; reduce the angle and beta arrays while developing a new study.

Each run saves:

- `settings.mat`: the full configuration.
- `beta_batches_refine/`: one raw result batch per completed beta.
- `anisotropy_study_refine.mat`: the raw angle/beta table.
- With `postprocess = true`: cleaned, filled and filtered tables, plus reports identifying removed/interpolated entries.
- `anisotropy_filter_refine.mat`: the final table and refinement settings, suitable as a new `cfg.library` input.
- `library_result.mat`: final table, configuration and MATLAB version.

Postprocessing removes isolated support regions and fills selected gaps by interpolation. Filled values are estimates, not additional solved samples. Endpoint refinement reruns the mechanics to move the detected state toward the prescribed endpoint; inspect its `correction_status` and iteration-history columns. The batches preserve completed work if a later solve fails, but automatic resume is not implemented.

The bundled library is never overwritten. A dense sweep, full-sheet validation and reproduction of all manuscript figures are outside the current validation record; see [the checks performed](validation.md).
