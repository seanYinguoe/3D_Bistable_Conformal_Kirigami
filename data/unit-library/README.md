# Supplied unit library

`anisotropy_filter_refine.mat` is unchanged from the table selected by the original 3D driver:

`Xiaoyuan_bistable_kirigami_unit/Inverse_design/output/ternary_refine_20260325_161516/anisotropy_filter_refine.mat`

Source commit: `afb8ffa`. SHA-256: `e4d1f140f81d38bc0348ea34fc16b924faefe981d796b7fbd16b28ade47b7cfe`.

The MAT file contains the 1,590-row table `anisotropy_filter_refine` and `opts_refine` metadata.

| Column | Meaning |
|---|---|
| `a1`, `a2`, `a3` | Ordered target-triangle angles describing directional stretch ratios; radians, sum to pi |
| `eps_bist` | Library estimate of the reference-edge strain for the deployed state |
| `eta_val` | Normalised energy-barrier measure |
| `beta` | Unit tilt in radians |
| `*_old`, `*_new` | Values before and after endpoint refinement |
| `alpha_bist`, `correction_status` | Detected deployment fraction and refinement outcome |
| `alpha_iter_history`, `eps_iter_history` | Recorded refinement iterations |

The saved refinement geometry uses edge length 15, flank length 12.75, flank width 0.75 and ligament width 0.225 (ratios 0.85, 0.05 and 0.015). `main_surface` filters to finite `eps_bist` and eta above 0.10 by default.

This is a historical, processed numerical table: not every entry is a freshly solved unit or a directly observed energy minimum. Inspect its correction statuses. The current sweep configuration retains a flank-length ratio of 0.80 from the source sweep script; matching this table's geometry requires explicitly choosing 0.85. A full regeneration and bit-for-bit equivalence have not been established.

`main_unit_library` always writes a new table to `results/`. Set `cfg.library` in a surface run to use a reviewed replacement.
