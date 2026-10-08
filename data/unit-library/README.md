# Selected unit library

`anisotropy_filter_refine.mat` is copied unchanged from the table selected by the original 3D driver: `Inverse_design/output/ternary_refine_20260325_161516/anisotropy_filter_refine.mat` at the baseline commit recorded in `docs/source-commit.txt`.

The MAT file contains the MATLAB table `anisotropy_filter_refine`:

| Column | Meaning |
|---|---|
| `a1`, `a2`, `a3` | Angles describing the directional stretch ratios, in radians; sum to pi |
| `eps_bist` | Sampled strain at the detected second stable configuration |
| `eta_val` | The model's normalized bistability measure |
| `beta` | Unit tilting angle, in radians |

The original surface driver filters to finite `eps_bist` and `eta_val > 0.10`. Keep this source table fixed when comparing runs. Library-generation routines are under `src/inverse_design/`; their current settings are not claimed to regenerate this historical table bit for bit.
