# Input data

Put your own files here. Everything in this folder except this README is ignored by git, so no data are published by accident.

Examples of each format are in `../templates/`.

| File | Format | Columns |
|---|---|---|
| `stations.csv` | comma-separated | `station`, `ctd_file` (file name of the CTD cast in this folder) |
| CTD files | Sea-Bird `.cnv` or processed bin `.txt` | pressure, temperature, salinity, depth, latitude, longitude (read with `oce`) |
| `abundance.csv` | comma-separated, long format | `station`, `min_depth`, `max_depth`, `species`, `stage` (may be empty), `abundance_m3` (ind m⁻³) |
| `taxon_lookup.csv` | comma-separated | `species`, `stage`, `taxon` (name in the conversion table; several separated by `\|` are averaged) |

The conversion table is not needed here: it is included in `../config/zooplankton_dry_mass_conversions.csv`. To use your own, keep the same format (**semicolon**-separated; `taxon`, `biogeographicOrigin` = Arctic / Atlantic / Both / Unknown, `conversionFactorForBiomassEstimation` in mg DW ind⁻¹, `planktonType` = holoplankton / meroplankton) and change the path in `run_baci.R`. Without a `planktonType` column, all taxa are treated as holoplankton.

Notes:
- `station` must be spelled the same in `stations.csv` and `abundance.csv`.
- Net layers must be contiguous per station (e.g. 0–50, 50–100, 100–265).
- Abundance sheets come in many layouts; reshape yours to the long format above (one row per station × layer × species × stage) before running, e.g. with `tidyr::pivot_longer()`.
