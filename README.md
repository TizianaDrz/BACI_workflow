# BACI – Borealization–Atlantification Coupling Index

R workflow to calculate the Borealization–Atlantification Coupling Index (BACI) from CTD profiles and depth-stratified zooplankton abundance, as described in:

> Durazzano T., Ntinou I. V., Skjæveland B. H., Kwasniewski S., Varpe Ø., Hop H., Daase M., Søreide J. E. (in prep.) Borealization–Atlantification Coupling Index used to assess ecosystem change and its effects on *Calanus* size in Svalbard fjords.

This repository documents **how the index is constructed** and provides general, reusable code to calculate it from your own data. It contains **no data**, so the manuscript results cannot be re-run from it directly. See [data/README.md](data/README.md) for the input formats.

## What BACI measures

BACI shows, for each net layer, whether the **water** and the **zooplankton** point the same way along the Arctic–Atlantic gradient:

- **AWMS** (Atlantic Water Mass Score): how Atlantic the water is. The dominant water mass in the layer gets a fixed score:

  | Water mass | AWMS |
  |---|---|
  | Arctic Water (ArW) | 0.0 |
  | Winter-Cooled Water (WCW) | 0.1 |
  | Surface Water (SW) | 0.3 |
  | Intermediate Water (IW) | 0.5 |
  | Transformed Atlantic Water (TAW) | 0.8 |
  | Atlantic Water (AW) | 1.0 |

- **Bmid**: how Atlantic the zooplankton is, as the share of holoplankton biomass from taxa with Atlantic affinity. Taxa found in both regions ("Both") make this uncertain, so two bounds are calculated:
  - Bmin = Atlantic / (Atlantic + Arctic + Both), counting "Both" as Arctic
  - Bmax = (Atlantic + Both) / (Atlantic + Arctic + Both), counting "Both" as Atlantic
  - Bmid = (Bmin + Bmax) / 2

  Unclassified taxa and meroplankton are excluded.

When the two components agree, the layer is *coupled*; when they disagree (e.g. Atlantic Water with an Arctic community), it is *decoupled*. Always report AWMS and Bmid next to BACI, because the same BACI can come from very different combinations.

### Two versions of the index

| | `BACI_relative` (manuscript) | `BACI_unscaled` |
|---|---|---|
| Formula | z-score AWMS and Bmid across all layers, average 50/50, rescale to 0–1 | 0.5 × AWMS + 0.5 × Bmid |
| Scale | Relative to the dataset: lowest layer = 0, highest = 1 | Absolute: both components are already on fixed 0–1 scales |
| Use for | Positioning layers/fjords within **one** survey | Comparing **years, cruises or studies** |
| Changes when data are added? | Yes | No |

⚠️ **Values of `BACI_relative` from different datasets cannot be compared.** For example, Kongsfjorden 2024 is 0.88 on the manuscript scale and 0.65 unscaled. Use `BACI_unscaled` for comparisons over time.

Fjord values (BACI_Fjord, Eq. 4) are the layer values weighted by layer thickness: Σ(BACI × Δz) / ΣΔz.

## Workflow

```
CTD files ──► TEOS-10 (SA, CT) ──► water mass per depth bin ──► dominant water mass per net layer ──► AWMS
                                                                                                        │
abundance (ind m⁻³) ──► taxon lookup ──► dry mass (mg DW m⁻³) + affinity ──► Bmin / Bmid / Bmax ─────────┤
                                                                                                        ▼
                                                                              BACI per layer ──► BACI_Fjord
```

| Step | Function | File |
|---|---|---|
| Read Sea-Bird CTD file | `read_ctd_sbe()` | `R/01_water_masses.R` |
| Absolute Salinity, Conservative Temperature | `add_teos10()` | `R/01_water_masses.R` |
| Water mass classification | `classify_water_mass()` | `R/01_water_masses.R` |
| Dominant water mass + AWMS per net layer | `layer_water_mass()` | `R/01_water_masses.R` |
| Conversion table | `read_conversion_table()` | `R/02_biomass.R` |
| Match counts to conversion taxa | `resolve_lookup()` | `R/02_biomass.R` |
| Biomass and meroplankton flag | `compute_biomass()` | `R/02_biomass.R` |
| Bmin / Bmid / Bmax | `affinity_shares()` | `R/03_baci.R` |
| BACI per layer (both versions) | `baci_layers()` | `R/03_baci.R` |
| BACI_Fjord | `baci_fjord()` | `R/03_baci.R` |
| 25/75 weighting sensitivity + Spearman ρ | `weighting_sensitivity()` | `R/03_baci.R` |

## How to run

1. Open `BACI_workflow.Rproj` in RStudio (or set the working directory to this folder).
2. Install the packages once:
   ```r
   install.packages(c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "oce", "gsw"))
   ```
3. Put your input files in `data/`, following [data/README.md](data/README.md) and the examples in `templates/`. The abundance file must be in long format (one row per station × layer × species × stage); reshape your own sheet to this first.
4. Copy `templates/taxon_lookup.csv` to `data/` and adapt it to your taxa names. The conversion table is already included (`config/`); check its `planktonType` column.
5. Run:
   ```r
   source("run_baci.R")
   ```
   Results are written to `outputs/`:
   - `ctd_water_masses.csv`
   - `biomass_by_taxon.csv`
   - `baci_layers.csv`
   - `baci_fjord.csv`
   - `baci_weighting_sensitivity.csv`

To check that the code works, run `source("tests/test_baci.R")`. It uses small made-up data with known answers and should end with "All checks passed".

## Choices you should know about

**Water mass classification.** Skogseth et al. (2020), with extensions to close gaps in the original table (manuscript Table S1). Thresholds are in Conservative Temperature (°C) and Absolute Salinity (g kg⁻¹):

| Water mass | Original | Used here |
|---|---|---|
| AW | T ≥ 3, S ≥ 35.07 | same, **plus 1 ≤ T < 3 with S > 35.07** |
| ArW | T ≤ 0, S 34.46–34.97 | **T < 1, S 34.16–34.97** |
| SW | T ≥ 1, S < 34.16 | same |
| WCW | T < −0.5, S ≥ 34.56 | same |
| IW | T > 1, S 34.16–34.87 | **T > 1, S 34.16–35.07** |
| TAW | 1 ≤ T ≤ 3, S 34.87–35.07 | same |
| LW | T < 1, S ≥ 32.15 | not used |

Rules are applied in the order listed in `classify_water_mass()`; the first match wins. Anything else is "Unknown" (no AWMS).

**Conversion factors.** `config/zooplankton_dry_mass_conversions.csv` gives, for 216 taxa and stages or size classes, the dry mass per individual (mg DW ind⁻¹), the biogeographic affinity and whether the taxon is holo- or meroplankton, with the source of each value in the `bibliographicCitation` column (mainly Hop et al. 2019b, *Zooplankton in Kongsfjorden (1996–2016) in relation to climate change*, and Wold et al. 2025; see the column for each value). Assumed values are explained in `occurrenceRemarks`. Please cite the original sources when you use it. Where counts have no stage or size, a rule is needed. `templates/taxon_lookup.csv` uses:
- the adult factor (AF/CVI) for copepods counted without stage;
- the smallest size class for pteropods and *Themisto*;
- several taxa separated by `|`, which averages their factors.

State the rules you use. Large taxa such as pteropods can dominate biomass.

**Meroplankton** (larvae of benthic taxa) are excluded from BACI. Which taxa count as meroplankton is set in the `planktonType` column of the conversion table (holoplankton / meroplankton), so it travels with each taxon. Change it there if you classify a taxon differently.

**Layers.** Water masses are summarised over the net layers, so coarse layers merge thin water masses. For example, a 0–50 m layer can hide a 0–20 m Surface Water layer. Compare datasets with similar layer structure where possible.

**Uncertainty.** Differences smaller than about ±0.05 in BACI are usually within the Bmin–Bmax range.

## Folder structure

```
BACI_workflow/
├── README.md
├── LICENSE
├── run_baci.R              main script
├── R/                      functions
├── config/zooplankton_dry_mass_conversions.csv   dry mass, affinity and holo/meroplankton per taxon
├── templates/              example input files (stations, abundance, taxon lookup)
├── tests/test_baci.R       checks on made-up data
├── data/                   your input files (ignored by git)
└── outputs/                results (ignored by git)
```

## Versions

| Version | Date | Changes |
|---|---|---|
| v1.0 | 2026 | Version as submitted with the manuscript. |

## Citation

If you use this workflow, please cite the manuscript above and this repository (DOI will be added on release via Zenodo).

## Licence

The code is released under the MIT licence (see `LICENSE`).

The conversion factors and biogeographic affinities in `config/zooplankton_dry_mass_conversions.csv` are compiled from published sources; the source of each value is given in the `bibliographicCitation` column. When you use them, please cite those original publications.
