# ==============================================================================
# BACI workflow: CTD + multinet abundance -> water masses, biomass, BACI
#
# Run from the BACI_workflow folder (open BACI_workflow.Rproj in RStudio, or setwd()).
# Put your input files in data/ (see data/README.md and templates/).
# Results are written to outputs/.
# ==============================================================================

if (!dir.exists("R") && dir.exists("BACI_workflow/R")) setwd("BACI_workflow")
if (!dir.exists("R")) stop("Run this from the BACI_workflow folder (open BACI_workflow.Rproj, or setwd()).")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

# ---- 1. Settings --------------------------------------------------------------
paths <- list(
  stations     = "data/stations.csv",          # station, ctd_file
  abundance    = "data/abundance.csv",         # station, min_depth, max_depth, species, stage, abundance_m3
  lookup       = "data/taxon_lookup.csv",      # species, stage, taxon (start from templates/)
  conversion   = "config/zooplankton_dry_mass_conversions.csv"   # included; also marks meroplankton
)
w_bio <- 0.5   # weight of the biological component (manuscript: 0.5)

dir.create("outputs", showWarnings = FALSE)

# ---- 2. Hydrography: water masses and AWMS per net layer ----------------------
stations <- read_csv(paths$stations, show_col_types = FALSE)

ctd <- map2_dfr(file.path("data", stations$ctd_file), stations$station, read_ctd_sbe) %>%
  add_teos10() %>%
  mutate(water_mass = classify_water_mass(CT, SA))

abundance <- read_csv(paths$abundance, col_types = cols(stage = col_character()))
layers    <- abundance %>% distinct(station, min_depth, max_depth)

awms <- layer_water_mass(ctd, layers)

# ---- 3. Biomass and biogeographic affinity ------------------------------------
conversion <- read_conversion_table(paths$conversion)
lookup     <- read_csv(paths$lookup, col_types = cols(.default = col_character()))

biomass <- compute_biomass(abundance, resolve_lookup(lookup, conversion))
shares  <- affinity_shares(biomass)

# ---- 4. BACI ------------------------------------------------------------------
baci       <- baci_layers(awms, shares, w_bio = w_bio)
baci_fj    <- baci_fjord(baci)
weighting  <- weighting_sensitivity(baci)

print(baci %>% select(station, min_depth, max_depth, water_mass, awms, Bmin, Bmid, Bmax,
                      BACI_unscaled, BACI_relative))
print(baci_fj)
print(weighting$spearman)

# ---- 5. Outputs ---------------------------------------------------------------
write_excel_csv(ctd,              "outputs/ctd_water_masses.csv")
write_excel_csv(biomass,          "outputs/biomass_by_taxon.csv")
write_excel_csv(baci,             "outputs/baci_layers.csv")
write_excel_csv(baci_fj,          "outputs/baci_fjord.csv")
write_excel_csv(weighting$layers, "outputs/baci_weighting_sensitivity.csv")
