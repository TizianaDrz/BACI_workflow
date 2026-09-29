# ==============================================================================
# Checks for the BACI workflow, on small made-up data (no real data needed).
# Run from the BACI_workflow folder:  source("tests/test_baci.R")
# ==============================================================================

if (!dir.exists("R") && dir.exists("BACI_workflow/R")) setwd("BACI_workflow")
if (!dir.exists("R")) stop("Run this from the BACI_workflow folder (open BACI_workflow.Rproj, or setwd()).")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

check <- function(desc, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", desc, call. = FALSE)
  cat("ok  -", desc, "\n")
}
near <- function(a, b, tol = 1e-9) all(abs(a - b) < tol)

# ---- Water mass classification --------------------------------------------------
check("AW (warm, saline)",            classify_water_mass(4, 35.2)    == "Atlantic Water (AW)")
check("AW gap rule (1-3 C, SA>35.07)", classify_water_mass(2, 35.15)  == "Atlantic Water (AW)")
check("TAW",                          classify_water_mass(2, 34.95)   == "transformed Atlantic Water (tAW)")
check("WCW",                          classify_water_mass(-1.8, 34.9) == "Winter-Cooled Water (WCW)")
check("ArW extension (T 0-1)",        classify_water_mass(0.5, 34.5)  == "Arctic Water (ArW)")
check("SW",                           classify_water_mass(6, 32.5)    == "Surface Water (SW)")
check("IW extension (T>3, SA 34.87-35.07)", classify_water_mass(4, 34.95) == "Intermediate Water (IW)")
check("cold fresh water is Unknown",  classify_water_mass(0, 33)      == "Unknown")

# ---- Dominant water mass per layer ----------------------------------------------------
toy_ctd <- tibble(
  station    = "A",
  depth      = c(5, 20, 40, 60, 70, 80, 90, 120),
  water_mass = c("Surface Water (SW)", "Intermediate Water (IW)", "Intermediate Water (IW)",
                 "Atlantic Water (AW)", "Atlantic Water (AW)", "Atlantic Water (AW)",
                 "Intermediate Water (IW)", "Atlantic Water (AW)")   # 120 m is below the last layer
)
toy_layers <- tibble(station = "A", min_depth = c(0, 50), max_depth = c(50, 100))
wm <- layer_water_mass(toy_ctd, toy_layers)
check("dominant water mass per layer", identical(wm$water_mass, c("Intermediate Water (IW)", "Atlantic Water (AW)")))
check("AWMS scores",                   near(wm$awms, c(0.5, 1)))
check("CTD bins below last layer go to last layer", wm$n_ctd_bins[2] == 5)

# ---- Lookup and biomass ------------------------------------------------------------------
toy_conv <- tibble(
  taxon     = c("Sp atl", "Sp arc", "Sp both CI", "Sp both CII", "Larva"),
  affinity  = c("Atlantic", "Arctic", "Both", "Both", "Unknown"),
  dw_mg_ind = c(1, 2, 0.1, 0.3, 5),
  planktonType = c("holoplankton", "holoplankton", "holoplankton", "holoplankton", "meroplankton")
)
toy_lookup <- tibble(
  species = c("Sp atl", "Sp arc", "Sp both", "Larva"),
  stage   = c("", "", "", ""),
  taxon   = c("Sp atl", "Sp arc", "Sp both CI|Sp both CII", "Larva")
)
lk <- resolve_lookup(toy_lookup, toy_conv)
check("'|' in lookup averages factors", near(lk$dw_mg_ind[lk$species == "Sp both"], 0.2))

toy_ab <- tibble(
  station = "A", min_depth = c(0, 0, 0, 0, 50, 50, 50),  max_depth = c(50, 50, 50, 50, 100, 100, 100),
  species = c("Sp atl", "Sp arc", "Sp both", "Larva", "Sp atl", "Sp arc", "Sp both"),
  stage   = "",
  abundance_m3 = c(1, 0.5, 10, 100, 3, 0.5, 0)
)
bm <- compute_biomass(toy_ab, lk)
check("meroplankton flagged from the conversion table", !bm$holoplankton[bm$species == "Larva"][1])
check("holoplankton kept", all(bm$holoplankton[bm$species != "Larva"]))

sh <- affinity_shares(bm)
# Layer 0-50: Atl 1, Arc 1, Both 2 -> Bmin 0.25, Bmax 0.75, Bmid 0.5 (Larva excluded)
# Layer 50-100: Atl 3, Arc 1, Both 0 -> Bmin 0.75, Bmax 0.75, Bmid 0.75
check("biomass per affinity", near(c(sh$B_atl[1], sh$B_arc[1], sh$B_both[1]), c(1, 1, 2)))
check("Bmin/Bmax/Bmid (Eq. 1-3)", near(c(sh$Bmin[1], sh$Bmax[1], sh$Bmid[1]), c(0.25, 0.75, 0.5)))

# ---- BACI --------------------------------------------------------------------------------
# Add a second station so the relative index has something to scale against
awms2 <- bind_rows(wm %>% select(station, min_depth, max_depth, water_mass, awms),
                   tibble(station = "B", min_depth = c(0, 50), max_depth = c(50, 150),
                          water_mass = c("Surface Water (SW)", "Winter-Cooled Water (WCW)"), awms = c(0.3, 0.1)))
sh2 <- bind_rows(sh, tibble(station = "B", min_depth = c(0, 50), max_depth = c(50, 150),
                            B_atl = c(1, 1), B_arc = c(1, 3), B_both = c(0, 0), B_unknown = 0) %>%
                   add_atlantic_shares())
bl <- baci_layers(awms2, sh2)

check("unscaled BACI = 0.5*AWMS + 0.5*Bmid", near(bl$BACI_unscaled, 0.5 * bl$awms + 0.5 * bl$Bmid))
check("relative BACI spans exactly 0-1", near(range(bl$BACI_relative), c(0, 1)))
check("unscaled min <= mid <= max", all(bl$BACI_unscaled_min <= bl$BACI_unscaled & bl$BACI_unscaled <= bl$BACI_unscaled_max))

fj <- baci_fjord(bl)
b  <- bl %>% filter(station == "B")
check("BACI_Fjord is thickness-weighted (Eq. 4)",
      near(fj$BACI_unscaled[fj$station == "B"], sum(b$BACI_unscaled * c(50, 100)) / 150))

ws <- weighting_sensitivity(bl)
check("weighting: 50/50 equals relative BACI", near(ws$layers$BACI_w50, bl$BACI_relative))

cat("\nAll checks passed.\n")
