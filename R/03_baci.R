# ==============================================================================
# Borealization-Atlantification Coupling Index (BACI)
# ==============================================================================

# Biomass per affinity and Atlantic share (Eq. 1-3) per station x layer, holoplankton only.
# Unclassified ("Unknown") taxa are reported but excluded from Bmin/Bmax/Bmid.
affinity_shares <- function(biomass) {
  biomass %>%
    filter(holoplankton) %>%
    group_by(station, min_depth, max_depth) %>%
    summarise(
      B_atl     = sum(biomass_mg_m3[affinity == "Atlantic"], na.rm = TRUE),
      B_arc     = sum(biomass_mg_m3[affinity == "Arctic"],   na.rm = TRUE),
      B_both    = sum(biomass_mg_m3[affinity == "Both"],     na.rm = TRUE),
      B_unknown = sum(biomass_mg_m3[!affinity %in% c("Atlantic", "Arctic", "Both")], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    add_atlantic_shares()
}

add_atlantic_shares <- function(df) {
  df %>%
    mutate(
      Bmin = B_atl / (B_atl + B_arc + B_both),            # "Both" counted as Arctic
      Bmax = (B_atl + B_both) / (B_atl + B_arc + B_both), # "Both" counted as Atlantic
      Bmid = (Bmin + Bmax) / 2
    )
}

# BACI per layer, in two versions:
#   BACI_unscaled = (1 - w_bio) * AWMS + w_bio * Bmid
#       Absolute: comparable between years, cruises and studies.
#   BACI_relative = z-score AWMS and Bmid across all layers, combine, rescale 0-1
#       Manuscript version: positions layers relative to each other within ONE dataset;
#       values change whenever layers are added or removed.
# *_min / *_max use Bmin / Bmax instead of Bmid (classification sensitivity).
# The relative min/max are standardised with the Bmid mean/sd and rescaled with the
# main index range, as in the manuscript.
baci_layers <- function(awms_layers, shares, w_bio = 0.5) {
  d <- inner_join(awms_layers, shares, by = c("station", "min_depth", "max_depth"))

  dropped <- d %>% filter(is.na(awms) | is.na(Bmid))
  if (nrow(dropped) > 0) {
    warning("Layers dropped (no AWMS or no classified biomass): ",
            paste(dropped$station, dropped$min_depth, sep = " ", collapse = ", "))
  }
  d <- d %>% filter(!is.na(awms), !is.na(Bmid))

  z <- function(x, ref) (x - mean(ref)) / sd(ref)
  rescale_to <- function(x, ref) (x - min(ref)) / (max(ref) - min(ref))

  w_env <- 1 - w_bio
  s     <- w_env * z(d$awms, d$awms) + w_bio * z(d$Bmid, d$Bmid)
  s_min <- w_env * z(d$awms, d$awms) + w_bio * z(d$Bmin, d$Bmid)
  s_max <- w_env * z(d$awms, d$awms) + w_bio * z(d$Bmax, d$Bmid)

  d %>%
    mutate(
      thickness         = max_depth - min_depth,
      BACI_unscaled     = w_env * awms + w_bio * Bmid,
      BACI_unscaled_min = w_env * awms + w_bio * Bmin,
      BACI_unscaled_max = w_env * awms + w_bio * Bmax,
      BACI_relative     = rescale_to(s, s),
      BACI_relative_min = rescale_to(s_min, s),
      BACI_relative_max = rescale_to(s_max, s)
    ) %>%
    arrange(station, min_depth)
}

# BACI_Fjord (Eq. 4): layer values weighted by layer thickness
baci_fjord <- function(layers) {
  layers %>%
    group_by(station) %>%
    summarise(
      sampled_depth = paste0(min(min_depth), "-", max(max_depth)),
      across(c(awms, Bmin, Bmid, Bmax, starts_with("BACI_")),
             ~ sum(.x * thickness) / sum(thickness)),
      .groups = "drop"
    )
}

# Weighting sensitivity of the relative BACI (manuscript Fig. S1B / Table S3):
# 25% and 75% biological weight, rescaled with the 50/50 range; Spearman rho between schemes.
weighting_sensitivity <- function(layers) {
  z <- function(x) (x - mean(x)) / sd(x)
  s50 <- 0.5 * z(layers$awms) + 0.5 * z(layers$Bmid)
  rescale_to <- function(x) (x - min(s50)) / (max(s50) - min(s50))
  out <- layers %>%
    mutate(
      BACI_w25 = rescale_to(0.75 * z(awms) + 0.25 * z(Bmid)),
      BACI_w50 = rescale_to(s50),
      BACI_w75 = rescale_to(0.25 * z(awms) + 0.75 * z(Bmid))
    )
  list(
    layers   = out,
    spearman = cor.test(out$BACI_w75, out$BACI_w25, method = "spearman", exact = FALSE)
  )
}
