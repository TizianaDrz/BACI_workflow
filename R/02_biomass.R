# ==============================================================================
# Zooplankton: abundance -> dry-mass biomass, biogeographic affinity,
# holoplankton / meroplankton flag
# ==============================================================================

# Conversion table: semicolon-separated with columns
#   taxon; biogeographicOrigin; conversionFactorForBiomassEstimation (mg DW ind-1);
#   planktonType (holoplankton / meroplankton)
# If planktonType is missing, all taxa are treated as holoplankton.
read_conversion_table <- function(path) {
  raw <- read_delim(path, delim = ";", locale = locale(encoding = "UTF-8"),
                    show_col_types = FALSE, name_repair = "minimal")
  if (!"planktonType" %in% names(raw)) {
    warning("No planktonType column in the conversion table: all taxa treated as holoplankton.")
    raw$planktonType <- "holoplankton"
  }
  raw %>%
    filter(!is.na(taxon)) %>%
    transmute(
      taxon        = str_squish(taxon),
      affinity     = biogeographicOrigin,
      dw_mg_ind    = suppressWarnings(as.numeric(conversionFactorForBiomassEstimation)),
      planktonType = coalesce(str_to_lower(str_squish(planktonType)), "holoplankton")
    )
}

# Link each species + stage in the counts to one conversion factor, affinity and plankton type.
# lookup: species, stage, taxon. `taxon` must match the conversion table; several
# taxa separated by "|" are averaged (e.g. all stages of a taxon counted without stage).
resolve_lookup <- function(lookup, conversion) {
  res <- lookup %>%
    mutate(stage = coalesce(as.character(stage), ""),
           conv_taxon = str_split(taxon, "\\|")) %>%
    unnest(conv_taxon) %>%
    mutate(conv_taxon = str_squish(conv_taxon)) %>%
    left_join(conversion, by = c("conv_taxon" = "taxon"))

  missing <- res %>% filter(is.na(dw_mg_ind)) %>% distinct(conv_taxon)
  if (nrow(missing) > 0) {
    warning("Not in the conversion table (or no numeric factor): ",
            paste(missing$conv_taxon, collapse = "; "))
  }

  out <- res %>%
    group_by(species, stage, taxon) %>%
    summarise(dw_mg_ind    = mean(dw_mg_ind),
              affinity     = paste(unique(affinity), collapse = "/"),
              planktonType = paste(unique(planktonType), collapse = "/"),
              .groups = "drop")

  mixed <- out %>% filter(str_detect(affinity, "/") | str_detect(planktonType, "/"))
  if (nrow(mixed) > 0) {
    warning("Averaged taxa with different affinity or plankton type (fix the lookup): ",
            paste(mixed$species, mixed$stage, collapse = "; "))
  }
  out
}

# Biomass per taxon and layer.
# abundance: station, min_depth, max_depth, species, stage, abundance_m3 (ind m-3)
# Taxa marked "meroplankton" in the conversion table get holoplankton = FALSE
# and are left out of BACI (see affinity_shares()).
compute_biomass <- function(abundance, lookup_resolved) {
  out <- abundance %>%
    mutate(stage = coalesce(as.character(stage), "")) %>%
    left_join(lookup_resolved, by = c("species", "stage")) %>%
    mutate(
      biomass_mg_m3 = abundance_m3 * dw_mg_ind,
      affinity      = coalesce(affinity, "Unknown"),
      holoplankton  = coalesce(planktonType, "holoplankton") != "meroplankton"
    )

  unmatched <- out %>% filter(abundance_m3 > 0, is.na(dw_mg_ind)) %>% distinct(species, stage)
  if (nrow(unmatched) > 0) {
    warning("Taxa with abundance > 0 but no conversion factor (excluded from biomass): ",
            paste(unmatched$species, unmatched$stage, collapse = "; "))
  }
  out
}
