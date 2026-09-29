# ==============================================================================
# Hydrography: CTD reading, TEOS-10 conversion, water mass classification and
# Atlantic Water Mass Score (AWMS) per net layer
# ==============================================================================

# AWMS per water mass (Table 1 of the manuscript)
awms_scores <- c(
  "Arctic Water (ArW)"               = 0,
  "Winter-Cooled Water (WCW)"        = 0.1,
  "Surface Water (SW)"               = 0.3,
  "Intermediate Water (IW)"          = 0.5,
  "transformed Atlantic Water (tAW)" = 0.8,
  "Atlantic Water (AW)"              = 1
)

# Read one Sea-Bird CTD file (.cnv or processed bin .txt) into a data frame.
# Sea-Bird bad-flag values (-9.990e-29) are set to NA.
read_ctd_sbe <- function(file, station) {
  ctd <- oce::read.oce(file)
  n <- length(ctd[["pressure"]])
  get_var <- function(...) {
    for (name in c(...)) {
      x <- ctd[[name]]
      if (!is.null(x)) {
        x <- as.numeric(x)
        return(if (length(x) == 1) rep(x, n) else x)
      }
    }
    rep(NA_real_, n)
  }
  df <- tibble(
    station     = station,
    pressure    = get_var("pressure", "prDM"),
    depth       = get_var("depth", "depSM"),
    temperature = get_var("temperature", "t090C"),
    salinity    = get_var("salinity", "sal00"),
    latitude    = get_var("latitude"),
    longitude   = get_var("longitude")
  )
  df <- mutate(df, across(where(is.numeric), ~ replace(.x, abs(.x) < 1e-20 & .x != 0, NA)))
  if (all(is.na(df$depth))) {
    df$depth <- oce::swDepth(df$pressure, mean(df$latitude, na.rm = TRUE))
  }
  df
}

# Practical salinity / in-situ temperature -> Absolute Salinity (SA) and Conservative Temperature (CT)
add_teos10 <- function(ctd) {
  ctd %>%
    mutate(
      SA = gsw::gsw_SA_from_SP(salinity, pressure, longitude, latitude),
      CT = gsw::gsw_CT_from_t(SA, temperature, pressure)
    )
}

# Water mass classification after Skogseth et al. (2020), with the extensions used
# in the manuscript to close gaps in the original table (see README):
#   - ArW: T < 1 (original <= 0) and SA 34.16-34.97 (original 34.46-34.97)
#   - IW:  SA up to 35.07 (original < 34.87) for T > 1
#   - AW:  also 1 <= T < 3 with SA > 35.07 (gap between TAW and AW)
#   - Local Water (LW) is not used
# Rules are evaluated top to bottom; the first match wins.
classify_water_mass <- function(CT, SA) {
  case_when(
    CT >= 3 & SA >= 35.07                         ~ "Atlantic Water (AW)",
    CT >= 1 & CT < 3 & SA > 35.07                 ~ "Atlantic Water (AW)",
    CT < -0.5 & SA >= 34.56                       ~ "Winter-Cooled Water (WCW)",
    CT < 1 & SA >= 34.16 & SA <= 34.97            ~ "Arctic Water (ArW)",
    CT >= 1 & SA < 34.16                          ~ "Surface Water (SW)",
    CT >= 1 & CT <= 3 & SA >= 34.87 & SA <= 35.07 ~ "transformed Atlantic Water (tAW)",
    CT > 1 & SA >= 34.16 & SA <= 35.07            ~ "Intermediate Water (IW)",
    TRUE                                          ~ "Unknown"
  )
}

# Dominant water mass and AWMS for each net layer.
# ctd:    station, depth, water_mass (one row per CTD depth bin)
# layers: station, min_depth, max_depth (contiguous layers per station)
# CTD bins shallower than the first layer or deeper than the last one are
# assigned to the nearest layer (as in the manuscript).
layer_water_mass <- function(ctd, layers, scores = awms_scores) {
  out <- ctd %>%
    filter(!is.na(depth), !is.na(water_mass)) %>%
    group_by(station) %>%
    group_modify(function(d, key) {
      lay <- layers %>% filter(station == key$station) %>% arrange(min_depth)
      if (nrow(lay) == 0) return(tibble())
      breaks <- c(lay$min_depth, max(lay$max_depth))
      idx <- findInterval(d$depth, breaks, all.inside = TRUE)
      mutate(d, min_depth = lay$min_depth[idx], max_depth = lay$max_depth[idx])
    }) %>%
    group_by(station, min_depth, max_depth) %>%
    summarise(
      water_mass = names(sort(table(water_mass), decreasing = TRUE))[1],
      n_ctd_bins = n(),
      .groups = "drop"
    ) %>%
    mutate(awms = unname(scores[water_mass]))

  if (any(is.na(out$awms))) {
    warning("Layers without an AWMS score (dominant water mass not in `scores`): ",
            paste(out$station[is.na(out$awms)], out$min_depth[is.na(out$awms)], sep = " ", collapse = ", "))
  }
  missing <- anti_join(layers, out, by = c("station", "min_depth", "max_depth"))
  if (nrow(missing) > 0) {
    warning("Net layers without CTD data: ",
            paste(missing$station, missing$min_depth, missing$max_depth, sep = " ", collapse = ", "))
  }
  out
}
