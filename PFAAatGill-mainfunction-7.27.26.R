##########################################################################################
# Last edit: 2026-07-26
##########################################################################################
rm(list = ls())

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(tictoc)

setwd("/Users/huzhiji/OneDrive - Harvard University/Research/PFAAatGill/")

# =========================================================
# SETTINGS
# =========================================================
excel_path <- "PFAA_inputs.xlsx"
sheet_name <- "Properties"
out_dir    <- "PBTK_lab_outputs"
n_iter     <- 10000
verbose    <- TRUE

# Choose uptake mode: "lab" or "natural"
uptake_mode <- "lab"    

# Flexible setting
initial_concentration <- 1e-3 # ug/cm3
days_uptake <- 28 # day
days_depuration <- 28  # day
target_pfas <- "PFOS" # default = NULL

# =========================================================
# 1) Model constants 
# =========================================================
build_constants <- function() {
  const <- list(
    # Simulation
    days_uptake     = days_uptake,  # day
    days_depuration = days_depuration,  # day
    dt_seconds      = 100, # second
    
    # Fish physiology
    FW = 8,         # g, body weight
    g  = 4.9e-3,    # g/day, linear growth rate constant
    
    # Membrane and transport
    A_gill   = 13,    # cm2, secondary lamallel surface area
    A_tissue = 166,   # cm2, main organ surface area
    X_mem    = 1e-4,  # cm, gill epithelial membrane thickens
    X_ABL    = 0.016, # cm, aqueous boundary layer thickness
    
    # Volume fractions (unitless)
    Vf_gill    = 0.02385,
    Vf_blood   = 0.365,
    Vf_liver   = 0.0114,
    Vf_kidney  = 0.0084,
    Vf_muscle  = 0.465,
    Vf_adipose = 0.0089,
    
    # Blood composition (unitless)
    Vf_lipid_blood             = 0.0092,
    Vf_phospholipids_blood     = 0.0092,
    Vf_protein_function_blood  = 0.0292051,
    Vf_protein_structure_blood = 0.0119949,
    Vf_nonlipid_blood          = 0.00433,
    Vf_water_blood             = 0.93697,
    
    # Water conditions
    V_water     = 45000,  # cm3 (45 L)
    C_DOC       = 1e-6,   # kg/L
    DOC_density = 1.8,    # kg DOC / L DOC
    k_des       = 0.1,    # /s
    D_DOC       = 4.765e-06, # cm/s 
    
    # Exposure water conconcentration
    C_water_default = initial_concentration, # ug/cm3
    
    # Tissue composition for partitioning (K_tissue/w)
    tissue = list(
      gill    = list(Vf_NL=0.002949197, Vf_PL=0.002280803, Vf_FP=0.0408,      Vf_SP=0.1632,      Vf_NLOM=0.00516, Vf_W=0.78561),
      liver   = list(Vf_NL=0.00400058,  Vf_PL=0.00352942,  Vf_FP=0.032985385, Vf_SP=0.090314615, Vf_NLOM=0.035,   Vf_W=0.83417),
      kidney  = list(Vf_NL=0.004884,    Vf_PL=0.001716,    Vf_FP=0.049571709, Vf_SP=0.135728291, Vf_NLOM=0.0018,  Vf_W=0.8063),
      muscle  = list(Vf_NL=0,           Vf_PL=0,           Vf_FP=0,           Vf_SP=1,           Vf_NLOM=0,       Vf_W=0),
      adipose = list(Vf_NL=1,           Vf_PL=0,           Vf_FP=0,           Vf_SP=0,           Vf_NLOM=0,       Vf_W=0),
      blood   = list(Vf_NL=0.0092,      Vf_PL=0.0092,      Vf_FP=0.029205063, Vf_SP=0.011994937, Vf_NLOM=0.00433, Vf_W=0.93607)
    ),
    
    # Gill membrane composition
    Vf_ph = 0.6,
    Vf_fp = 0.4
  )
  
  const$Vf_tissue <- 1 - const$Vf_blood
  const
}

# =========================================================
# 2) Import PFAA properties 
# =========================================================
standardize_pfas_table <- function(df) {
  
  df2 <- df %>%
    rename(
      chemical_name = chemical_name,
      MW            = MW,
      PFAA_density  = PFAA_density,
      
      logKow_mean   = logKOW_mean,
      logKow_sd     = logKOW_sd,
      
      logKmw_mean   = logKPLW_mean,
      logKmw_sd     = logKPLW_sd,
      
      logKfpw_mean  = logKFPW_mean,
      logKfpw_sd    = logKFPW_sd,
      
      logKspw_mean  = logKSPW_mean,
      logKspw_sd    = logKSPW_sd,
      
      logKnw_mean   = logKNLOMW_mean,
      logKnw_sd     = logKNLOMW_sd,
      
      logKoc_mean   = logKOC_mean,
      logKoc_sd     = logKOC_sd,
      
      k_elim_kidney = k_elim_kidney
    )
  
  required <- c(
    "chemical_name","MW","PFAA_density",
    "logKow_mean","logKow_sd",
    "logKmw_mean","logKmw_sd",
    "logKfpw_mean","logKfpw_sd",
    "logKspw_mean","logKspw_sd",
    "logKnw_mean","logKnw_sd",
    "logKoc_mean","logKoc_sd",
    "k_elim_kidney"
  )
  
  missing <- setdiff(required, names(df2))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
  
  df2
}

read_pfas_table <- function(excel_path, sheet_name = 1) {
  df_raw <- read.xlsx(excel_path, sheet = sheet_name)
  standardize_pfas_table(df_raw)
}

# =========================================================
# 3) Monte Carlo sampling of partition parameters 
# =========================================================
sample_params <- function(pfas_row, const) {
  
  logKow  <- rnorm(1, mean = pfas_row$logKow_mean,  sd = pfas_row$logKow_sd)
  logKmw  <- rnorm(1, mean = pfas_row$logKmw_mean,  sd = pfas_row$logKmw_sd)
  logKfpw <- rnorm(1, mean = pfas_row$logKfpw_mean, sd = pfas_row$logKfpw_sd)
  logKspw <- rnorm(1, mean = pfas_row$logKspw_mean, sd = pfas_row$logKspw_sd)
  logKnw  <- rnorm(1, mean = pfas_row$logKnw_mean,  sd = pfas_row$logKnw_sd)
  logKoc  <- rnorm(1, mean = pfas_row$logKoc_mean,  sd = pfas_row$logKoc_sd)
  
  K_ow  <- 10^logKow
  K_mw  <- 10^logKmw
  K_fpw <- 10^logKfpw
  K_spw <- 10^logKspw
  K_nw  <- 10^logKnw
  K_oc  <- 10^logKoc
  
  K_tissue <- function(tiss) {
    K_ow * tiss$Vf_NL +
      K_mw * tiss$Vf_PL +
      K_fpw * tiss$Vf_FP +
      K_spw * tiss$Vf_SP +
      0.035 * K_nw * tiss$Vf_NLOM +
      tiss$Vf_W
  }
  
  K_gill_water    <- K_tissue(const$tissue$gill)
  K_liver_water   <- K_tissue(const$tissue$liver)
  K_kidney_water  <- K_tissue(const$tissue$kidney)
  K_muscle_water  <- K_tissue(const$tissue$muscle)
  K_adipose_water <- K_tissue(const$tissue$adipose)
  K_blood_water   <- K_tissue(const$tissue$blood)
  
  list(
    K_lipid_water        = K_ow,
    K_membrane_water     = K_mw,
    K_protein_blood      = K_fpw,
    K_protein_structural = K_spw,
    K_nonlipid_water     = K_nw,
    K_DOC_water          = K_oc,
    
    K_gill_water         = K_gill_water,
    K_liver_water        = K_liver_water,
    K_kidney_water       = K_kidney_water,
    K_muscle_water       = K_muscle_water,
    K_adipose_water      = K_adipose_water,
    K_blood_water        = K_blood_water,
    
    k_elim_kidney        = pfas_row$k_elim_kidney
  )
}

# =========================================================
# 4) Derived parameters
# =========================================================
calculate_derived_params <- function(sampled, pfas_row, const) {
  
  MW  <- pfas_row$MW
  rho <- pfas_row$PFAA_density
  
  K_mem <- sampled$K_membrane_water * const$Vf_ph + sampled$K_protein_blood * const$Vf_fp
  
  # L/L assume fish density = 1 kg/L
  K_tissue_water <- (
    sampled$K_liver_water    * const$Vf_liver +
      sampled$K_kidney_water * const$Vf_kidney +
      sampled$K_muscle_water * const$Vf_muscle +
      sampled$K_adipose_water * const$Vf_adipose +
      sampled$K_gill_water   * const$Vf_gill) / const$Vf_tissue
  
  # Diffusivities (cm2/s)
  D_aq     <- 10^(-4.13 - 0.453 * log10(MW)) * 1.348
  D_mem    <- D_aq / 10
  D_free   <- 1.52e-4 * (MW / rho)^(-0.64)
  
  # Permeabilities (cm/2)
  P_mem    <- (K_mem * D_mem) / const$X_mem
  P_ABL    <- (D_free + const$D_DOC * sampled$K_DOC_water * const$C_DOC) / const$X_ABL
  P_gill   <- 1 / (1 / P_ABL + 1 / P_mem)
  P_tissue <- (sampled$K_membrane_water * D_mem) / const$X_mem
  
  V_DOC    <- const$C_DOC * (const$V_water / 1000) * (1 / const$DOC_density) * 1000
  
  C_water  <- const$C_water_default
  C_free_water_initial <- C_water / (1 + sampled$K_DOC_water * const$C_DOC)
  C_DOC_water_initial  <- C_water - C_free_water_initial
  
  list(
    K_mem                = K_mem,
    K_tissue_water       = K_tissue_water,
    K_blood_water        = sampled$K_blood_water,
    P_gill               = P_gill,
    P_tissue             = P_tissue,
    V_DOC                = V_DOC,
    C_free_water_initial = C_free_water_initial,
    C_DOC_water_initial  = C_DOC_water_initial,
    k_elim_kidney        = sampled$k_elim_kidney,
    K_DOC_water          = sampled$K_DOC_water,
    C_water              = C_water
  )
}

# =========================================================
# 5) ODE models
# =========================================================
uptake_model <- function(t, y, parms) {
  const <- parms$const
  
  C_free_water   <- y[1]
  C_DOC_water    <- y[2]
  C_free_blood   <- y[3]
  C_bound_blood  <- y[4]
  C_tissue       <- y[5]
  
  P_gill         <- parms$P_gill
  P_tissue       <- parms$P_tissue
  K_tissue_water <- parms$K_tissue_water
  K_blood_water  <- parms$K_blood_water
  K_DOC_water    <- parms$K_DOC_water
  V_DOC          <- parms$V_DOC
  k_e            <- parms$k_elim_kidney
  
  t_days    <- t / 86400
  FW_t      <- const$FW + const$g * t_days
  V_blood   <- FW_t * const$Vf_blood
  V_tissue  <- FW_t * const$Vf_tissue
  
  growth_dilution_rate <- (const$g / FW_t) / 86400
  
  V_lipid         <- V_blood * const$Vf_lipid_blood
  V_phospholipids <- V_blood * const$Vf_phospholipids_blood
  V_protein       <- V_blood * const$Vf_protein_function_blood + V_blood * const$Vf_protein_structure_blood
  V_nonlipid      <- V_blood * const$Vf_nonlipid_blood
  V_colloids      <- V_lipid + V_phospholipids + V_protein + V_nonlipid
  
  M_water_blood   <- P_gill * const$A_gill * (C_free_water - C_free_blood)
  M_blood_tissue  <- P_tissue * const$A_tissue * (C_free_blood - C_tissue / K_tissue_water)
  M_renal_el      <- k_e * V_blood * const$Vf_water_blood * C_free_blood
  M_colloids_free <- const$k_des * V_colloids * (C_bound_blood - C_free_blood * K_blood_water)
  M_DOC_water     <- const$k_des * V_DOC * (C_DOC_water - C_free_water * K_DOC_water)
  
  dC_free_water   <- (-M_water_blood + M_DOC_water + M_renal_el) / const$V_water
  dC_DOC_water    <- (-M_DOC_water) / V_DOC
  dC_free_blood   <- (M_water_blood + M_colloids_free - M_blood_tissue - M_renal_el) / (V_blood * const$Vf_water_blood) - C_free_blood * growth_dilution_rate
  dC_bound_blood  <- (-M_colloids_free) / V_colloids - C_bound_blood * growth_dilution_rate
  dC_tissue       <- (M_blood_tissue) / V_tissue - C_tissue * growth_dilution_rate
  
  list(c(dC_free_water, dC_DOC_water, dC_free_blood, dC_bound_blood, dC_tissue))
}

# Natural uptake: water is held constant at initial values; state excludes water/DOC
uptake_natural_model <- function(t, y, parms) {
  const <- parms$const
  
  C_free_blood  <- y[1]
  C_bound_blood <- y[2]
  C_tissue      <- y[3]
  
  P_gill          <- parms$P_gill
  P_tissue        <- parms$P_tissue
  K_tissue_water  <- parms$K_tissue_water
  K_blood_water   <- parms$K_blood_water
  k_e             <- parms$k_elim_kidney
  
  # Pull fixed water concentrations from parms
  C_free_water <- parms$C_free_water_initial
  C_DOC_water  <- parms$C_DOC_water_initial
  
  t_days    <- t / 86400
  FW_t      <- const$FW + const$g * t_days
  V_blood   <- FW_t * const$Vf_blood
  V_tissue  <- FW_t * const$Vf_tissue
  
  growth_dilution_rate <- (const$g / FW_t) / 86400
  
  V_lipid         <- V_blood * const$Vf_lipid_blood
  V_phospholipids <- V_blood * const$Vf_phospholipids_blood
  V_protein       <- V_blood * const$Vf_protein_function_blood + V_blood * const$Vf_protein_structure_blood
  V_nonlipid      <- V_blood * const$Vf_nonlipid_blood
  V_colloids      <- V_lipid + V_phospholipids + V_protein + V_nonlipid
  
  M_water_blood   <- P_gill * const$A_gill * (C_free_water - C_free_blood)
  M_blood_tissue  <- P_tissue * const$A_tissue * (C_free_blood - C_tissue / K_tissue_water)
  M_renal_el      <- k_e * V_blood * const$Vf_water_blood * C_free_blood
  M_colloids_free <- const$k_des * V_colloids * (C_bound_blood - C_free_blood * K_blood_water)
  
  dC_free_blood  <- (M_water_blood + M_colloids_free - M_blood_tissue - M_renal_el) / (V_blood * const$Vf_water_blood) - C_free_blood * growth_dilution_rate
  dC_bound_blood <- (-M_colloids_free) / V_colloids - C_bound_blood * growth_dilution_rate
  dC_tissue      <- (M_blood_tissue) / V_tissue - C_tissue * growth_dilution_rate
  
  list(c(dC_free_blood, dC_bound_blood, dC_tissue))
}

depuration_model <- function(t, y, parms) {
  const <- parms$const
  
  C_free_water  <- y[1]
  C_DOC_water   <- y[2]
  C_free_blood  <- y[3]
  C_bound_blood <- y[4]
  C_tissue      <- y[5]
  
  P_gill          <- parms$P_gill
  P_tissue        <- parms$P_tissue
  K_tissue_water  <- parms$K_tissue_water
  K_blood_water   <- parms$K_blood_water
  K_DOC_water     <- parms$K_DOC_water
  V_DOC           <- parms$V_DOC
  k_e             <- parms$k_elim_kidney
  
  t_days_total  <- const$days_uptake + t / 86400
  FW_t          <- const$FW + const$g * t_days_total
  V_blood       <- FW_t * const$Vf_blood
  V_tissue      <- FW_t * const$Vf_tissue
  
  growth_dilution_rate <- (const$g / FW_t) / 86400
  
  V_lipid         <- V_blood * const$Vf_lipid_blood
  V_phospholipids <- V_blood * const$Vf_phospholipids_blood
  V_protein       <- V_blood * const$Vf_protein_function_blood + V_blood * const$Vf_protein_structure_blood
  V_nonlipid      <- V_blood * const$Vf_nonlipid_blood
  V_colloids      <- V_lipid + V_phospholipids + V_protein + V_nonlipid
  
  M_water_blood   <- P_gill * const$A_gill * (C_free_water - C_free_blood)
  M_blood_tissue  <- P_tissue * const$A_tissue * (C_free_blood - C_tissue / K_tissue_water)
  M_renal_el      <- k_e * V_blood * const$Vf_water_blood * C_free_blood
  M_DOC_water     <- const$k_des * V_DOC * (C_DOC_water - C_free_water * K_DOC_water)
  M_colloids_free <- const$k_des * V_colloids * (C_bound_blood - C_free_blood * K_blood_water)
  
  dC_free_water  <- (-M_water_blood + M_DOC_water + M_renal_el) / const$V_water
  dC_DOC_water   <- (-M_DOC_water) / V_DOC
  dC_free_blood  <- (M_water_blood + M_colloids_free - M_blood_tissue - M_renal_el) / (V_blood * const$Vf_water_blood) - C_free_blood * growth_dilution_rate
  dC_bound_blood <- (-M_colloids_free) / V_colloids - C_bound_blood * growth_dilution_rate
  dC_tissue      <- (M_blood_tissue) / V_tissue - C_tissue * growth_dilution_rate
  
  list(c(dC_free_water, dC_DOC_water, dC_free_blood, dC_bound_blood, dC_tissue))
}

# =========================================================
# 6) Total concentrations 
# =========================================================
calculate_total_conc <- function(out_data, phase, const) {
  
  if (phase == "uptake") {
    FW_t <- const$FW + const$g * (out_data[, "time"] / 86400)
    time_days <- out_data[, "time"] / 86400
  } else {
    FW_t <- const$FW + const$g * (const$days_uptake + out_data[, "time"] / 86400)
    time_days <- const$days_uptake + out_data[, "time"] / 86400
  }
  
  V_blood   <- FW_t * const$Vf_blood
  V_tissue  <- FW_t * const$Vf_tissue
  
  V_lipid         <- V_blood * const$Vf_lipid_blood
  V_phospholipids <- V_blood * const$Vf_phospholipids_blood
  V_protein       <- V_blood * const$Vf_protein_function_blood + V_blood * const$Vf_protein_structure_blood
  V_nonlipid      <- V_blood * const$Vf_nonlipid_blood
  V_colloids      <- V_lipid + V_phospholipids + V_protein + V_nonlipid
  
  C_blood <- (out_data[, "C_free_blood"] * V_blood * const$Vf_water_blood +
                out_data[, "C_bound_blood"] * V_colloids) / V_blood
  
  V_DOC_local   <- const$C_DOC * (const$V_water / 1000) * (1 / const$DOC_density) * 1000
  C_water_total <- out_data[, "C_free_water"] + out_data[, "C_DOC_water"] * V_DOC_local / const$V_water
  
  C_fish        <- (C_blood * V_blood + out_data[, "C_tissue"] * V_tissue) / (V_blood + V_tissue)
  
  data.frame(
    time_days = time_days,
    C_blood   = C_blood,
    C_tissue  = out_data[, "C_tissue"],
    C_fish    = C_fish,
    C_water   = C_water_total,
    phase     = phase
  )
}

# Total concentrations for NATURAL uptake output 
calculate_total_conc_natural_uptake <- function(out_data, const, derived) {
  FW_t <- const$FW + const$g * (out_data[, "time"] / 86400)
  time_days <- out_data[, "time"] / 86400
  
  V_blood  <- FW_t * const$Vf_blood
  V_tissue <- FW_t * const$Vf_tissue
  
  V_lipid         <- V_blood * const$Vf_lipid_blood
  V_phospholipids <- V_blood * const$Vf_phospholipids_blood
  V_protein       <- V_blood * const$Vf_protein_function_blood + V_blood * const$Vf_protein_structure_blood
  V_nonlipid      <- V_blood * const$Vf_nonlipid_blood
  V_colloids      <- V_lipid + V_phospholipids + V_protein + V_nonlipid
  
  C_blood <- (out_data[, "C_free_blood"] * V_blood * const$Vf_water_blood + out_data[, "C_bound_blood"] * V_colloids) / V_blood
  
  C_fish <- (C_blood * V_blood + out_data[, "C_tissue"] * V_tissue) / (V_blood + V_tissue)
  
  # In natural uptake, water total is constant at exposure C_water
  C_water_total <- rep(derived$C_water, length(time_days))
  
  data.frame(
    time_days = time_days,
    C_blood   = C_blood,
    C_tissue  = out_data[, "C_tissue"],
    C_fish    = C_fish,
    C_water   = C_water_total,
    phase     = "uptake"
  )
}

# =========================================================
# 7) Single PBTK run 
# =========================================================
run_single_pbtk <- function(derived, const, uptake_mode = c("lab", "natural")) {
  
  uptake_mode <- match.arg(uptake_mode)
  
  t_end_up  <- const$days_uptake * 86400
  t_end_dep <- const$days_depuration * 86400
  
  # Put derived + const into parms so ODEs can access derived fields too
  parms <- c(derived, list(const = const))
  
  # --- Uptake ---
  times_up <- seq(0, t_end_up, length.out = 1000)
  
  if (uptake_mode == "lab") {
    yini_up <- c(
      C_free_water  = derived$C_free_water_initial,
      C_DOC_water   = derived$C_DOC_water_initial,
      C_free_blood  = 0,
      C_bound_blood = 0,
      C_tissue      = 0
    )
    out_uptake <- tryCatch(
      ode(y = yini_up, times = times_up, func = uptake_model, parms = parms, method = "lsoda", hmax = 0),
      error = function(e) NULL
    )
    if (is.null(out_uptake)) return(NULL)
    
    results_uptake <- calculate_total_conc(out_uptake, "uptake", const)
    final_uptake <- tail(out_uptake, 1)
    
  } else {
    # NATURAL: water held constant; state excludes water/DOC
    yini_up <- c(
      C_free_blood  = 0,
      C_bound_blood = 0,
      C_tissue      = 0
    )
    out_uptake <- tryCatch(
      ode(y = yini_up, times = times_up, func = uptake_natural_model, parms = parms, method = "lsoda", hmax = 0),
      error = function(e) NULL
    )
    if (is.null(out_uptake)) return(NULL)
    
    # Rename columns to match downstream expectations
    colnames(out_uptake)[colnames(out_uptake) == "1"] <- "C_free_blood"
    colnames(out_uptake)[colnames(out_uptake) == "2"] <- "C_bound_blood"
    colnames(out_uptake)[colnames(out_uptake) == "3"] <- "C_tissue"
    
    results_uptake <- calculate_total_conc_natural_uptake(out_uptake, const, derived)
    
    # Construct "final_uptake" structure needed for depuration ICs
    final_uptake <- tail(out_uptake, 1)
  }
  
  # --- Depuration ---
  yini_dep <- c(
    C_free_water  = 0,
    C_DOC_water   = 0,
    C_free_blood  = final_uptake[, "C_free_blood"],
    C_bound_blood = final_uptake[, "C_bound_blood"],
    C_tissue      = final_uptake[, "C_tissue"]
  )
  times_dep <- seq(0, t_end_dep, length.out = 1000)
  
  out_depuration <- tryCatch(
    ode(y = yini_dep, times = times_dep, func = depuration_model, parms = parms, method = "lsoda", hmax = 0),
    error = function(e) NULL
  )
  if (is.null(out_depuration)) return(NULL)
  
  results_depuration <- calculate_total_conc(out_depuration, "depuration", const)
  
  # --- Metrics ---
  BAF <- results_uptake$C_fish[nrow(results_uptake)] / derived$C_water
  
  dep_data <- results_depuration %>%
    mutate(log_C_fish = log(C_fish + 1e-10))
  
  elim_model <- tryCatch(
    lm(log_C_fish ~ time_days, data = dep_data %>% filter(time_days > const$days_uptake)),
    error = function(e) NULL
  )
  if (is.null(elim_model)) return(NULL)
  
  k_d       <- -coef(elim_model)[2]
  half_life <- log(2) / k_d
  
  initial_slope <- (results_uptake$C_fish[2] - results_uptake$C_fish[1]) / (results_uptake$time_days[2] - results_uptake$time_days[1])
  k_u_initial <- (initial_slope * k_d) / derived$C_water
  
  nls_fit <- tryCatch(
    nls(
      C_fish ~ (k_u * derived$C_water / k_d) * (1 - exp(-k_d * time_days)),
      data = results_uptake,
      start = list(k_u = k_u_initial)
    ),
    error = function(e) NULL
  )
  
  if (is.null(nls_fit)) {
    BCF <- NA_real_
  } else {
    k_u <- coef(nls_fit)["k_u"]
    BCF <- as.numeric(k_u / k_d)
  }
  
  list(
    BAF = BAF,
    BCF = BCF,
    half_life = half_life,
    k_d = k_d,
    results_uptake = results_uptake,
    results_depuration = results_depuration
  )
}

# =========================================================
# 8) Monte Carlo 
# =========================================================
run_mc_simulation <- function(pfas_row, const, n_iter = 10000, verbose = TRUE, uptake_mode = c("lab","natural")) {
  
  uptake_mode <- match.arg(uptake_mode)
  
  BAF_values <- BCF_values <- half_life_values <- k_d_values <- rep(NA_real_, n_iter)
  all_uptake_results <- vector("list", n_iter)
  all_depuration_results <- vector("list", n_iter)
  successful_runs <- 0
  
  for (i in seq_len(n_iter)) {
    if (verbose && i %% 10 == 0) cat("  Iter", i, "/", n_iter, "\n")
    
    sampled <- sample_params(pfas_row, const)
    derived <- calculate_derived_params(sampled, pfas_row, const)
    result  <- run_single_pbtk(derived, const, uptake_mode = uptake_mode)
    
    if (!is.null(result)) {
      successful_runs <- successful_runs + 1
      BAF_values[i]       <- result$BAF
      BCF_values[i]       <- result$BCF
      half_life_values[i] <- result$half_life
      k_d_values[i]       <- result$k_d
      all_uptake_results[[i]]     <- result$results_uptake
      all_depuration_results[[i]] <- result$results_depuration
    }
  }
  
  list(
    BAF = BAF_values,
    BCF = BCF_values,
    half_life = half_life_values,
    k_d = k_d_values,
    uptake_results = all_uptake_results,
    depuration_results = all_depuration_results,
    n_successful = successful_runs
  )
}

# =========================================================
# 9) Summary table
# =========================================================
summarize_mc_results <- function(mc_results) {
  
  BAF_clean       <- na.omit(mc_results$BAF)
  BCF_clean       <- na.omit(mc_results$BCF)
  half_life_clean <- na.omit(mc_results$half_life)
  
  data.frame(
    Metric = c("BAF", "BCF", "log10(BAF)", "log10(BCF)", "Half-life (days)"),
    Mean = c(mean(BAF_clean), mean(BCF_clean),
             mean(log10(BAF_clean)), mean(log10(BCF_clean)),
             mean(half_life_clean)),
    SD = c(sd(BAF_clean), sd(BCF_clean),
           sd(log10(BAF_clean)), sd(log10(BCF_clean)),
           sd(half_life_clean)),
    CI_lower_2.5 = c(quantile(BAF_clean, 0.025), quantile(BCF_clean, 0.025),
                     quantile(log10(BAF_clean), 0.025), quantile(log10(BCF_clean), 0.025),
                     quantile(half_life_clean, 0.025)),
    CI_upper_97.5 = c(quantile(BAF_clean, 0.975), quantile(BCF_clean, 0.975),
                      quantile(log10(BAF_clean), 0.975), quantile(log10(BCF_clean), 0.975),
                      quantile(half_life_clean, 0.975)),
    Median = c(median(BAF_clean), median(BCF_clean),
               median(log10(BAF_clean)), median(log10(BCF_clean)),
               median(half_life_clean)),
    N = c(length(BAF_clean), length(BCF_clean),
          length(BAF_clean), length(BCF_clean),
          length(half_life_clean))
  )
}

# =========================================================
# 10) Time-course plot with uncertainty
# =========================================================
create_time_course_plot <- function(mc_results, const, chemical_name) {
  
  valid_uptake <- mc_results$uptake_results[!sapply(mc_results$uptake_results, is.null)]
  valid_dep    <- mc_results$depuration_results[!sapply(mc_results$depuration_results, is.null)]
  if (length(valid_uptake) == 0) return(NULL)
  
  time_uptake       <- valid_uptake[[1]]$time_days
  time_dep          <- valid_dep[[1]]$time_days
  
  # Blood compartment
  C_blood_uptake    <- sapply(valid_uptake, \(x) x$C_blood)
  C_blood_dep       <- sapply(valid_dep,    \(x) x$C_blood)
  C_water_uptake    <- sapply(valid_uptake, \(x) x$C_water)
  
  mean_C_blood_up   <- rowMeans(C_blood_uptake, na.rm = TRUE)
  ci_lower_up_blood <- apply(C_blood_uptake, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_up_blood <- apply(C_blood_uptake, 1, quantile, 0.975, na.rm = TRUE)
  
  mean_C_blood_dep  <- rowMeans(C_blood_dep, na.rm = TRUE)
  ci_lower_dep_blood<- apply(C_blood_dep, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_dep_blood<- apply(C_blood_dep, 1, quantile, 0.975, na.rm = TRUE)
  
  # Tissue compartment
  C_tissue_uptake    <- sapply(valid_uptake, \(x) x$C_tissue)
  C_tissue_dep       <- sapply(valid_dep,    \(x) x$C_tissue)
  C_water_uptake     <- sapply(valid_uptake, \(x) x$C_water)
  
  mean_C_tissue_up   <- rowMeans(C_tissue_uptake, na.rm = TRUE)
  ci_lower_up_tissue <- apply(C_tissue_uptake, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_up_tissue<- apply(C_tissue_uptake, 1, quantile, 0.975, na.rm = TRUE)
  
  mean_C_tissue_dep  <- rowMeans(C_tissue_dep, na.rm = TRUE)
  ci_lower_dep_tissue<- apply(C_tissue_dep, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_dep_tissue<- apply(C_tissue_dep, 1, quantile, 0.975, na.rm = TRUE)
   
  # Whole fish
  C_fish_uptake     <- sapply(valid_uptake, \(x) x$C_fish)
  C_fish_dep        <- sapply(valid_dep,    \(x) x$C_fish)
  C_water_uptake    <- sapply(valid_uptake, \(x) x$C_water)
  
  mean_C_fish_up    <- rowMeans(C_fish_uptake, na.rm = TRUE)
  ci_lower_up       <- apply(C_fish_uptake, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_up       <- apply(C_fish_uptake, 1, quantile, 0.975, na.rm = TRUE)
  
  mean_C_fish_dep   <- rowMeans(C_fish_dep, na.rm = TRUE)
  ci_lower_dep      <- apply(C_fish_dep, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper_dep      <- apply(C_fish_dep, 1, quantile, 0.975, na.rm = TRUE)
  
  # Water
  mean_C_water      <- rowMeans(C_water_uptake, na.rm = TRUE)
  ci_lower_up_water <- apply(C_water_uptake, 1, quantile, 0.025, na.rm = TRUE)
  ci_lower_up_water <- apply(C_water_uptake, 1, quantile, 0.975, na.rm = TRUE)
  
  plot_data_fish_blood <- rbind(
    data.frame(time = time_uptake, mean = mean_C_blood_up, ci_lower = ci_lower_up_blood, ci_upper = ci_upper_up_blood, phase = "uptake"),
    data.frame(time = time_dep,    mean = mean_C_blood_dep, ci_lower = ci_lower_dep_blood, ci_upper = ci_upper_dep_blood, phase = "depuration")
  )
  
  plot_data_fish_tissue <- rbind(
    data.frame(time = time_uptake, mean = mean_C_tissue_up, ci_lower = ci_lower_up_tissue, ci_upper = ci_upper_up_tissue, phase = "uptake"),
    data.frame(time = time_dep,    mean = mean_C_tissue_dep, ci_lower = ci_lower_dep_tissue, ci_upper = ci_upper_dep_tissue, phase = "depuration")
  )
  
  plot_data_fish <- rbind(
    data.frame(time = time_uptake, mean = mean_C_fish_up, ci_lower = ci_lower_up, ci_upper = ci_upper_up, phase = "uptake"),
    data.frame(time = time_dep,    mean = mean_C_fish_dep, ci_lower = ci_lower_dep, ci_upper = ci_upper_dep, phase = "depuration")
  )
  
  plot_data_water <- data.frame(time = time_uptake, mean = mean_C_water, ci_lower = ci_lower_up_water, ci_upper = ci_lower_up_water, phase = "uptake")
  
  p <- ggplot() +
    geom_ribbon(data = plot_data_fish, 
                aes(x = time, ymin = ci_lower * 1000, ymax = ci_upper * 1000),
                fill = "red", alpha = 0.3) +
    geom_line(data = plot_data_fish, 
              aes(x = time, y = mean * 1000), color = "red", size = 1.2) +
    geom_line(data = plot_data_water, 
              aes(x = time, y = mean * 1000), color = "blue", size = 1.2) +
    geom_vline(xintercept = const$days_uptake, linetype = "dashed", alpha = 0.5) +
    labs(x = "Time (days)",
         y = "Concentration (ng/mL)",
         title = paste("Fish PBTK Monte Carlo Simulation:", chemical_name),
         subtitle = "Red = Fish (mean ± 95% CI), Blue = Water") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)) +
    annotate("text", x = const$days_uptake/2, y = max(plot_data_fish$ci_upper) * 1100, 
             label = "Uptake", size = 4) +
    annotate("text", x = const$days_uptake + const$days_depuration/2, 
             y = max(plot_data_fish$ci_upper) * 1100, 
             label = "Depuration", size = 4)
  
  list(plot = p, data_fish_blood = plot_data_fish_blood, data_fish_tissue = plot_data_fish_tissue, data_fish = plot_data_fish, data_water = plot_data_water)
}

# =========================================================
# 11) Distribution plots
# =========================================================
create_distribution_plots <- function(mc_results) {
  
  dist_data <- data.frame(
    BAF = mc_results$BAF,
    BCF = mc_results$BCF,
    log10_BAF = log10(mc_results$BAF),
    log10_BCF = log10(mc_results$BCF),
    half_life = mc_results$half_life
    ) %>%
    filter(!is.na(BAF) & !is.na(BCF))
  
  # BAF distribution
  p_baf <- ggplot(dist_data, aes(x = log10_BAF)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    geom_vline(xintercept = mean(dist_data$log10_BAF, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    labs(x = "log10(BAF)", y = "Density", title = "Distribution of log10(BAF)") +
    theme_minimal()
  
  # BCF distribution
  p_bcf <- ggplot(dist_data, aes(x = log10_BCF)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    geom_vline(xintercept = mean(dist_data$log10_BCF, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    labs(x = "log10(BCF)", y = "Density", title = "Distribution of log10(BCF)") +
    theme_minimal()
  
  # Half-life distribution
  p_hl <- ggplot(dist_data, aes(x = half_life)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    geom_vline(xintercept = mean(dist_data$half_life, na.rm = TRUE), 
               color = "red", linetype = "dashed", size = 1) +
    labs(x = "Half-life (days)", y = "Density", title = "Distribution of Half-life") +
    theme_minimal()
  
  list(BAF = p_baf, BCF = p_bcf, half_life = p_hl, dist_data = dist_data)
}

# =========================================================
# 12) Full workflow: loop over PFAS, save outputs
# =========================================================
run_pfas_workflow <- function(excel_path, sheet_name, out_dir, n_iter, verbose, uptake_mode = c("lab","natural"), target = NULL) {
  
  uptake_mode <- match.arg(uptake_mode)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  const <- build_constants()
  pfas_df <- read_pfas_table(excel_path, sheet_name)
  
  # Pick specific PFAS if target is provided
  if (!is.null(target)) {
    pfas_df <- pfas_df %>% filter(chemical_name == target)
  }
  
  summary_all <- list()
  
  for (r in seq_len(nrow(pfas_df))) {
    pfas_row <- pfas_df[r, , drop = FALSE]
    chemical_name <- as.character(pfas_row$chemical_name)
    safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", chemical_name)
    
    cat("\n=========================================\n")
    cat("Running PFAS:", chemical_name, " (row ", r, "/", nrow(pfas_df), ")\n", sep = "")
    cat("Uptake mode:", uptake_mode, "\n")
    cat("=========================================\n")
    
    mc <- run_mc_simulation(pfas_row, const, n_iter = n_iter, verbose = verbose, uptake_mode = uptake_mode)
    summ <- summarize_mc_results(mc)
    
    tc <- create_time_course_plot(mc, const, chemical_name)
    dp <- create_distribution_plots(mc)
    
    raw_results <- data.frame(
      BAF = mc$BAF,
      BCF = mc$BCF,
      log10_BAF = log10(mc$BAF),
      log10_BCF = log10(mc$BCF),
      half_life = mc$half_life,
      k_d = mc$k_d
    )
    
    time_course_fish_blood_df <- if (!is.null(tc)) tc$data_fish_blood else data.frame()
    time_course_fish_tissue_df <- if (!is.null(tc)) tc$data_fish_tissue else data.frame()
    time_course_fish_df <- if (!is.null(tc)) tc$data_fish else data.frame()
    time_course_water_df <- if (!is.null(tc)) tc$data_water else data.frame()
    
    out_xlsx <- file.path(out_dir, paste0("Fish_PBTK_MC_Results_", safe_name, ".xlsx"))
    
    write.xlsx(
      list(
        PFAS_Input_Row = pfas_row,
        Summary = summ,
        Raw_Results = raw_results,
        Time_Course_Fish_Blood = time_course_fish_blood_df,
        Time_Course_Fish_Tissue = time_course_fish_tissue_df,
        Time_Course_Fish = time_course_fish_df,
        Time_Course_Water = time_course_water_df,
        Dist_Data = dp$dist_data
      ),
      file = out_xlsx,
      rowNames = FALSE
    )
    
    if (!is.null(tc)) {
      ggsave(file.path(out_dir, paste0("TimeCourse_", safe_name, ".png")),
             tc$plot, width = 8, height = 5, dpi = 300)
    }
    ggsave(file.path(out_dir, paste0("Dist_log10BAF_", safe_name, ".png")),
           dp$BAF, width = 6, height = 4, dpi = 300)
    ggsave(file.path(out_dir, paste0("Dist_log10BCF_", safe_name, ".png")),
           dp$BCF, width = 6, height = 4, dpi = 300)
    ggsave(file.path(out_dir, paste0("Dist_halflife_", safe_name, ".png")),
           dp$half_life, width = 6, height = 4, dpi = 300)
    
    summary_all[[chemical_name]] <- summ %>% mutate(chemical_name = chemical_name)
    
    cat("Saved: ", out_xlsx, "\n", sep = "")
    cat("Successful runs: ", mc$n_successful, " / ", n_iter, "\n", sep = "")
  }
  
  summary_all_df <- bind_rows(summary_all)
  write.xlsx(summary_all_df, file = file.path(out_dir, paste0("Summary_ALL_PFAS_", uptake_mode, ".xlsx")), rowNames = FALSE)
  
  invisible(list(summary_all = summary_all_df, constants = const, pfas_table = pfas_df))
}

# =========================================================
# RUN for lab
# =========================================================
tic()

results <- run_pfas_workflow(
  excel_path = excel_path,
  sheet_name = sheet_name,
  out_dir    = "PBTK_lab_outputs",
  n_iter     = n_iter,
  verbose    = verbose,
  uptake_mode = "lab"
)

toc()


# =========================================================
# RUN for nat
# =========================================================
tic()

results <- run_pfas_workflow(
  excel_path = excel_path,
  sheet_name = sheet_name,
  out_dir    = "PBTK_nat_outputs",
  n_iter     = n_iter,
  verbose    = verbose,
  uptake_mode = "natural"
)

toc()

# =========================================================
# RUN for ABL = 0.016 cm vs 0.082 cm
# =========================================================
# X_ABL = 0.082

tic()

results <- run_pfas_workflow(
  excel_path = excel_path,
  sheet_name = sheet_name,
  out_dir    = "PBTK_ABL_outputs",
  n_iter     = n_iter,
  verbose    = verbose,
  uptake_mode = "lab"
)

toc()


# =========================================================
# RUN for C_D0C = 0.0001 g mL-1 vs. 0.000001 g mL-1
# =========================================================
# C_DOC = 1e-4

tic()

results <- run_pfas_workflow(
  excel_path = excel_path,
  sheet_name = sheet_name,
  out_dir    = "PBTK_DOC_outputs",
  n_iter     = n_iter,
  verbose    = verbose,
  uptake_mode = "lab"
)

toc()


tic()

results <- run_pfas_workflow(
  excel_path = excel_path,
  sheet_name = sheet_name,
  out_dir    = "PBTK_DOC_nat_outputs",
  n_iter     = n_iter,
  verbose    = verbose,
  uptake_mode = "nat"
)

toc()
