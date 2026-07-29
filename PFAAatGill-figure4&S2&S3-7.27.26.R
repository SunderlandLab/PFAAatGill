##########################################################################################
# Last edit: 2026-07-11
##########################################################################################
rm(list = ls())

library(tidyverse) 
library(readxl)   
library(ggpubr)
library(RColorBrewer) 

setwd("/Users/huzhiji/OneDrive - Harvard University/Research/PFAAatGill/")

# =========================================================
## Load properties
# =========================================================
PFAS_log <- read_excel('PFAA_inputs.xlsx', sheet = "Log", col_names = TRUE) %>% column_to_rownames('...1')

# =========================================================
## Main workflow
# =========================================================
distribution_workflow = function(target){
  
  # partitioning coefficient (used in calculation)
  # names(PFAS_log) = c("CF4 PFSA" ,"CF6 PFSA", "CF8 PFSA", "CF5 PFCA", "CF6 PFCA", "CF7 PFCA", "CF8 PFCA", "CF9 PFCA", "CF10 PFCA")
  
  PFAS = data.frame(lapply(PFAS_log, function(x) 10^(x)))
  rownames(PFAS) = c("KOW", "KMW", "KFPW", "KSPW", "KNW")
  
  beta_NLOM = 0.035
  
  # fish tissue volume fraction (default suggested by Amritage)
  Vf_neutral_lipids = 0.04                      #L neutral lipids / L fish
  Vf_phospholipids = 0.01                      #L phospholipids / L fish
  Vf_protein = 0.003                           #L proteins / L fish
  Vf_nonlipid_organic_matter = 0.15            #L nonlipid organic matter / L fish
  Vf_water = 1 - Vf_neutral_lipids - Vf_phospholipids - Vf_protein - Vf_nonlipid_organic_matter  #L water / L fish
  
  # fish organ volume fraction 
  Vf_gill_ori = 0.02385                   
  Vf_liver_ori = 0.0114                   
  Vf_kidney_ori = 0.0084
  Vf_muscle_ori = 0.465
  Vf_adipose_ori = 0.0089
  Vf_blood_ori = 0.365
  Vf_rest_ori = 1- Vf_gill_ori - Vf_liver_ori - Vf_kidney_ori - Vf_muscle_ori - Vf_adipose_ori - Vf_blood_ori
  
  Vf_total_withoutblood = 1 - Vf_blood_ori
  
  Vf_gill = Vf_gill_ori / Vf_total_withoutblood                
  Vf_liver = Vf_liver_ori / Vf_total_withoutblood          
  Vf_kidney = Vf_kidney_ori / Vf_total_withoutblood 
  Vf_muscle = Vf_muscle_ori / Vf_total_withoutblood
  Vf_adipose = Vf_adipose_ori / Vf_total_withoutblood
  Vf_rest = Vf_rest_ori / Vf_total_withoutblood
  
  if (target == "gill without blood"){
    organ = data.frame(Vf_gill, Vf_liver, Vf_kidney, Vf_muscle, Vf_adipose, Vf_rest)
  } else if (target == "gill with blood") {
    organ = data.frame(Vf_gill_ori, Vf_blood_ori, Vf_liver_ori, Vf_kidney_ori, Vf_muscle_ori, Vf_adipose_ori, Vf_rest_ori)
    names(organ) = c("Vf_gill", "Vf_blood","Vf_liver", "Vf_kidney", "Vf_muscle", "Vf_adipose", "Vf_rest")
  } else if (target == "within gill") {
    organ = NULL }
  
  # fish gill
  M_gill_lipid = 0.523 # mg/100 mg wt
  M_gill_protein = 20.4 # mg/100 mg wt
  M_gill_glycogen = 0.516 # mg/100 mg wt
  
  percent_gill_phospholipids = 0.4361
  percent_gill_albumin = 0.2
  
  Vf_gill_phospholipids = (M_gill_lipid * percent_gill_phospholipids) / 100
  Vf_gill_neutral_lipids = (M_gill_lipid - M_gill_lipid * percent_gill_phospholipids) /100
  
  Vf_gill_functional_proteins = (M_gill_protein * percent_gill_albumin) / 100
  Vf_gill_structural_proteins = (M_gill_protein - M_gill_protein * percent_gill_albumin) / 100
  
  Vf_gill_nonlipid_organic_matter = M_gill_glycogen / 100
  
  Vf_gill_water = 1 - Vf_gill_phospholipids - Vf_gill_neutral_lipids - 
    Vf_gill_functional_proteins - Vf_gill_structural_proteins - 
    Vf_gill_nonlipid_organic_matter
  
  gill = data.frame(Vf_gill_neutral_lipids, Vf_gill_phospholipids, 
                    Vf_gill_functional_proteins, Vf_gill_structural_proteins,
                    Vf_gill_nonlipid_organic_matter, Vf_gill_water)
  
  # fish blood
  M_blood_lipid = 1.84 # mg/100 mg wt
  M_blood_protein = 4.12 # mg/100 mg wt
  
  M_blood_glycogen = 0.433 # mg/100 mg wt
  
  percent_blood_phospholipids = 0.5
  percent_blood_albumin = 0.7088607
  
  Vf_blood_neutral_lipids = (M_blood_lipid - M_blood_lipid * percent_blood_phospholipids) / 100
  Vf_blood_phospholipids = (M_blood_lipid * percent_blood_phospholipids) / 100
  
  Vf_blood_functional_proteins = (M_blood_protein * percent_blood_albumin) / 100
  Vf_blood_structural_proteins = (M_blood_protein - M_blood_protein * percent_blood_albumin) / 100
  
  Vf_blood_nonlipid_organic_matter = M_blood_glycogen / 100
  
  Vf_blood_water = 1 - Vf_blood_phospholipids - Vf_blood_neutral_lipids - 
    Vf_blood_functional_proteins - Vf_blood_structural_proteins - 
    Vf_blood_nonlipid_organic_matter
  
  blood = data.frame(Vf_blood_neutral_lipids, Vf_blood_phospholipids, 
                     Vf_blood_functional_proteins, Vf_blood_structural_proteins,
                     Vf_blood_nonlipid_organic_matter, Vf_blood_water)
  
  # fish liver
  M_liver_lipid = 0.753 # mg/100 mg wt
  M_liver_protein = 12.33 # mg/100 mg wt
  M_liver_glycogen = 3.5 # mg/100 mg wt
  
  percent_liver_phospholipids = 24.72 / (28.02 + 24.72)
  percent_liver_albumin = 0.313 / 1.17
  
  Vf_liver_phospholipids = (M_liver_lipid * percent_liver_phospholipids) / 100
  Vf_liver_neutral_lipids = (M_liver_lipid - M_liver_lipid * percent_liver_phospholipids) /100
  
  Vf_liver_functional_proteins = (M_liver_protein * percent_liver_albumin) / 100
  Vf_liver_structural_proteins = (M_liver_protein - M_liver_protein * percent_liver_albumin) / 100
  
  Vf_liver_nonlipid_organic_matter = M_liver_glycogen / 100
  
  Vf_liver_water = 1 - Vf_liver_phospholipids - Vf_liver_neutral_lipids - 
    Vf_liver_functional_proteins - Vf_liver_structural_proteins - 
    Vf_liver_nonlipid_organic_matter
  
  liver = data.frame(Vf_liver_neutral_lipids, Vf_liver_phospholipids, 
                     Vf_liver_functional_proteins, Vf_liver_structural_proteins,
                     Vf_liver_nonlipid_organic_matter, Vf_liver_water)
  
  # fish kidney
  M_kidney_lipid = 0.66 # mg/100 mg wt
  M_kidney_protein = 18.53 # mg/100 mg wt
  M_kidney_glycogen = 0.18 # mg/100 mg wt
  
  percent_kidney_phospholipids = 0.26
  percent_kidney_albumin = 0.313 / 1.17
  
  Vf_kidney_phospholipids = (M_kidney_lipid * percent_kidney_phospholipids) / 100
  Vf_kidney_neutral_lipids = (M_kidney_lipid - M_kidney_lipid * percent_kidney_phospholipids) /100
  
  Vf_kidney_functional_proteins = (M_kidney_protein * percent_kidney_albumin) / 100
  Vf_kidney_structural_proteins = (M_kidney_protein - M_kidney_protein * percent_kidney_albumin) / 100
  
  Vf_kidney_nonlipid_organic_matter = M_kidney_glycogen / 100
  
  Vf_kidney_water = 1 - Vf_kidney_phospholipids - Vf_kidney_neutral_lipids - 
    Vf_kidney_functional_proteins - Vf_kidney_structural_proteins - 
    Vf_kidney_nonlipid_organic_matter
  
  kidney = data.frame(Vf_kidney_neutral_lipids, Vf_kidney_phospholipids, 
                      Vf_kidney_functional_proteins, Vf_kidney_structural_proteins,
                      Vf_kidney_nonlipid_organic_matter, Vf_kidney_water)
  
  # Korgan/water
  BAF_organ <- function(organ, PFAS, beta_NLOM){
    result_all = data.frame()
    for (i in 1:ncol(PFAS)){
      Exact_PFAS = PFAS[i]
      
      # Tissue parameter
      Vf_neutral_lipids = organ[grepl("neutral_lipids", names(organ))]
      Vf_phospholipids = organ[grepl("phospholipids", names(organ))]
      Vf_functional_proteins = organ[grepl("functional_proteins", names(organ))]
      Vf_structural_proteins = organ[grepl("structural_proteins", names(organ))]
      Vf_nonlipid_organic_matter = organ[grepl("nonlipid_organic_matter", names(organ))]
      Vf_water = organ[grepl("water", names(organ))]
      
      # PFAS parameter
      KOW =  Exact_PFAS['KOW',]
      KMW = Exact_PFAS['KMW',]
      KFPW = Exact_PFAS['KFPW',]
      KSPW = Exact_PFAS['KSPW',]
      KNW = Exact_PFAS['KNW',]
      
      # calculation
      BAF = Vf_neutral_lipids * KOW + Vf_phospholipids * KMW + 
        Vf_functional_proteins *  KFPW + Vf_structural_proteins * KSPW + 
        Vf_nonlipid_organic_matter * KNW * beta_NLOM + Vf_water
      
      species = colnames(Exact_PFAS)
      logBAF = log10(BAF)
      Organ = deparse(substitute(organ))
      result_all = rbind(result_all, cbind(species, BAF, logBAF, Organ))
    }
    result = result_all
    names(result) = c("species", "BAF", "logBAF", "Organ")
    return(result)
  }
  
  gill.output = BAF_organ(gill, PFAS, beta_NLOM)
  liver.output = BAF_organ(liver, PFAS, beta_NLOM)
  kidney.output = BAF_organ(kidney, PFAS, beta_NLOM)
  blood.output = BAF_organ(blood, PFAS, beta_NLOM)
  muscle.output = data.frame(t(PFAS[4,]))
  names(muscle.output) = "BAF"
  
  adipose.output = data.frame(t(PFAS[1,]))
  names(adipose.output) = "BAF"
  
  if (target == "gill without blood") {
    organ.output = cbind(gill.output[c('species')], gill.output[c('BAF')], liver.output[c('BAF')], 
                         kidney.output[c('BAF')], muscle.output[c('BAF')], adipose.output[c('BAF')])
    names(organ.output) = c('species', 'K_gill_water', 'K_liver_water', 
                            'K_kidney_water', 'K_muscle_water', 'K_adipose_water')
  } else if (target == "gill with blood") {
    organ.output = cbind(gill.output[c('species')], gill.output[c('BAF')], blood.output[c('BAF')], liver.output[c('BAF')], 
                         kidney.output[c('BAF')], muscle.output[c('BAF')], adipose.output[c('BAF')])
    names(organ.output) = c('species', 'K_gill_water', 'K_blood_water', 'K_liver_water', 
                            'K_kidney_water', 'K_muscle_water', 'K_adipose_water')
  } else if (target == "within gill") {
    organ.output = NULL
  }
  
  fraction_organ_noblood <- function(organ, organ.partition){
    result_all = data.frame()
    for (i in 1:nrow(organ.partition)){
      Exact_PFAS = organ.partition[i,]
      
      # organ parameter
      Vf_gill = Vf_gill                  
      Vf_liver = Vf_liver                
      Vf_kidney = Vf_kidney
      Vf_muscle = Vf_muscle
      Vf_adipose = Vf_adipose
      Vf_water =  Vf_rest
      
      # PFAS parameter
      K_gill_water =  Exact_PFAS[,'K_gill_water']
      K_liver_water = Exact_PFAS[,'K_liver_water']
      K_kidney_water = Exact_PFAS[,'K_kidney_water']
      K_muscle_water = Exact_PFAS[,'K_muscle_water']
      K_adipose_water = Exact_PFAS[,'K_adipose_water']
      
      # calculation
      f_gill = 1 / (1 + (K_liver_water/K_gill_water) * (Vf_liver/Vf_gill) +
                      (K_kidney_water/K_gill_water) * (Vf_kidney/Vf_gill) +
                      (K_muscle_water/K_gill_water) * (Vf_muscle/Vf_gill) + 
                      (K_adipose_water/K_gill_water) * (Vf_adipose/Vf_gill) +
                      (1/K_gill_water) * (Vf_water/Vf_gill))
      
      f_liver = 1 / (1 + (K_gill_water/K_liver_water) * (Vf_gill/Vf_liver) +
                       (K_kidney_water/K_liver_water) * (Vf_kidney/Vf_liver) +
                       (K_muscle_water/K_liver_water) * (Vf_muscle/Vf_liver) + 
                       (K_adipose_water/K_liver_water) * (Vf_adipose/Vf_liver) +
                       (1/K_liver_water) * (Vf_water/Vf_liver))
      
      f_kidney = 1 / (1 + (K_gill_water/K_kidney_water) * (Vf_gill/Vf_kidney) +
                        (K_liver_water/K_kidney_water) * (Vf_liver/Vf_kidney) +
                        (K_muscle_water/K_kidney_water) * (Vf_muscle/Vf_kidney) + 
                        (K_adipose_water/K_kidney_water) * (Vf_adipose/Vf_kidney) +
                        (1/K_kidney_water) * (Vf_water/Vf_kidney))
      
      f_muscle = 1 / (1 + (K_gill_water/K_muscle_water) * (Vf_gill/Vf_muscle) +
                        (K_liver_water/K_muscle_water) * (Vf_liver/Vf_muscle) +
                        (K_kidney_water/K_muscle_water) * (Vf_kidney/Vf_muscle) + 
                        (K_adipose_water/K_muscle_water) * (Vf_adipose/Vf_muscle) +
                        (1/K_muscle_water) * (Vf_water/Vf_muscle))
      
      f_adipose = 1 / (1 + (K_gill_water/K_adipose_water) * (Vf_gill/Vf_adipose) +
                         (K_liver_water/K_adipose_water) * (Vf_liver/Vf_adipose) +
                         (K_kidney_water/K_adipose_water) * (Vf_kidney/Vf_adipose) + 
                         (K_muscle_water/K_adipose_water) * (Vf_muscle/Vf_adipose) +
                         (1/K_adipose_water) * (Vf_water/Vf_adipose))
      
      f_rest = 1 / (1 + K_gill_water * (Vf_gill/Vf_water) +
                      K_liver_water * (Vf_liver/Vf_water) +
                      K_kidney_water * (Vf_kidney/Vf_water) + 
                      K_muscle_water * (Vf_muscle/Vf_water) +
                      K_adipose_water * (Vf_adipose/Vf_water))
      
      species = rownames(Exact_PFAS)
      result_all = rbind(result_all, cbind(species, f_gill,
                                           f_liver, f_kidney, 
                                           f_muscle, f_adipose, f_rest))
    }
    result = result_all
    return(result)
  }
  
  fraction_organ <- function(organ, organ.partition){
    result_all = data.frame()
    for (i in 1:nrow(organ.partition)){
      Exact_PFAS = organ.partition[i,]
      
      # organ parameter
      Vf_gill = Vf_gill_ori                  
      Vf_blood = Vf_blood_ori                  
      Vf_liver = Vf_liver_ori                
      Vf_kidney = Vf_kidney_ori
      Vf_muscle = Vf_muscle_ori
      Vf_adipose = Vf_adipose_ori
      Vf_water =  Vf_rest_ori
      
      # PFAS parameter
      K_gill_water =  Exact_PFAS[,'K_gill_water']
      K_blood_water = Exact_PFAS[,'K_blood_water']
      K_liver_water = Exact_PFAS[,'K_liver_water']
      K_kidney_water = Exact_PFAS[,'K_kidney_water']
      K_muscle_water = Exact_PFAS[,'K_muscle_water']
      K_adipose_water = Exact_PFAS[,'K_adipose_water']
      
      # calculation
      f_gill = 1 / (1 + (K_blood_water/K_gill_water) * (Vf_blood/Vf_gill) +
                      (K_liver_water/K_gill_water) * (Vf_liver/Vf_gill) +
                      (K_kidney_water/K_gill_water) * (Vf_kidney/Vf_gill) +
                      (K_muscle_water/K_gill_water) * (Vf_muscle/Vf_gill) + 
                      (K_adipose_water/K_gill_water) * (Vf_adipose/Vf_gill) +
                      (1/K_gill_water) * (Vf_water/Vf_gill))
      
      
      f_blood = 1 / (1 + (K_gill_water/K_blood_water) * (Vf_gill/Vf_blood) +
                       (K_liver_water/K_blood_water) * (Vf_liver/Vf_blood) +
                       (K_kidney_water/K_blood_water) * (Vf_kidney/Vf_blood) +
                       (K_muscle_water/K_blood_water) * (Vf_muscle/Vf_blood) + 
                       (K_adipose_water/K_blood_water) * (Vf_adipose/Vf_blood) +
                       (1/K_blood_water) * (Vf_water/Vf_blood))
      
      
      f_liver = 1 / (1 + (K_gill_water/K_liver_water) * (Vf_gill/Vf_liver) +
                       (K_blood_water/K_liver_water) * (Vf_blood/Vf_liver) +
                       (K_kidney_water/K_liver_water) * (Vf_kidney/Vf_liver) +
                       (K_muscle_water/K_liver_water) * (Vf_muscle/Vf_liver) + 
                       (K_adipose_water/K_liver_water) * (Vf_adipose/Vf_liver) +
                       (1/K_liver_water) * (Vf_water/Vf_liver))
      
      f_kidney = 1 / (1 + (K_gill_water/K_kidney_water) * (Vf_gill/Vf_kidney) +
                        (K_blood_water/K_kidney_water) * (Vf_blood/Vf_kidney) +
                        (K_liver_water/K_kidney_water) * (Vf_liver/Vf_kidney) +
                        (K_muscle_water/K_kidney_water) * (Vf_muscle/Vf_kidney) + 
                        (K_adipose_water/K_kidney_water) * (Vf_adipose/Vf_kidney) +
                        (1/K_kidney_water) * (Vf_water/Vf_kidney))
      
      f_muscle = 1 / (1 + (K_gill_water/K_muscle_water) * (Vf_gill/Vf_muscle) +
                        (K_blood_water/K_muscle_water) * (Vf_blood/Vf_muscle) +
                        (K_liver_water/K_muscle_water) * (Vf_liver/Vf_muscle) +
                        (K_kidney_water/K_muscle_water) * (Vf_kidney/Vf_muscle) + 
                        (K_adipose_water/K_muscle_water) * (Vf_adipose/Vf_muscle) +
                        (1/K_muscle_water) * (Vf_water/Vf_muscle))
      
      f_adipose = 1 / (1 + (K_gill_water/K_adipose_water) * (Vf_gill/Vf_adipose) +
                         (K_blood_water/K_adipose_water) * (Vf_blood/Vf_adipose) +
                         (K_liver_water/K_adipose_water) * (Vf_liver/Vf_adipose) +
                         (K_kidney_water/K_adipose_water) * (Vf_kidney/Vf_adipose) + 
                         (K_muscle_water/K_adipose_water) * (Vf_muscle/Vf_adipose) +
                         (1/K_adipose_water) * (Vf_water/Vf_adipose))
      
      f_rest = 1 / (1 + K_gill_water * (Vf_gill/Vf_water) +
                      K_blood_water * (Vf_blood/Vf_water) +
                      K_liver_water * (Vf_liver/Vf_water) +
                      K_kidney_water * (Vf_kidney/Vf_water) + 
                      K_muscle_water * (Vf_muscle/Vf_water) +
                      K_adipose_water * (Vf_adipose/Vf_water))
      
      species = rownames(Exact_PFAS)
      result_all = rbind(result_all, cbind(species, f_gill, f_blood, 
                                           f_liver, f_kidney, 
                                           f_muscle, f_adipose, f_rest))
    }
    result = result_all
    return(result)
  }
  
  fraction_gill <- function(gill, PFAS){
    result_all = data.frame()
    for (i in 1:ncol(PFAS)){
      Exact_PFAS = PFAS[i]
      
      Vf_neutral_lipids = gill['Vf_gill_neutral_lipids'] %>% as.numeric()
      Vf_phospholipids = gill['Vf_gill_phospholipids'] %>% as.numeric()
      Vf_functional_proteins = gill['Vf_gill_functional_proteins'] %>% as.numeric()
      Vf_structural_proteins = gill['Vf_gill_structural_proteins'] %>% as.numeric()
      Vf_nonlipid_organic_matter = gill['Vf_gill_nonlipid_organic_matter'] %>% as.numeric()
      Vf_water = gill['Vf_gill_water'] %>% as.numeric()
      
      # parameter
      KOW =  Exact_PFAS['KOW',]
      KMW = Exact_PFAS['KMW',]
      KFPW = Exact_PFAS['KFPW',]
      KSPW = Exact_PFAS['KSPW',]
      KNW = Exact_PFAS['KNW',]
      
      # calculation
      #BAF = Vf_neutral_lipids * KOW + Vf_phospholipids * KMW + 
      #Vf_functional_proteins *  KFPW + Vf_structural_proteins * KSPW + 
      #Vf_nonlipid_organic_matter * KNW * beta_NLOM + Vf_water
      
      f_neutral_lipids = 1 / (1 + (KMW/KOW) * (Vf_phospholipids/Vf_neutral_lipids) +
                                (KFPW/KOW) * (Vf_functional_proteins/Vf_neutral_lipids) +
                                (KSPW/KOW) * (Vf_structural_proteins/Vf_neutral_lipids) +
                                ((KNW*beta_NLOM)/KOW) * (Vf_nonlipid_organic_matter/Vf_neutral_lipids) +
                                (1/KOW) * (Vf_water/Vf_neutral_lipids))
      
      f_phospholipids = 1 / (1 + (KOW/KMW) * (Vf_neutral_lipids/Vf_phospholipids) +
                               (KFPW/KMW) * (Vf_functional_proteins/Vf_phospholipids) +
                               (KSPW/KMW) * (Vf_structural_proteins/Vf_phospholipids) +
                               ((KNW*beta_NLOM)/KMW) * (Vf_nonlipid_organic_matter/Vf_phospholipids) +
                               (1/KMW) * (Vf_water/Vf_phospholipids))
      
      f_functional_proteins = 1 / (1 + (KOW/KFPW) * (Vf_neutral_lipids/Vf_functional_proteins) +
                                     (KMW/KFPW) * (Vf_phospholipids/Vf_functional_proteins) +
                                     (KSPW/KFPW) * (Vf_structural_proteins/Vf_functional_proteins) +
                                     ((KNW*beta_NLOM)/KFPW) * (Vf_nonlipid_organic_matter/Vf_functional_proteins) +
                                     (1/KFPW) * (Vf_water/Vf_functional_proteins))
      
      f_structural_proteins = 1 / (1 + (KOW/KSPW) * (Vf_neutral_lipids/Vf_structural_proteins) +
                                     (KMW/KSPW) * (Vf_phospholipids/Vf_structural_proteins) +
                                     (KFPW/KSPW) * (Vf_functional_proteins/Vf_structural_proteins) +
                                     ((KNW*beta_NLOM)/KSPW) * (Vf_nonlipid_organic_matter/Vf_structural_proteins) +
                                     (1/KSPW) * (Vf_water/Vf_structural_proteins))
      
      f_nonlipid_organic_matter = 1 / (1 + (KOW/(KNW*beta_NLOM)) * (Vf_neutral_lipids/Vf_nonlipid_organic_matter) +
                                         (KMW/(KNW*beta_NLOM)) * (Vf_phospholipids/Vf_nonlipid_organic_matter) +
                                         (KFPW/(KNW*beta_NLOM)) * (Vf_functional_proteins/Vf_nonlipid_organic_matter) +
                                         (KSPW/(KNW*beta_NLOM)) * (Vf_structural_proteins/Vf_nonlipid_organic_matter) +
                                         (1/(KNW*beta_NLOM)) * (Vf_water/Vf_nonlipid_organic_matter))
      
      f_water = 1 / (1 + KOW * (Vf_neutral_lipids/Vf_water) +
                       KMW * (Vf_phospholipids/Vf_water) +
                       KFPW * (Vf_functional_proteins/Vf_water) +
                       KSPW * (Vf_structural_proteins/Vf_water) +
                       (KNW*beta_NLOM) * (Vf_nonlipid_organic_matter/Vf_water))
      
      species = colnames(Exact_PFAS)
      #logBAF = log10(BAF)
      result_all = rbind(result_all, cbind(species, #logBAF, 
                                           f_neutral_lipids, f_phospholipids, 
                                           f_functional_proteins, f_structural_proteins, 
                                           f_nonlipid_organic_matter, f_water))
    }
    result = result_all
    return(result)
  }
  
  if (target == "gill without blood") {
    # fraction
    fraction = fraction_organ_noblood(organ, organ.output)
    
    # concentration
    concentration = data.frame(fraction[,'species'],
                               as.numeric(fraction[,'f_gill']) / organ[,'Vf_gill'],
                               as.numeric(fraction[,'f_liver']) / organ[,'Vf_liver'],
                               as.numeric(fraction[,'f_kidney']) / organ[,'Vf_kidney'],
                               as.numeric(fraction[,'f_muscle']) / organ[,'Vf_muscle'],
                               as.numeric(fraction[,'f_adipose']) / organ[,'Vf_adipose'],
                               as.numeric(fraction[,'f_rest']) / organ[,'Vf_rest'])
    names(concentration) = c('species', 'C_gill', 'C_liver', 'C_kidney', 'C_muscle', 'C_adipose', 'C_rest')
    
    result = list(fraction = fraction, concentration = concentration)
    
  } else if (target == "gill with blood") {
    # fraction
    fraction = fraction_organ(organ, organ.output)
    
    # concentration
    concentration = data.frame(fraction[,'species'],
                               as.numeric(fraction[,'f_gill']) / organ[,'Vf_gill'],
                               as.numeric(fraction[,'f_blood']) / organ[,'Vf_blood'],
                               as.numeric(fraction[,'f_liver']) / organ[,'Vf_liver'],
                               as.numeric(fraction[,'f_kidney']) / organ[,'Vf_kidney'],
                               as.numeric(fraction[,'f_muscle']) / organ[,'Vf_muscle'],
                               as.numeric(fraction[,'f_adipose']) / organ[,'Vf_adipose'],
                               as.numeric(fraction[,'f_rest']) / organ[,'Vf_rest'])
    names(concentration) = c('species', 'C_gill', 'C_blood', 'C_liver',
                             'C_kidney', 'C_muscle', 'C_adipose', 'C_rest')
    result = list(fraction = fraction, concentration = concentration)
    
  } else if (target == "within gill") {
    # fraction
    fraction = fraction_gill(gill, PFAS)
    
    # concentration
    concentration = data.frame(fraction[,'species'],
                               as.numeric(fraction[,'f_neutral_lipids']) / gill[,'Vf_gill_neutral_lipids'],
                               as.numeric(fraction[,'f_phospholipids']) / gill[,'Vf_gill_phospholipids'],
                               as.numeric(fraction[,'f_functional_proteins']) / gill[,'Vf_gill_functional_proteins'],
                               as.numeric(fraction[,'f_structural_proteins']) / gill[,'Vf_gill_structural_proteins'],
                               as.numeric(fraction[,'f_nonlipid_organic_matter']) / gill[,'Vf_gill_nonlipid_organic_matter'],
                               as.numeric(fraction[,'f_water']) / gill[,'Vf_gill_water'])
    names(concentration) = c('species', 'C_neutral_lipids' , 'C_phospholipids','C_functional_proteins', 
                             'C_structural_proteins', 'C_nonlipid_organic_matter', 'C_water')
    result = list(fraction = fraction, concentration = concentration)
  }
  return(result)
}


# =========================================================
## Run
# =========================================================
noblood = distribution_workflow(target = "gill without blood")
withblood = distribution_workflow(target = "gill with blood")
gilllonly = distribution_workflow(target = "within gill")


# =========================================================
## Figure 4 (composition)
# =========================================================
# color scheme
organ.color = c(f_adipose = "#A6CEE3", f_gill = "#1F78B4", f_kidney = "#B2DF8A", f_liver = "#33A02C",
                f_blood = "#E31A1C", f_muscle = "#CAB2D6", f_rest = "#6A3D9A",
                C_adipose = "#A6CEE3", C_gill = "#1F78B4", C_kidney = "#B2DF8A", C_liver = "#33A02C",
                C_blood = "#E31A1C", C_muscle = "#CAB2D6", C_rest = "#6A3D9A",
                Vf_neutral_lipids = "#8DD3C7",
                Vf_phospholipids = "#CCEBC5",
                Vf_functional_proteins = "#FDB462",
                Vf_structural_proteins = "#FB8072",
                Vf_nonlipid_organic_matter = "#80B1D3",
                Vf_water = "#FCCDE5",
                f_neutral_lipids = "#8DD3C7",
                f_phospholipids = "#CCEBC5",
                f_functional_proteins = "#FDB462",
                f_structural_proteins = "#FB8072",
                f_nonlipid_organic_matter = "#80B1D3",
                f_water = "#FCCDE5",
                C_neutral_lipids = "#8DD3C7",
                C_phospholipids = "#CCEBC5",
                C_functional_proteins = "#FDB462",
                C_structural_proteins = "#FB8072",
                C_nonlipid_organic_matter = "#80B1D3",
                C_water = "#FCCDE5")

fraction_p = distribution_workflow(target = "gill without blood")$fraction %>% 
  pivot_longer(!species, 
               names_to = "organ",
               values_to = "fraction")

name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

fraction_p <- fraction_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
fraction_p$species = factor(fraction_p$species, levels = name_order)

fraction_pp = ggplot(fraction_p, aes(fill = organ, y = as.numeric(fraction), x = species)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.6, color = "black") +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(size = 0),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = NULL, fill = "Tissue fraction", x = NULL) +
  scale_fill_manual(values = organ.color)


concentration_p =  distribution_workflow(target = "gill without blood")$concentration %>% 
  pivot_longer(-species, 
               names_to = "organ",
               values_to = "concentration")

name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

concentration_p <- concentration_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
concentration_p$species <- factor(concentration_p$species, levels = name_order)

concentration_pp = ggplot(data = concentration_p, aes(y = as.numeric(concentration), x = species, fill = organ, group = organ)) + 
  #geom_line(color='grey', size = 1, linetype='dashed') +
  geom_point(shape = 21, size = 6, alpha = 0.8) +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(colour = "white"),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = "Concentration (mol/L)", fill = "Tissue concentration", x = NULL) +
  scale_x_discrete() +
  scale_fill_manual(values = organ.color)

cairo_pdf("final_figure/figure3.pdf", height = 4, width = 16)
ggarrange(fraction_pp, concentration_pp, labels = c("A.", "B."), nrow = 1, align = "h", font.label = list(size = 20), widths = c(1.1, 1.2)) # add legend later
dev.off()


# =========================================================
## Figure S2
# =========================================================
fraction_p = distribution_workflow(target = "gill with blood")$fraction %>% 
  pivot_longer(!species, 
               names_to = "organ",
               values_to = "fraction")

# name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

fraction_p <- fraction_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
fraction_p$species = factor(fraction_p$species, levels = name_order)

fraction_pp = ggplot(fraction_p, aes(fill = organ, y = as.numeric(fraction), x = species)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.6, color = "black", linewidth = 0.5) +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(size = 0),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = "Fraction (%)", fill = "Tissue fraction", x = NULL) +
  scale_fill_manual(values = organ.color)


concentration_p =  distribution_workflow(target = "gill with blood")$concentration %>% 
  pivot_longer(-species, 
               names_to = "organ",
               values_to = "concentration")

# name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

# concentration_p$species = gsub(".", " ", concentration_p$species, fixed = T)
concentration_p <- concentration_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
concentration_p$species <- factor(concentration_p$species, levels = name_order)

concentration_pp = ggplot(data = concentration_p, aes(y = as.numeric(concentration), x = species, fill = organ, group = organ)) + 
  #geom_line(color='grey', size = 1, linetype='dashed') +
  geom_point(shape = 21, size = 6, alpha = 0.8, stroke = 0.5) +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(colour = "white"),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = "Concentration (mol/L)", fill = "Tissue concentration", x = NULL) +
  scale_x_discrete() +
  scale_fill_manual(values = organ.color)

cairo_pdf("final_figure/figureS2.pdf", height = 4, width = 16)
ggarrange(fraction_pp, concentration_pp, labels = c("A.", "B."), nrow = 1, align = "h", font.label = list(size = 20), widths = c(1.1, 1.2)) # add legend later
dev.off()

# =========================================================
## Figure S3
# =========================================================
fraction_p = distribution_workflow(target = "within gill")$fraction %>% 
  pivot_longer(!species, 
               names_to = "organ",
               values_to = "fraction")

# name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

# fraction_p$species = gsub(".", " ", fraction_p$species, fixed = T)
fraction_p <- fraction_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
fraction_p$species = factor(fraction_p$species, levels = name_order)

fraction_pp = ggplot(fraction_p, aes(fill = organ, y = as.numeric(fraction), x = species)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.6, color = "black") +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(size = 0),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = NULL, fill = "Tissue composition\nfraction", x = NULL) +
  scale_fill_manual(values = organ.color)


concentration_p =  distribution_workflow(target = "within gill")$concentration %>% 
  pivot_longer(-species, 
               names_to = "organ",
               values_to = "concentration")

# name_order = c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA")

# concentration_p$species = gsub(".", " ", concentration_p$species, fixed = T)
concentration_p <- concentration_p %>% mutate(species = ifelse(species == "PFUnDA", "PFUDA", species))
concentration_p$species <- factor(concentration_p$species, levels = name_order)

concentration_pp = ggplot(data = concentration_p, aes(y = as.numeric(concentration), x = species, fill = organ, group = organ)) + 
  #geom_line(color='grey', size = 1, linetype='dashed') +
  geom_point(shape = 21, size = 6, alpha = 0.8) +
  #geom_vline(xintercept = 3.5, linetype = "dashed", color = "black", alpha = 1) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5), 
    legend.title = element_text(face = "bold"), 
    legend.key = element_rect(colour = "white"),
    plot.title = element_text(face = "bold"), 
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)) +
  labs(y = "Concentration (mol/L)", fill = "Tissue compositon\nconcentration", x = NULL) +
  scale_x_discrete() +
  scale_fill_manual(values = organ.color)

cairo_pdf("final_figure/figureS3.pdf", height = 4, width = 17)
ggarrange(fraction_pp, concentration_pp, labels = c("A.", "B."), nrow = 1, align = "h", font.label = list(size = 20), widths = c(1, 1)) # add legend later
dev.off()


