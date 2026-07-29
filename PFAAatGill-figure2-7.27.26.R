##########################################################################################
# Last edit: 2026-07-27
##########################################################################################
rm(list = ls())

library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(Cairo)

setwd("/Users/huzhiji/OneDrive - Harvard University/Research/PFAAatGill/")

# =========================================================
# SETTINGS
# =========================================================
lab_dir <- "PBTK_lab_outputs"   # BCF source

name_order <- c("PFBS","PFHxA","PFHpA","PFHxS","PFOA","PFOS","PFNA","PFDA","PFUDA")

# Choose model uncertainty display: SD or 95% CI
use_ci <- T  # TRUE = 95% CI, FALSE = +/- SD

# =========================================================
# Published data
# =========================================================
PFAS_df <- data.frame(name = name_order)

PFAS_df$logBCF_buck = c(1.06, 0.98, 1.26,
                        2.07, 1.38, 3.01,
                        2.78, 3.79, 3.57)
PFAS_df$logBCF_buck_sd = c(0.49, 0.30, NA,
                           0.25, 0.61, 0.66,
                           0.51, 0.48, 0.31)

PFAS_df$name <- factor(PFAS_df$name, levels = name_order)

# =========================================================
# Find corresponding output
# =========================================================
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

read_summary_sheet <- function(chemical_name, dir_path) {
  xlsx_path <- file.path(dir_path, paste0("Fish_PBTK_MC_Results_", safe_name(chemical_name), ".xlsx"))
  if (!file.exists(xlsx_path)) {
    warning("Missing model output file: ", xlsx_path)
    return(NULL)
  }
  read.xlsx(xlsx_path, sheet = "Summary")
}

# Read log10(BCF) from a folder (lab)
read_logBCF_from_dir <- function(name_order, dir_path) {
  rows <- lapply(name_order, function(chem) {
    summ <- read_summary_sheet(chem, dir_path)
    if (is.null(summ)) return(NULL)
    r <- summ %>% filter(Metric == "log10(BCF)")
    if (nrow(r) == 0) {
      warning("No log10(BCF) row in Summary for ", chem, " in ", dir_path)
      return(NULL)
    }
    data.frame(
      name = chem,
      logBCF_model_median = r$Median[1],
      logBCF_model_sd   = r$SD[1],
      logBCF_model_lo   = r$CI_lower_2.5[1],
      logBCF_model_hi   = r$CI_upper_97.5[1]
    )
  })
  bind_rows(rows)
}

# =========================================================
# Read model outputs
# =========================================================
bcf_lab <- read_logBCF_from_dir(name_order, lab_dir)

PFAS_df <- PFAS_df %>% left_join(bcf_lab, by = "name")

PFAS_df$name = factor(c("PFBS" ,"PFHxA", "PFHpA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFDA", "PFUDA"), levels = name_order)

# Quick sanity checks
if (any(is.na(PFAS_df$logBCF_model_median))) message("Some LAB BCF outputs missing in: ", lab_dir)

# =========================================================
# Plot
#   - BCF points from LAB folder
# =========================================================
BCF_long <- bind_rows(
  PFAS_df %>%
    transmute(
      name,
      source = "Buck",
      est    = logBCF_buck,                  # median
      ymin   = logBCF_buck - logBCF_buck_sd, # SD
      ymax   = logBCF_buck + logBCF_buck_sd
    ),
  PFAS_df %>%
    transmute(
      name,
      source = "Model",
      est    = logBCF_model_median,         # median
      ymin   = logBCF_model_lo,             # CI
      ymax   = logBCF_model_hi            
    )
) %>%
  mutate(source = factor(source, levels = c("Buck","Model")))


pd <- position_dodge(width = 0.6)

BCF_p <- ggplot(BCF_long, aes(x = name, y = est, group = source)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = source), width = 0.35, linewidth = 1, position = pd, show.legend = F) +
  geom_point(aes(fill = source), shape = 21, size = 3, color = "black", stroke = 0.5, position = pd,  show.legend = F) +
  scale_fill_manual(values = c("Buck" = "#992224", "Model" = "#f1a805")) +
  scale_color_manual(values = c("Buck" = "#EF8B67", "Model" = "#f2d6a1")) +
  theme_bw(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title  = element_text(face = "bold"),
    panel.grid  = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 1)
  ) +
  labs(x = NULL, y = "Apparent logBCF") +
  coord_cartesian(ylim = c(-0.3, 5))

# =========================================================
# Save PDF (add legend later in Adobe illustrator)
# =========================================================
cairo_pdf("./final_figure/figure2.pdf", height = 4, width = 5)
BCF_p
dev.off()

