##########################################################################################
# Last edit: 2026-07-26
##########################################################################################
rm(list = ls())

library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(Cairo)

setwd("/Users/huzhiji/OneDrive - Harvard University/Research/PFAAatGill/")

# =========================================================
# RUN for Martin
# update = 12 days, depuration = 38 days, uptake mode = lab
# PFHxS initial concentration = 1.4 * 10^-3
# PFOS initial concentration = 0.35 * 10^-3
# PFOA initial concentration = 1.5 * 10^-3
# PFDA initial concentration = 0.71 * 10^-3
# PFUDA initial concentration = 0.48 * 10^-3
# =========================================================
days_uptake <- 12 # day
days_depuration <- 38  # day

# =========================================================
# PFHxS
initial_concentration <- 1.4 * 10^-3 # ug/cm3

tic()

results <- run_pfas_workflow(
  excel_path  = excel_path,
  sheet_name  = sheet_name,
  out_dir     = "PBTK_PFHxS_outputs",
  n_iter      = n_iter,
  verbose     = verbose,
  uptake_mode = "lab",
  target      = "PFHxS"
)

toc()
# =========================================================

# =========================================================
# PFOS
initial_concentration <- 0.35 * 10^-3 # ug/cm3

tic()

results <- run_pfas_workflow(
  excel_path  = excel_path,
  sheet_name  = sheet_name,
  out_dir     = "PBTK_PFOS_outputs",
  n_iter      = n_iter,
  verbose     = verbose,
  uptake_mode = "lab",
  target      = "PFOS"
)

toc()
# =========================================================

# =========================================================
# PFOA
initial_concentration <- 1.5 * 10^-3 # ug/cm3

tic()

results <- run_pfas_workflow(
  excel_path  = excel_path,
  sheet_name  = sheet_name,
  out_dir     = "PBTK_PFOA_outputs",
  n_iter      = n_iter,
  verbose     = verbose,
  uptake_mode = "lab",
  target      = "PFOA"
)

toc()
# =========================================================

# =========================================================
# PFDA
initial_concentration <- 0.71 * 10^-3 # ug/cm3

tic()

results <- run_pfas_workflow(
  excel_path  = excel_path,
  sheet_name  = sheet_name,
  out_dir     = "PBTK_PFDA_outputs",
  n_iter      = n_iter,
  verbose     = verbose,
  uptake_mode = "lab",
  target      = "PFDA"
)

toc()
# =========================================================

# =========================================================
# PFUDA
initial_concentration <- 0.48 * 10^-3 # ug/cm3

tic()

results <- run_pfas_workflow(
  excel_path  = excel_path,
  sheet_name  = sheet_name,
  out_dir     = "PBTK_PFUDA_outputs",
  n_iter      = n_iter,
  verbose     = verbose,
  uptake_mode = "lab",
  target      = "PFUDA"
)

toc()
# =========================================================

# =========================================================
## Figure S1
# =========================================================
# scheme
color.id = c(C_water = "blue", C_blood = "red")
fill.id = c(C_water = "darkblue", C_blood = "darkred")

scientific_10x <- function(x) {parse(text = paste0("10^", formatC(log10(x), format = "f", digits = 0)))}

# =========================================================
# PFHxS
# =========================================================
martin.fish.PFHxS = data.frame(
  medium = rep("C_blood", 15),
  time = c(9, 18, 36, 72, 144, 288, 4.5+288, 9+288, 18+288, 36+288, 72+288, 144+288, 288+288, 456+288, 792+288),
  concentration = c(0.00360, 0.00475, 0.0125, 0.0197, 0.0344, 0.0888, 0.0742, 0.0544, 0.0620, 0.0600, 0.0562, 0.0418, 0.0165, 0.0157, 0.00856)
)

martin.water.PFHxS = data.frame(
  medium = rep("C_water", 10),
  time = c(0.25, 4.5, 12, 18, 36, 72, 144, 197, 244, 288),
  concentration = c(0.00168, 0.00138, 0.00138, 0.00127, 0.00127, 0.00122, 0.00131, 0.00161, 0.00136, 0.00156)
)

martin.PFHxS = rbind(martin.fish.PFHxS, martin.water.PFHxS)

PFHxS_model_fish = read.xlsx("PBTK_PFHxS_outputs/Fish_PBTK_MC_Results_PFHxS.xlsx", sheet = "Time_Course_Fish_Blood") %>% 
  select(time = time, C_blood = mean, C_blood_low = ci_lower, C_blood_hi = ci_upper) %>%
  mutate(hour = time*24) 

PFHxS_model_water = read.xlsx("PBTK_PFHxS_outputs/Fish_PBTK_MC_Results_PFHxS.xlsx", sheet = "Time_Course_Water") %>% 
  select(time = time, C_water = mean, C_water_low = ci_lower, C_water_hi = ci_upper) %>% 
  mutate(hour = time*24)

# day 12 comparison
PFHxS_model_fish %>% filter(time == 12)

# Save PDF
pdf("final_figure/PFHxS_evaluation.pdf", height = 5, width = 6)
ggplot() +
  geom_vline(xintercept = 288, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 1080, linetype = "dashed", color = "black", alpha = 0.5) +
  # fish blood
  geom_line(data = PFHxS_model_fish, aes(x = as.numeric(hour), y = C_blood), color = "red") +
  geom_ribbon(data = PFHxS_model_fish, 
              aes(x = as.numeric(hour), ymin = C_blood_low, ymax = C_blood_hi), 
              fill = "red", alpha = 0.2) +
  # water
  geom_line(data = PFHxS_model_water, aes(x = as.numeric(hour), y = C_water), color = "blue", show.legend = F) +
  geom_ribbon(data = PFHxS_model_water, 
              aes(x = as.numeric(hour), ymin = C_water_low, ymax = C_water_hi), 
              fill = "blue", alpha = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_point(data = martin.PFHxS, aes(x = time, y = concentration, fill = medium), shape = 21, color = "black", size = 3, show.legend = F) +
  annotate("text", x = 144, y = 10, label = "uptake", size = 8, fontface = "bold") +
  annotate("text", x = 684, y = 10, label = "depuration", size = 8, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_log10(limits = c(1e-07, 10), 
                breaks = c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10), 
                labels = scientific_10x) + 
  labs(x = "Time (h)", y = "Concentration (ug/mL)", 
       color = "Our Model", 
       fill = "Martin et al. (2003)", 
       title = "PFHxS") + # (ηpfc = 6)
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    plot.title = element_text(size = 25, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = fill.id)
dev.off()

# =========================================================
# PFOS
# =========================================================
martin.fish.PFOS = data.frame(
  medium = rep("C_blood", 15),
  time = c(9, 18, 36, 72, 144, 288, 4.5+288, 9+288, 18+288, 36+288, 72+288, 144+288, 288+288, 456+288, 792+288),
  concentration = c(0.0274, 0.0452, 0.132, 0.210, 0.459, 1.22, 1.08, 0.924, 1.03, 0.967, 0.918, 0.785, 0.270, 0.320, 0.206)
)

martin.water.PFOS = data.frame(
  medium = rep("C_water", 10),
  time = c(0.25, 4.5, 12, 18, 36, 72, 144, 197, 244, 288),
  concentration = c(0.000502, 0.000467, 0.000297, 0.000239, 0.000300, 0.000268, 0.000291, 0.000388, 0.000303, 0.000439)
)

martin.PFOS = rbind(martin.fish.PFOS, martin.water.PFOS)

PFOS_model_fish = read.xlsx("PBTK_PFOS_outputs/Fish_PBTK_MC_Results_PFOS.xlsx", sheet = "Time_Course_Fish_Blood") %>% 
  select(time = time, C_blood = mean, C_blood_low = ci_lower, C_blood_hi = ci_upper) %>%
  mutate(hour = time*24) 

PFOS_model_water = read.xlsx("PBTK_PFOS_outputs/Fish_PBTK_MC_Results_PFOS.xlsx", sheet = "Time_Course_Water") %>% 
  select(time = time, C_water = mean, C_water_low = ci_lower, C_water_hi = ci_upper) %>% 
  mutate(hour = time*24)

# day 12 comparison
PFOS_model_fish %>% filter(time == 12)

# Save PDF
pdf("final_figure/PFOS_evaluation.pdf", height = 5, width = 6)
ggplot() +
  geom_vline(xintercept = 288, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 1080, linetype = "dashed", color = "black", alpha = 0.5) +
  # fish blood
  geom_line(data = PFOS_model_fish, aes(x = as.numeric(hour), y = C_blood), color = "red") +
  geom_ribbon(data = PFOS_model_fish, 
              aes(x = as.numeric(hour), ymin = C_blood_low, ymax = C_blood_hi), 
              fill = "red", alpha = 0.2) +
  # water
  geom_line(data = PFOS_model_water, aes(x = as.numeric(hour), y = C_water), color = "blue", show.legend = F) +
  geom_ribbon(data = PFOS_model_water, 
              aes(x = as.numeric(hour), ymin = C_water_low, ymax = C_water_hi), 
              fill = "blue", alpha = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_point(data = martin.PFOS, aes(x = time, y = concentration, fill = medium), shape = 21, color = "black", size = 3, show.legend = F) +
  annotate("text", x = 144, y = 10, label = "uptake", size = 8, fontface = "bold") +
  annotate("text", x = 684, y = 10, label = "depuration", size = 8, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_log10(limits = c(1e-07, 10), 
                breaks = c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10), 
                labels = scientific_10x) + 
  labs(x = "Time (h)", y = "Concentration (ug/mL)", 
       color = "Our Model", 
       fill = "Martin et al. (2003)", 
       title = "PFOS") + # (ηpfc = 8)
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    plot.title = element_text(size = 25, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = fill.id)
dev.off()



# =========================================================
# PFOA
# =========================================================
martin.fish.PFOA = data.frame(
  medium = rep("C_blood", 14),
  time = c(4.5, 9, 18, 36, 72, 144, 288, 4.5+288, 9+288, 18+288, 36+288, 72+288, 144+288, 288+288),
  concentration = c(0.00493, 0.00641, 0.00800, 0.0121, 0.0167, 0.0207, 0.0398, 0.0365, 0.0306, 0.0273, 0.0203, 0.0185, 0.0106, 0.00253)
)

martin.water.PFOA = data.frame(
  medium = rep("C_water", 10),
  time = c(0.25, 4.5, 12, 18, 36, 72, 144, 197, 244, 288),
  concentration = c(0.00174, 0.00131, 0.00146, 0.00138, 0.00131, 0.00148, 0.00149, 0.00187, 0.00149, 0.00177)
)

martin.PFOA = rbind(martin.fish.PFOA, martin.water.PFOA)

PFOA_model_fish = read.xlsx("PBTK_PFOA_outputs/Fish_PBTK_MC_Results_PFOA.xlsx", sheet = "Time_Course_Fish_Blood") %>% 
  select(time = time, C_blood = mean, C_blood_low = ci_lower, C_blood_hi = ci_upper) %>%
  mutate(hour = time*24) 

PFOA_model_water = read.xlsx("PBTK_PFOA_outputs/Fish_PBTK_MC_Results_PFOA.xlsx", sheet = "Time_Course_Water") %>% 
  select(time = time, C_water = mean, C_water_low = ci_lower, C_water_hi = ci_upper) %>% 
  mutate(hour = time*24)

# day 12 comparison
PFOA_model_fish %>% filter(time == 12)

# Save PDF
pdf("final_figure/PFOA_evaluation.pdf", height = 5, width = 6)
ggplot() +
  geom_vline(xintercept = 288, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 1080, linetype = "dashed", color = "black", alpha = 0.5) +
  # fish blood
  geom_line(data = PFOA_model_fish, aes(x = as.numeric(hour), y = C_blood), color = "red") +
  geom_ribbon(data = PFOA_model_fish, 
              aes(x = as.numeric(hour), ymin = C_blood_low, ymax = C_blood_hi), 
              fill = "red", alpha = 0.2) +
  # water
  geom_line(data = PFOA_model_water, aes(x = as.numeric(hour), y = C_water), color = "blue", show.legend = F) +
  geom_ribbon(data = PFOA_model_water, 
              aes(x = as.numeric(hour), ymin = C_water_low, ymax = C_water_hi), 
              fill = "blue", alpha = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_point(data = martin.PFOA, aes(x = time, y = concentration, fill = medium), shape = 21, color = "black", size = 3, show.legend = F) +
  annotate("text", x = 144, y = 10, label = "uptake", size = 8, fontface = "bold") +
  annotate("text", x = 684, y = 10, label = "depuration", size = 8, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_log10(limits = c(1e-07, 10), 
                breaks = c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10), 
                labels = scientific_10x) + 
  labs(x = "Time (h)", y = "Concentration (ug/mL)", 
       color = "Our Model", 
       fill = "Martin et al. (2003)", 
       title = "PFOA") + #  (ηpfc = 7)
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    plot.title = element_text(size = 25, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = fill.id)
dev.off()



# =========================================================
# PFDA
# =========================================================
martin.fish.PFDA = data.frame(
  medium = rep("C_blood", 16),
  time = c(4.5, 9, 18, 36, 72, 144, 288, 4.5+288, 9+288, 18+288, 36+288, 72+288, 144+288, 288+288, 456+288, 792+288),
  concentration = c(7.075962325290E-3, 3.178320748968E-2, 6.874833276451E-2, 1.487220666641E-1, 2.849006688257E-1,
                    5.763901790937E-1, 1.509360712268E0, 1.354981086906E0, 1.256219734267E0, 1.407502900709E0,
                    1.187179939081E0, 1.210415384165E0, 1.081366431300E0, 3.477216282441E-1, 4.373643196939E-1,
                    2.395838008604E-1))

martin.water.PFDA = data.frame(
  medium = rep("C_water", 10),
  time = c(0.25, 4.5, 12, 18, 36, 72, 144, 197, 244, 288),
  concentration = c(8.121432482833E-4, 7.814148212151E-4, 7.374560936922E-4, 5.963631195890E-4, 7.164782029116E-4,
                    8.781488037422E-4, 6.391578740194E-4, 9.501555105388E-4, 6.039538261584E-4, 8.805448978580E-4))

martin.PFDA = rbind(martin.fish.PFDA, martin.water.PFDA)

PFDA_model_fish = read.xlsx("PBTK_PFDA_outputs/Fish_PBTK_MC_Results_PFDA.xlsx", sheet = "Time_Course_Fish_Blood") %>% 
  select(time = time, C_blood = mean, C_blood_low = ci_lower, C_blood_hi = ci_upper) %>%
  mutate(hour = time*24) 

PFDA_model_water = read.xlsx("PBTK_PFDA_outputs/Fish_PBTK_MC_Results_PFDA.xlsx", sheet = "Time_Course_Water") %>% 
  select(time = time, C_water = mean, C_water_low = ci_lower, C_water_hi = ci_upper) %>% 
  mutate(hour = time*24)

# day 12 comparison
PFDA_model_fish %>% filter(time == 12)

# Save PDF
pdf("final_figure/PFDA_evaluation.pdf", height = 5, width = 6)
ggplot() +
  geom_vline(xintercept = 288, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 1080, linetype = "dashed", color = "black", alpha = 0.5) +
  # fish blood
  geom_line(data = PFDA_model_fish, aes(x = as.numeric(hour), y = C_blood), color = "red") +
  geom_ribbon(data = PFDA_model_fish, 
              aes(x = as.numeric(hour), ymin = C_blood_low, ymax = C_blood_hi), 
              fill = "red", alpha = 0.2) +
  # water
  geom_line(data = PFDA_model_water, aes(x = as.numeric(hour), y = C_water), color = "blue", show.legend = F) +
  geom_ribbon(data = PFDA_model_water, 
              aes(x = as.numeric(hour), ymin = C_water_low, ymax = C_water_hi), 
              fill = "blue", alpha = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_point(data = martin.PFDA, aes(x = time, y = concentration, fill = medium), shape = 21, color = "black", size = 3, show.legend = F) +
  annotate("text", x = 144, y = 10, label = "uptake", size = 8, fontface = "bold") +
  annotate("text", x = 684, y = 10, label = "depuration", size = 8, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_log10(limits = c(1e-07, 10), 
                breaks = c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10), 
                labels = scientific_10x) + 
  labs(x = "Time (h)", y = "Concentration (ug/mL)", 
       color = "Our Model", 
       fill = "Martin et al. (2003)", 
       title = "PFDA") + #  (ηpfc = 9)
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    plot.title = element_text(size = 25, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = fill.id)

dev.off()


# =========================================================
# PFUDA
# =========================================================
martin.fish.PFUDA = data.frame(
  medium = rep("C_blood", 16),
  time = c(4.5, 9, 18, 36, 72, 144, 288, 4.5+288, 9+288, 18+288, 36+288, 72+288, 144+288, 288+288, 456+288, 792+288),
  concentration = c(3.264648818548E-2, 7.030044351550E-2, 1.965858681621E-1, 3.655787413278E-1, 6.911923138282E-1,
                    1.265418126404E0, 2.866850583951E0, 2.634742506518E0, 2.948244957796E0, 2.515826892764E0,
                    2.483583699674E0, 2.743646962894E0, 2.237375769520E0, 9.492630991346E-1, 1.122163231481E0,
                    7.333343361731E-1))

martin.water.PFUDA = data.frame(
  medium = rep("C_water", 10),
  time = c(0.25, 4.5, 12, 18, 36, 72, 144, 197, 244, 288),
  concentration = c(3.953044549936E-4, 4.400936582785E-4, 4.898963699052E-4, 3.317929627165E-4, 3.552877062059E-4,
                    6.015592091577E-4, 4.198647670597E-4, 5.412718831272E-4, 3.776671377149E-4, 5.634212447714E-4))

martin.PFUDA = rbind(martin.fish.PFUDA, martin.water.PFUDA)


PFUDA_model_fish = read.xlsx("PBTK_PFUDA_outputs/Fish_PBTK_MC_Results_PFUDA.xlsx", sheet = "Time_Course_Fish_Blood") %>% 
  select(time = time, C_blood = mean, C_blood_low = ci_lower, C_blood_hi = ci_upper) %>%
  mutate(hour = time*24) 

PFUDA_model_water = read.xlsx("PBTK_PFUDA_outputs/Fish_PBTK_MC_Results_PFUDA.xlsx", sheet = "Time_Course_Water") %>% 
  select(time = time, C_water = mean, C_water_low = ci_lower, C_water_hi = ci_upper) %>% 
  mutate(hour = time*24)

# day 12 comparison
PFUDA_model_fish %>% filter(time == 12)

# Save PDF
pdf("final_figure/PFUDA_evaluation.pdf", height = 5, width = 6)
ggplot() +
  geom_vline(xintercept = 288, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 1080, linetype = "dashed", color = "black", alpha = 0.5) +
  # fish blood
  geom_line(data = PFUDA_model_fish, aes(x = as.numeric(hour), y = C_blood), color = "red") +
  geom_ribbon(data = PFUDA_model_fish, 
              aes(x = as.numeric(hour), ymin = C_blood_low, ymax = C_blood_hi), 
              fill = "red", alpha = 0.2) +
  # water
  geom_line(data = PFUDA_model_water, aes(x = as.numeric(hour), y = C_water), color = "blue", show.legend = F) +
  geom_ribbon(data = PFUDA_model_water, 
              aes(x = as.numeric(hour), ymin = C_water_low, ymax = C_water_hi), 
              fill = "blue", alpha = 0.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_point(data = martin.PFUDA, aes(x = time, y = concentration, fill = medium), shape = 21, color = "black", size = 3, show.legend = F) +
  annotate("text", x = 144, y = 10, label = "uptake", size = 8, fontface = "bold") +
  annotate("text", x = 684, y = 10, label = "depuration", size = 8, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_log10(limits = c(1e-07, 10), 
                breaks = c(1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1, 10), 
                labels = scientific_10x) + 
  labs(x = "Time (h)", y = "Concentration (ug/mL)", 
       color = "Our Model", 
       fill = "Martin et al. (2003)", 
       title = "PFUDA") + # (ηpfc = 10)
  theme_bw() +
  theme(
    axis.text = element_text(size = 20, color = "black"), 
    legend.text = element_text(size = 20), 
    legend.title = element_text(size = 20, face = "bold"), 
    plot.title = element_text(size = 25, face = "bold"), 
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1)) +
  scale_fill_manual(values = fill.id)

dev.off()

