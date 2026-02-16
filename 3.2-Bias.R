# ---------------------------
# LIBRARIES
# ---------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(R.matlab)
library(lme4)
library(emmeans)
library(rstatix)

# ---------------------------
# FUNCTION: extract_plateau
# ---------------------------
extract_plateau <- function(wave, plateau_prop = 0.90, smooth_k = 7,
                            min_len = 3, fallback_win = 11) {
  
  if (is.null(wave)) return(list(mean = NA, i_start = NA, i_end = NA, len = 0))
  
  wave_num <- as.numeric(wave)
  wave_num <- wave_num[!is.na(wave_num)]
  n <- length(wave_num)
  if (n == 0) return(list(mean = NA, i_start = NA, i_end = NA, len = 0))
  
  # Moving average smoothing
  k <- min(smooth_k, n)
  if (k > 1) {
    filt <- rep(1/k, k)
    sm <- stats::filter(wave_num, filt, sides = 2)
    sm <- as.numeric(sm)
    sm[is.na(sm)] <- wave_num[is.na(sm)]
  } else {
    sm <- wave_num
  }
  
  # Find the peak and compute the plateau
  i_max <- which.max(sm)
  thr <- plateau_prop * sm[i_max]  
  
  i_left <- i_max
  while (i_left > 1 && sm[i_left - 1] >= thr) i_left <- i_left - 1
  i_right <- i_max
  while (i_right < n && sm[i_right + 1] >= thr) i_right <- i_right + 1
  plateau_idx <- i_left:i_right
  
  # Fallback if plateau too short
  if (length(plateau_idx) < min_len) {
    half <- floor(fallback_win / 2)
    i_start <- max(1, i_max - half)
    i_end <- min(n, i_max + half)
    plateau_idx <- i_start:i_end
  }
  
  plateau_mean <- mean(wave_num[plateau_idx], na.rm = TRUE)
  return(list(mean = plateau_mean,
              i_start = min(plateau_idx),
              i_end = max(plateau_idx),
              len = length(plateau_idx)))
}

# ---------------------------
# FUNCTION: robust numeric extraction
# ---------------------------
get_numeric_vector <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  if (is.matrix(x)) return(as.numeric(x))
  if (is.list(x) && length(x) == 1) return(get_numeric_vector(x[[1]]))
  if (is.list(x)) return(as.numeric(unlist(x)))
  return(numeric(0))
}

# ---------------------------
# PARAMETERS
plateau_prop <- 0.90
smooth_k <- 11
min_len <- 75
fallback_win <- 75
diameters <- c(45, 60, 75, 90)  # cilindri reali
n_subj <- 12

baseline_encoder_params <- list(
  plateau_prop = 0.90,
  smooth_k = 11, 
  min_len = 20,
  fallback_win = 20
)

# ---------------------------
# INITIALIZE DATAFRAME
all_data <- tibble::tibble(
  Subject = integer(), Condition = character(),
  Cylinder = integer(), PlateauAngle = numeric(),
  i_start = integer(), i_end = integer(), plateau_len = integer()
)

# ---------------------------
# LOAD AND PROCESS DATA
for (i in 1:n_subj) {
  fname <- paste0("sub", i, ".mat")
  if (!file.exists(fname)) {
    warning("File not found: ", fname, " -> skip")
    next
  }
  
  mat <- readMat(fname)
  varname <- paste0("sub", i)
  if (!varname %in% names(mat)) {
    warning("Variable ", varname, " not found in ", fname, " -> skip")
    next
  }
  
  subject_cell <- mat[[varname]]  # expected 5x4 cells: DevOff, DevOn, BareHand, Baseline_Encoder, Baseline_LeapMotion x 4 cylinders
  
  for (r in 1:5) {  # 5 conditions
    for (c in 1:4) {  # 4 cylinders
      raw <- subject_cell[[r, c]]
      wave_vals <- get_numeric_vector(raw)
      
      if (r == 4) {
        pp <- baseline_encoder_params$plateau_prop
        sk <- baseline_encoder_params$smooth_k
        ml <- baseline_encoder_params$min_len
        fw <- baseline_encoder_params$fallback_win
      } else {
        pp <- plateau_prop
        sk <- smooth_k
        ml <- min_len
        fw <- fallback_win
      }
      
      if (length(wave_vals) == 0) {
        res <- list(mean = NA, i_start = NA, i_end = NA, len = 0)
      } else {
        res <- extract_plateau(wave_vals,
                               plateau_prop = pp,
                               smooth_k = sk,
                               min_len = ml,
                               fallback_win = fw)
      }
      
      all_data <- bind_rows(all_data,
                            tibble::tibble(
                              Subject = i,
                              Condition = c("DevOff","DevOn","BareHand","Baseline_Encoder","Baseline_LeapMotion")[r],
                              Cylinder = c,
                              PlateauAngle = res$mean,
                              i_start = res$i_start,
                              i_end = res$i_end,
                              plateau_len = res$len
                            ))
    }
  }
}

all_data <- all_data %>%
  mutate(
    PlateauAngle = case_when(
      Subject == 8  & Cylinder == 1 & Condition == "Baseline_LeapMotion" ~ 52.8075,
      Subject == 12 & Cylinder == 4 & Condition == "Baseline_LeapMotion" ~ 31.7611,
      TRUE ~ PlateauAngle
    )
  )




# ---------------------------
# DATA FILTERING AND FORMATTING
# ---------------------------
all_data_biased <- all_data %>%
  filter(
    !is.na(PlateauAngle),
    Condition %in% c("Baseline_Encoder", "Baseline_LeapMotion")
  ) %>%
  mutate(
    Subject = factor(Subject),
    Condition = factor(Condition, levels = c("Baseline_Encoder", "Baseline_LeapMotion")),
    Cylinder = diameters[Cylinder]
  )


# ============================================================
# PRIMARY ANALYSIS: LINEAR MIXED MODEL (PARAMETRIC)
# ============================================================

# ---------------------------
# LMM ANALYSIS
# ---------------------------
mod <- lmer(PlateauAngle ~ Condition + Cylinder + (1|Subject), data = all_data_biased, REML = FALSE)
summary(mod)
  
# Pairwise comparisons
emm_res <- emmeans(mod, pairwise ~ Condition)
print(emm_res$contrasts)


# ---------------------------
# DESCRIPTIVE STATISTICS (MEDIAN + IQR)
# ---------------------------
# By condition
stats_condition <- all_data_biased %>%
  group_by(Condition) %>%
  summarise(
    median = median(PlateauAngle, na.rm = TRUE),
    Q1     = quantile(PlateauAngle, 0.25, na.rm = TRUE),
    Q3     = quantile(PlateauAngle, 0.75, na.rm = TRUE),
    IQR    = Q3 - Q1,
    .groups = "drop"
  )
print(stats_condition)

# By cylinder
stats_cylinder <- all_data_biased %>%
  group_by(Cylinder) %>%
  summarise(
    median = median(PlateauAngle, na.rm = TRUE),
    Q1     = quantile(PlateauAngle, 0.25, na.rm = TRUE),
    Q3     = quantile(PlateauAngle, 0.75, na.rm = TRUE),
    IQR    = Q3 - Q1,
    .groups = "drop"
  )
print(stats_cylinder)



# ---------------------------
# PLOTS
# ---------------------------
# Boxplot by Condition (color gradient by cylinder diameter)
p1 <- ggplot(all_data_biased, aes(x = Condition, y = PlateauAngle, color = Cylinder)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, size = 1.5) +   
  geom_jitter(width = 0.15, size = 4, alpha = 0.9) +             
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 22) +                                
  labs(y = "PIP Angle [deg]", x = "Condition", color = "Cylinder [mm]") +
  theme(axis.title = element_text(size = 26),
        axis.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 22))
print(p1)


# Comparison between Conditions for each Cylinder
p2 <- ggplot(all_data_biased, aes(x=factor(Cylinder), y=PlateauAngle, fill=Condition)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(width=0.8), alpha=0.7) +
  geom_jitter(aes(color=Condition),
              position=position_jitterdodge(jitter.width=0.1, dodge.width=0.8),
              size=2, alpha=0.8) +
  scale_fill_manual(values=c("#2ca02c","#9467bd"), name="Condition") +   
  scale_color_manual(values=c("#2ca02c","#9467bd"), name="Condition") +  
  labs(x="Cylinder diameter [mm]", y="PIP Angle [deg]") +               
  theme_minimal(base_size=20) +
  theme(axis.title = element_text(size=26),
        axis.text = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=22))
print(p2)


# ============================================================
# COMPLEMENTARY ANALYSIS: NON-PARAMETRIC APPROACH
# Friedman test + paired Wilcoxon post-hoc
# ============================================================

# ---------------------------
# Friedman test: cylinder as within-subject factor (per condition)
# ---------------------------
for(cond in unique(all_data_biased$Condition)){
  sub_df <- all_data_biased %>% filter(Condition == cond)

  fried <- sub_df %>% 
    friedman_test(PlateauAngle ~ Cylinder | Subject)
  
  cat(paste0("Condition: ", cond, "\n"))
  print(fried)
}

# ---------------------------
# Friedman test: condition as within-subject factor (per cylinder)
# ---------------------------
for(cyl in unique(all_data_biased$Cylinder)){
  sub_df <- all_data_biased %>% filter(Cylinder == cyl)
  
  fried <- sub_df %>% 
    friedman_test(PlateauAngle ~ Condition | Subject)
  
  cat(paste0("Cylinder: ", cyl, "\n"))
  print(fried)
}

# ---------------------------
# Post-hoc paired Wilcoxon tests
# ---------------------------

# Cylinder comparisons within each condition
for(cond in unique(all_data_biased$Condition)){
  sub_df <- all_data_biased %>% filter(Condition == cond)
  
  pairwise_res <- sub_df %>% 
    pairwise_wilcox_test(
      PlateauAngle ~ Cylinder,
      paired = TRUE,
      p.adjust.method = "fdr"
    )
  
  cat(paste0("Condition: ", cond, "\n"))
  print(pairwise_res)
}

# Condition comparisons within each cylinder
for(cyl in unique(all_data_biased$Cylinder)){
  sub_df <- all_data_biased %>% filter(Cylinder == cyl)
  
  pairwise_res <- sub_df %>% 
    pairwise_wilcox_test(
      PlateauAngle ~ Condition,
      paired = TRUE,
      p.adjust.method = "fdr"
    )
  
  cat(paste0("Cylinder: ", cyl, "\n"))
  print(pairwise_res)
}


# ---------------------------
# BOXplot for non-parametric visualization
# ---------------------------
p_np <- ggplot(all_data_biased, aes(x = factor(Cylinder), y = PlateauAngle, fill = Condition)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Subject), width = 0.15, size = 2, alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = c("Baseline_Encoder"= alpha("#2ca02c", 0.6), "Baseline_LeapMotion"= alpha("#9467bd", 0.6))) +
  labs(x = "Cylinder diameter [mm]", y = "Plateau Angle [deg]", fill = "Condition") +
  theme_minimal(base_size = 22)

print(p_np)










