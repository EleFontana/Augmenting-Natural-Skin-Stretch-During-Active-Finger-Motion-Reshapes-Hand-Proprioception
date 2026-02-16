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
# PREPARE DATA
# ---------------------------
df_wide <- all_data %>%
  filter(Condition %in% c(
    "BareHand",
    "Baseline_LeapMotion",
    "Baseline_Encoder",
    "DevOn"
  )) %>%
  select(Subject, Cylinder, Condition, PlateauAngle) %>%
  mutate(
    Subject = factor(Subject),
    Cylinder_mm = factor(
      diameters[Cylinder],
      levels = diameters,
      labels = diameters
    )
  ) %>%
  select(-Cylinder) %>%
  pivot_wider(
    names_from = Condition,
    values_from = PlateauAngle
  )

# ---------------------------
# ERROR COMPUTATION
# ---------------------------
df_wide <- df_wide %>%
  mutate(
    Accuracy    = Baseline_LeapMotion - BareHand,
    SkinStretchEffect = Baseline_Encoder - DevOn,
    Diff              = SkinStretchEffect - Accuracy
  )

# ---------------------------
# Accuracy per subject
# ---------------------------
accuracy_subject <- df_wide %>%
  group_by(Subject) %>%
  summarise(
    mean_error = mean(Accuracy, na.rm = TRUE),
    sd_error   = sd(Accuracy, na.rm = TRUE),
    .groups = "drop"
  )

print(accuracy_subject)

# ---------------------------
# Skin stretch effect per subject
# ---------------------------
skin_subject <- df_wide %>%
  group_by(Subject) %>%
  summarise(
    mean_skin_effect = mean(SkinStretchEffect, na.rm = TRUE),
    sd_skin_effect   = sd(SkinStretchEffect, na.rm = TRUE),
    .groups = "drop"
  )

print(skin_subject)

# ---------------------------
# NORMALITY CHECK (per cylinder)
# ---------------------------
normality_check <- df_wide %>%
  group_by(Cylinder_mm) %>%
  summarise(
    n = sum(!is.na(Diff)),
    shapiro_p = ifelse(n >= 3, shapiro.test(Diff)$p.value, NA),
    .groups = "drop"
  )

print(normality_check)

# ---------------------------
# PAIRED T-TEST PER CYLINDER
# ---------------------------
cylinder_list <- levels(df_wide$Cylinder_mm)
results <- tibble()

for (cyl in cylinder_list) {
  sub <- df_wide %>% filter(Cylinder_mm == cyl)
  
  t_res <- t.test(
    sub$SkinStretchEffect,
    sub$Accuracy,
    paired = TRUE
  )
  
  results <- bind_rows(
    results,
    tibble(
      Cylinder_mm = cyl,
      Mean_Accuracy   = mean(sub$Accuracy, na.rm = TRUE),
      SD_Accuracy     = sd(sub$Accuracy, na.rm = TRUE),
      Mean_SkinStretch = mean(sub$SkinStretchEffect, na.rm = TRUE),
      SD_SkinStretch   = sd(sub$SkinStretchEffect, na.rm = TRUE),
      t_value = as.numeric(t_res$statistic),
      df      = as.numeric(t_res$parameter),
      p_value = t_res$p.value
    )
  )
}

print(results)

# ---------------------------
# MEAN ± SD PLOT PER CYLINDER
# ---------------------------
plot_data <- df_wide %>%
  select(Cylinder_mm, Accuracy, SkinStretchEffect) %>%
  pivot_longer(
    cols = c(Accuracy, SkinStretchEffect),
    names_to = "Measure",
    values_to = "Angle_deg"
  ) %>%
  group_by(Cylinder_mm, Measure) %>%
  summarise(
    Mean_deg = mean(Angle_deg, na.rm = TRUE),
    SD_deg   = sd(Angle_deg, na.rm = TRUE),
    .groups = "drop"
  )

p_cylinder <- ggplot(
  plot_data,
  aes(
    x = Cylinder_mm,
    y = Mean_deg,
    fill = Measure
  )
) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_errorbar(
    aes(
      ymin = Mean_deg - SD_deg,
      ymax = Mean_deg + SD_deg
    ),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(
    values = c(
      Accuracy    = "#3BBFB3", #"#2ca02c",
      SkinStretchEffect = "#D62728" #"#9467bd"
    ),
    labels = c(
      Accuracy =
        "Hand matching error\nBaseline_LeapMotion − BareHand",
      SkinStretchEffect =
        "Skin stretch effect\nBaseline_Encoder − DevOn"
    )
  ) +
  labs(
    x = "Cylinder diameter [mm]",
    y = "Variation [deg]",
    fill = "",
    title = "Hand matching error vs skin stretch effect per cylinder"
  ) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "top")

print(p_cylinder)
