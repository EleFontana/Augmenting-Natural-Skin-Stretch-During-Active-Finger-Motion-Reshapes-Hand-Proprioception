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
# Extracts the mean value of a signal plateau around its peak
# ---------------------------
extract_plateau <- function(wave, plateau_prop = 0.90, smooth_k = 7,
                            min_len = 3, fallback_win = 11) {
  
  if (is.null(wave)) return(list(mean = NA, i_start = NA, i_end = NA, len = 0))
  
  wave_num <- as.numeric(wave)
  wave_num <- wave_num[!is.na(wave_num)]
  n <- length(wave_num)
  if (n == 0) return(list(mean = NA, i_start = NA, i_end = NA, len = 0))
  
  # Moving-average smoothing
  k <- min(smooth_k, n)
  if (k > 1) {
    filt <- rep(1 / k, k)
    sm <- stats::filter(wave_num, filt, sides = 2)
    sm <- as.numeric(sm)
    sm[is.na(sm)] <- wave_num[is.na(sm)]
  } else {
    sm <- wave_num
  }
  
  # Identify peak and define plateau threshold
  i_max <- which.max(sm)
  thr <- plateau_prop * sm[i_max]
  
  # Expand plateau around the peak
  i_left <- i_max
  while (i_left > 1 && sm[i_left - 1] >= thr) i_left <- i_left - 1
  i_right <- i_max
  while (i_right < n && sm[i_right + 1] >= thr) i_right <- i_right + 1
  
  plateau_idx <- i_left:i_right
  
  # Fallback window if plateau is too short
  if (length(plateau_idx) < min_len) {
    half <- floor(fallback_win / 2)
    i_start <- max(1, i_max - half)
    i_end <- min(n, i_max + half)
    plateau_idx <- i_start:i_end
  }
  
  plateau_mean <- mean(wave_num[plateau_idx], na.rm = TRUE)
  
  list(mean = plateau_mean,
       i_start = min(plateau_idx),
       i_end = max(plateau_idx),
       len = length(plateau_idx))
}

# ---------------------------
# FUNCTION: robust numeric extraction from MATLAB structures
# ---------------------------
get_numeric_vector <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  if (is.matrix(x)) return(as.numeric(x))
  if (is.list(x) && length(x) == 1) return(get_numeric_vector(x[[1]]))
  if (is.list(x)) return(as.numeric(unlist(x)))
  numeric(0)
}

# ---------------------------
# PARAMETERS
# ---------------------------
plateau_prop <- 0.90
smooth_k <- 11
min_len <- 75
fallback_win <- 75

diameters <- c(45, 60, 75, 90)   # Real cylinder diameters [mm]
n_subj <- 12

# Specific parameters for Baseline_Encoder
baseline_encoder_params <- list(
  plateau_prop = 0.90,
  smooth_k = 11, 
  min_len = 20,
  fallback_win = 20
)

# ---------------------------
# INITIALIZE DATAFRAME
# ---------------------------
all_data <- tibble(
  Subject = integer(),
  Condition = character(),
  Cylinder = integer(),
  PlateauAngle = numeric(),
  i_start = integer(),
  i_end = integer(),
  plateau_len = integer()
)

# ---------------------------
# LOAD AND PROCESS DATA
# Baseline conditions are extracted but excluded later from analysis and plots
# ---------------------------
for (i in 1:n_subj) {
  
  fname <- paste0("sub", i, ".mat")
  if (!file.exists(fname)) {
    warning("File not found: ", fname, " -> skip")
    next
  }
  
  mat <- readMat(fname)
  varname <- paste0("sub", i)
  
  if (!varname %in% names(mat)) {
    warning("Variable ", varname, " not found in ", fname, " -> skipped")
    next
  }
  
  subject_cell <- mat[[varname]]  # 5 conditions × 4 cylinders
  
  for (r in 1:5) {
    for (c in 1:4) {
      
      raw <- subject_cell[[r, c]]
      wave_vals <- get_numeric_vector(raw)
      
      # Use specific parameters for Baseline_Encoder
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
      
      res <- if (length(wave_vals) == 0) {
        list(mean = NA, i_start = NA, i_end = NA, len = 0)
      } else {
        extract_plateau(wave_vals, pp, sk, ml, fw)
      }
      
      all_data <- bind_rows(
        all_data,
        tibble(
          Subject = i,
          Condition = c("DevOff", "DevOn", "BareHand",
                        "Baseline_Encoder", "Baseline_LeapMotion")[r],
          Cylinder = c,
          PlateauAngle = res$mean,
          i_start = res$i_start,
          i_end = res$i_end,
          plateau_len = res$len
        )
      )
    }
  }
}

# ---------------------------
# MANUAL CORRECTIONS (KNOWN OUTLIERS)
# ---------------------------
all_data <- all_data %>%
  mutate(
    PlateauAngle = case_when(
      Subject == 8  & Cylinder == 1 & Condition == "Baseline_LeapMotion" ~ 52.8075,
      Subject == 12 & Cylinder == 4 & Condition == "Baseline_LeapMotion" ~ 31.7611,
      TRUE ~ PlateauAngle
    )
  )

# ---------------------------
# FILTER DATA FOR ANALYSIS AND PLOTTING
# Baseline conditions are explicitly excluded here
# ---------------------------
all_data_first <- all_data %>%
  filter(
    !is.na(PlateauAngle),
    Condition %in% c("DevOff", "BareHand", "DevOn")
  ) %>%
  mutate(
    Subject = factor(Subject),
    Condition = factor(Condition, levels = c("DevOff", "BareHand", "DevOn")),
    Cylinder = diameters[Cylinder]
  )

# ---------------------------
# LINEAR MIXED MODEL
# ---------------------------
mod <- lmer(
  PlateauAngle ~ Condition + Cylinder + (1 | Subject),
  data = all_data_first,
  REML = FALSE
)
summary(mod)

# Pairwise comparisons between conditions
emm_res <- emmeans(mod, pairwise ~ Condition)
print(emm_res$contrasts)

# ---------------------------
# PLOTS
# ---------------------------

# Boxplot by condition (only DevOff, BareHand, DevOn)
p1 <- ggplot(all_data_first,
             aes(x = Condition, y = PlateauAngle, color = Cylinder)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, size = 1.5) +
  geom_jitter(width = 0.15, size = 4, alpha = 0.9) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 22) +
  labs(
    y = "PIP Angle [deg]",
    x = "Condition",
    color = "Cylinder diameter [mm]"
  ) +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22)
  )
print(p1)

# Comparison between conditions for each cylinder
p2 <- ggplot(all_data_first,
             aes(x = factor(Cylinder), y = PlateauAngle, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar",
               position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(
    aes(color = Condition),
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
    size = 2, alpha = 0.8
  ) +
  scale_fill_manual(values = c("#63b8f3", "violet", "orange")) +
  scale_color_manual(values = c("#63b8f3", "violet", "orange")) +
  labs(
    x = "Cylinder diameter [mm]",
    y = "PIP Angle [deg]"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22)
  )
print(p2)


# ---------------------------
# Mean and STD
# ---------------------------
# ---------------------------
# Descriptive statistics per condition
# Mean ± SE and Median + IQR
# ---------------------------
summary_all_data <- all_data_first %>%
  group_by(Condition) %>%
  summarise(
    n = sum(!is.na(PlateauAngle)),
    
    Mean = mean(PlateauAngle, na.rm = TRUE),
    SD   = sd(PlateauAngle, na.rm = TRUE),
    SE   = SD / sqrt(n),
    
    Median = median(PlateauAngle, na.rm = TRUE),
    Q1     = quantile(PlateauAngle, 0.25, na.rm = TRUE),
    Q3     = quantile(PlateauAngle, 0.75, na.rm = TRUE),
    IQR    = Q3 - Q1,
    
    .groups = "drop"
  )

print(summary_all_data)


