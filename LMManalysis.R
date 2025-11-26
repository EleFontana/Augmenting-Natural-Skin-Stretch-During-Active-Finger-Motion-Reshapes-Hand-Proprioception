library(R.matlab)
library(dplyr)
library(lme4)
library(emmeans)
library(ggplot2)

# Extract plateau around peak
extract_plateau <- function(wave, plateau_prop = 0.90, smooth_k = 7,
                            min_len = 3, fallback_win = 11) {
  
  if (is.null(wave)) return(list(mean = NA, i_start = NA, i_end = NA, len = 0))
  
  wave_num <- as.numeric(unlist(wave))
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
  
  # If plateau is too short, fallback to a fixed window around the peak
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
# Parameters
plateau_prop <- 0.90  # Threshold relative to the peak 
smooth_k <- 11  # Smoothing window 
min_len <- 100  # Minimum samples to consider a plateau valid
fallback_win <- 100  # Samples around the peak
# ---------------------------

# Load all .mat files
n_subj <- 12
all_data <- tibble::tibble(
  Subject = integer(), Condition = character(),
  Cylinder = integer(), PlateauAngle = numeric(),
  i_start = integer(), i_end = integer(), plateau_len = integer()
)

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
  subject_cell <- mat[[varname]]  

  # Compute plateau for each cell
  for (r in 1:3) {
    for (c in 1:4) {
      raw <- subject_cell[[r, c]]
      if (is.null(raw) || length(raw) == 0) {
        res <- list(mean = NA, i_start = NA, i_end = NA, len = 0)
      } else {
        
        wave_vals <- as.numeric(unlist(raw))
        
        if (length(wave_vals) == 0 || all(is.na(wave_vals))) {
          res <- list(mean = NA, i_start = NA, i_end = NA, len = 0)
        } else {
          res <- extract_plateau(wave_vals,
                                 plateau_prop = plateau_prop,
                                 smooth_k = smooth_k,
                                 min_len = min_len,
                                 fallback_win = fallback_win)
        }
      }
      
      all_data <- dplyr::bind_rows(all_data,
                                   tibble::tibble(
                                     Subject = i,
                                     Condition = c("DevOff","BareHand","DevOn")[r],
                                     Cylinder = c,
                                     PlateauAngle = res$mean,
                                     i_start = res$i_start,
                                     i_end = res$i_end,
                                     plateau_len = res$len
                                   ))
    }
  }
}



# ---------------------------
# LMM on PlateauAngle (one row per trial)
diameters <- c(45, 60, 75, 90)

all_data <- all_data %>%
  mutate(
    Subject = factor(Subject),
    Condition = factor(Condition, levels = c("DevOff","BareHand","DevOn")),
    Cylinder = diameters[Cylinder]
  )

mod <- lmer(PlateauAngle ~ Condition + Cylinder + (1|Subject), data = all_data, REML = FALSE)
summary(mod)

# contrasts AvsB, AvsC, BvsC
emm_res <- emmeans(mod, pairwise ~ Condition)
print(emm_res$contrasts)


# ---------------------------
# PLOTS

# ============================
# Boxplot by Condition (color gradient by cylinder diameter)
# ============================

p1 <- ggplot(all_data, aes(x = Condition, y = PlateauAngle, color = Cylinder)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, size = 1.5) +   
  geom_jitter(width = 0.15, size = 4, alpha = 0.9) +             
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 18) +                                
  labs(y = "PIP Angle [deg]", x = "Condition", color = "Cylinder [mm]") +
  theme(
    axis.title = element_text(size = 26),         
    axis.text = element_text(size = 22),                         
    legend.title = element_text(size = 24),       
    legend.text = element_text(size = 22)                        
  )

print(p1)


# ============================
# Comparison between Conditions for each Cylinder
# ============================

bar_data <- all_data %>%
  group_by(Condition, Cylinder) %>%
  summarise(mean_plateau = mean(PlateauAngle, na.rm=TRUE),
            sd_plateau = sd(PlateauAngle, na.rm=TRUE),
            .groups="drop")

p2 <- ggplot(all_data, aes(x=factor(Cylinder), y=PlateauAngle, fill=Condition)) +
  stat_summary(fun=mean, geom="bar", position=position_dodge(width=0.8), alpha=0.7) +
  geom_jitter(aes(color=Condition),
              position=position_jitterdodge(jitter.width=0.1, dodge.width=0.8),
              size=2, alpha=0.8) +
  scale_fill_manual(values=c("#63b8f3","violet","orange"), name="Condition") +   
  scale_color_manual(values=c("#63b8f3","violet","orange"), name="Condition") +  
  labs(x="Cylinder diameter [mm]", y="PIP Angle [deg]") +               
  theme_minimal(base_size=14) +
  theme(
    axis.title = element_text(size=26),      
    axis.text = element_text(size=22),       
    legend.title = element_text(size=24),   
    legend.text = element_text(size=22)      
  )

print(p2)

