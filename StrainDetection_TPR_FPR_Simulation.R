# Required packages
# --- Load Libraries ---
install.packages(c("glmmTMB", "future.apply", "dplyr", "ggplot2", "tidyr", "pROC"))
library(glmmTMB)
library(future.apply)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pROC)

# --- Parallel Setup ---
plan(multisession, workers = parallel::detectCores())

# --- Parameters ---
set.seed(123)
n_sim <- 100  # reduce for testing, increase for final
piglets_per_litter <- 12
baseline_prevalence <- 0.166
baseline_log_odds <- log(baseline_prevalence / (1 - baseline_prevalence))

icc <- 0.1
residual_var <- pi^2 / 3
between_litter_var <- (icc * residual_var) / (1 - icc)
between_litter_sd <- sqrt(between_litter_var)

n_features <- 1000
n_true_positive <- 100
true_positive_idx <- sample(1:n_features, n_true_positive)
alpha <- 0.05

# --- Simulation Function ---
simulate_glmmTMB_power <- function(n_pairs, log_odds_shift_true, n_sim = 100) {
  future_lapply(1:n_sim, function(sim) {
    # Simulate piglet and litter data
    litter_effects <- rnorm(n_pairs * 2, mean = 0, sd = between_litter_sd)
    df <- data.frame(
      PairID = rep(1:n_pairs, each = piglets_per_litter * 2),
      LitterID = 1:(n_pairs * 2),
      Group = rep(c("Control", "Case"), each = piglets_per_litter, times = n_pairs)
    ) %>%
      mutate(
        PWD = ifelse(Group == "Case", 1, 0),
        LitterRE = rep(litter_effects, each = piglets_per_litter)
      )
    
    n_piglets <- nrow(df)
    Y <- matrix(0, nrow = n_piglets, ncol = n_features)
    
    # Simulate strain presence
    for (j in 1:n_features) {
      shift <- ifelse(j %in% true_positive_idx, log_odds_shift_true, 0)
      lp <- baseline_log_odds + df$PWD * shift + df$LitterRE
      prob <- plogis(lp)
      Y[, j] <- rbinom(n_piglets, 1, prob)
    }
    
    # Fit GLMM per feature
    pvals <- numeric(n_features)
    for (j in 1:n_features) {
      data_j <- df
      data_j$Feature <- Y[, j]
      
      model <- tryCatch({
        glmmTMB(Feature ~ PWD + (1 | LitterID), data = data_j, family = binomial)
      }, error = function(e) NULL)
      
      pvals[j] <- if (!is.null(model)) {
        p <- tryCatch(summary(model)$coefficients$cond["PWD", "Pr(>|z|)"], error = function(e) NA)
        ifelse(is.na(p), 1, p)
      } else {
        1
      }
    }
    
    # Calculate detection metrics
    detected <- which(pvals < alpha)
    tpr <- sum(detected %in% true_positive_idx) / n_true_positive
    fpr <- sum(!(detected %in% true_positive_idx)) / (n_features - n_true_positive)
    
    # ROC + AUC
    labels <- rep(0, n_features)
    labels[true_positive_idx] <- 1
    roc_obj <- roc(response = labels, predictor = -pvals, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    
    return(data.frame(Sim = sim, TPR = tpr, FPR = fpr, AUC = auc_val))
  }) %>% bind_rows()
}

# Simulate
results <- simulate_glmmTMB_power(
  n_pairs = 50,
  log_odds_shift_true = log(1.2),
  n_sim = n_sim
)

# Pivot longer for Metric type
results_long <- results %>%
  select(Sim, TPR, FPR) %>%
  pivot_longer(cols = c(TPR, FPR), names_to = "Metric", values_to = "Rate")

# Apply smoothing with a check
smoothed_df <- results_long %>%
  group_by(Metric) %>%
  arrange(Sim) %>%
  mutate(
    Smoothed = tryCatch(
      predict(loess(Rate ~ Sim, span = 0.2)),  # increased span
      error = function(e) Rate  # fallback: use original Rate if loess fails
    )
  )

# Plot
ggplot(smoothed_df, aes(x = Sim, y = Smoothed, fill = Metric)) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line(aes(color = Metric), size = 1.2) +
  scale_fill_manual(values = c(TPR = "#7A7C94", FPR = "#CC6677")) +
  scale_color_manual(values = c(TPR = "#7A7C94", FPR = "#CC6677")) +
  labs(
    title = "Smoothed TPR and FPR Across Simulations",
    x = "Simulation",
    y = "Smoothed Rate",
    fill = "Metric",
    color = "Metric"
  ) +
  theme_minimal(base_size = 14)
