library(lme4)
library(dplyr)
library(ggplot2)

# ---- Core Parameters ----
set.seed(123)
n_sim <- 20  # For production-level power estimates
piglets_per_litter <- 12
baseline_prevalence <- 0.2
icc <- 0.1
residual_var <- pi^2 / 3  # Approx. variance for logistic distribution
between_litter_var <- (icc * residual_var) / (1 - icc)
between_litter_sd <- sqrt(between_litter_var)

# Function to simulate power at given OR and number of matched litter pairs
simulate_power_for_design <- function(n_pairs, log_odds_shift, n_sim = 20) {
  power_vec <- replicate(n_sim, {
    
    # Random effects per litter (2 litters per matched pair: 1 case, 1 control)
    litter_effects <- rnorm(n_pairs * 2, mean = 0, sd = between_litter_sd)
    
    df <- data.frame(
      PairID = rep(1:n_pairs, each = piglets_per_litter * 2),
      LitterID = 1:(n_pairs * 2),
      Group = rep(c("Control", "Case"), each = piglets_per_litter, times = n_pairs)
    )
    
    df <- df %>%
      arrange(LitterID) %>%
      mutate(
        PWD = ifelse(Group == "Case", 1, 0),
        LitterRE = rep(litter_effects, each = piglets_per_litter)
      )
    
    # Linear predictor
    lp <- log(baseline_prevalence / (1 - baseline_prevalence)) +
      df$PWD * log_odds_shift + df$LitterRE
    
    df$prob <- plogis(lp)
    df$E_coli <- rbinom(nrow(df), 1, df$prob)
    
    # Fit GLMM
    model <- tryCatch({
      glmer(
        E_coli ~ PWD + (1 | LitterID),  # PairID could be added too if needed
        data = df,
        family = binomial,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
      )
    }, error = function(e) NULL)
    
    # Extract p-value and return TRUE if significant
    if (!is.null(model)) {
      p_val <- tryCatch({
        summary(model)$coefficients["PWD", "Pr(>|z|)"]
      }, error = function(e) NA)
      return(!is.na(p_val) && p_val < 0.05)
    } else {
      return(FALSE)
    }
  })
  
  # Summary power results
  power <- mean(power_vec, na.rm = TRUE)
  se <- sqrt(power * (1 - power) / sum(!is.na(power_vec)))
  ci_low <- max(0, power - 1.96 * se)
  ci_high <- min(1, power + 1.96 * se)
  return(c(power = power, ci_low = ci_low, ci_high = ci_high))
}

# ---- Grid of Parameters ----
or_values <- seq(1.2, 2.0, by = 0.2)
log_odds_values <- log(or_values)
litter_range <- seq(10, 200, by = 20)  # Number of matched litter pairs

# Run simulations
results_list <- list()

for (or in or_values) {
  for (n in litter_range) {
    cat("Running simulation for OR =", or, "with", n, "matched pairs...\n")
    res <- simulate_power_for_design(n_pairs = n, log_odds_shift = log(or), n_sim = n_sim)
    results_list[[length(results_list) + 1]] <- data.frame(
      OR = or,
      MatchedLitterPairs = n,
      Power = res["power"],
      CI_Low = res["ci_low"],
      CI_High = res["ci_high"]
    )
  }
}

results_df <- bind_rows(results_list)

# ---- Plot Results ----
ggplot(results_df, aes(x = MatchedLitterPairs, y = Power, color = as.factor(OR))) +
  geom_ribbon(aes(ymin = CI_Low, ymax = CI_High, fill = as.factor(OR)), alpha = 0.1, color = NA) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray30") +
  scale_color_viridis_d(name = "Odds Ratio") +
  scale_fill_viridis_d(name = "Odds Ratio") +
  labs(
    x = "Number of Matched Litter Pairs",
    y = "Statistical Power",
    title = "Power Curve for Detecting PWD-Associated Microbial Features"
  ) +
  theme_minimal(base_size = 14)
