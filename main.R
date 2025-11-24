library(dplyr)
library(readr)
library(tidyverse)
library(monomvn)
library(BoomSpikeSlab)
library(car)
library(MASS)
# ===============================================
#                 DATA PREPROCESS
# ===============================================
FOOD_DATA_GROUP1 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP1.csv")
FOOD_DATA_GROUP2 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP2.csv")
FOOD_DATA_GROUP3 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP3.csv")
FOOD_DATA_GROUP4 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP4.csv")
FOOD_DATA_GROUP5 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP5.csv")

FOOD_DATA <- bind_rows(
  dplyr::mutate(FOOD_DATA_GROUP1, group = "group1"),
  dplyr::mutate(FOOD_DATA_GROUP2, group = "group2"),
  dplyr::mutate(FOOD_DATA_GROUP3, group = "group3"),
  dplyr::mutate(FOOD_DATA_GROUP4, group = "group4"),
  dplyr::mutate(FOOD_DATA_GROUP5, group = "group5")
) %>%
  # remove first two columns
  dplyr::select(-1, -2) %>%
  # replace spaces with underscores in colnames
  dplyr::rename_with(~ gsub(" ", "_", .x))

df <- FOOD_DATA %>%
  dplyr::select(-food, -group, -Nutrition_Density) %>%
  drop_na()

y <- df$Caloric_Value
X <- df %>% dplyr::select(-Caloric_Value)

X_scaled <- scale(X)

# Fat multicollinearity
plot(unlist(FOOD_DATA[, "Saturated_Fats"] + FOOD_DATA[, "Monounsaturated_Fats"] + FOOD_DATA[, "Polyunsaturated_Fats"]),
     unlist(FOOD_DATA[, "Fat"]))
abline(0, 1)


# ===============================================
#                 EDA PLOTS
# ===============================================
# Density Plot for the Response
df %>% 
  filter(Caloric_Value <= quantile(Caloric_Value, 0.975)) %>% 
  ggplot(aes(x = Caloric_Value)) +
  geom_density(alpha = 0.5, fill = "skyblue") +
  labs(
    x = "Caloric Value",
    y = "Density"
  )

# Density Plots for some regressors
some_regressors <- c(
  "Fat", "Carbohydrates", "Protein",
  "Saturated_Fats", "Monounsaturated_Fats", "Polyunsaturated_Fats",
  "Sugars", "Vitamin_A", "Vitamin_E",
  "Water", "Cholesterol", "Magnesium")

df_long <- df[c(some_regressors, "Caloric_Value")] %>%
  filter(if_all(all_of(some_regressors), ~ . <= quantile(., 0.975))) %>% # . >= quantile(., 0.025) &
  pivot_longer(
    cols = all_of(some_regressors), # or cols = c(A, B, C)
    names_to = "variable",
    values_to = "value"
  )

ggplot(df_long, aes(x = value)) +
  geom_density(alpha = 0.5, fill = "skyblue") +
  facet_wrap(~ variable, scales = "free") + # Creates separate plots for each variable
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Regressor Value",
    y = "Density"
  )

# Correlation Plot for some regressors
library(corrplot)

M = cor(df[c("Caloric_Value", "Fat", "Saturated_Fats", "Monounsaturated_Fats",
             "Polyunsaturated_Fats", "Carbohydrates", "Sugars", "Protein")])
colnames(M) <- c("Caloric Value", "Fats", "Saturated Fats", "Monounsaturated Fats",
                 "Polyunsaturated Fats", "Carbohydrates", "Sugars", "Protein")
rownames(M) <- colnames(M)
corrplot.mixed(M, diag="l", tl.pos = "lt", tl.col = "#000011", tl.srt = 45)

# Regressors Vs Response
num_labels <- c(
  "Caloric_Value" = "Caloric Value",
  "Saturated_Fats" = "Saturated Fats",
  "Monounsaturated_Fats" = "Monounsaturated Fats",
  "Polyunsaturated_Fats" = "Polyunsaturated Fats",
  "Fat" = "Fat",
  "Carbohydrates" = "Carbohydrates",
  "Sugars" = "Sugars",
  "Protein" = "Protein"
)

main_regressors <- c(
  "Fat", "Carbohydrates", "Protein",
  "Saturated_Fats", "Monounsaturated_Fats", "Polyunsaturated_Fats")

df_long <- df[c(main_regressors, "Caloric_Value")] %>%
  pivot_longer(
    cols = all_of(main_regressors),
    names_to = "regressor_variable",
    values_to = "regressor_value"
  )


ggplot(df_long, aes(x = regressor_value, y = Caloric_Value)) +
  geom_point(alpha = 0.5, color = "skyblue") +
  facet_wrap(~ factor(regressor_variable, 
                      levels = main_regressors), 
             scales = "free",
             labeller = as_labeller(num_labels)) +
  labs(
    x = "Regressor Value",
    y = "Caloric Value (Response)") +
  theme_minimal()


# =========================================================
#                        BAYESIAN LASSO
# ---------------------------------------------------------
#                          With FAT
# =========================================================
set.seed(123)
bl_fit <- blasso(
  X = X_scaled,
  y = y,
  T = 10000,
  thin = 5,
  verb = 0
)

# Posterior mean of each coefficient
beta_mean <- apply(bl_fit$beta, 2, mean)
names(beta_mean) <- colnames(X_scaled)

# Rank nutrients by absolute effect size
ranked_lasso <- sort(abs(beta_mean), decreasing = TRUE)
head(ranked_lasso, 20)

beta_ci <- t(apply(bl_fit$beta, 2, quantile, probs = c(0.025, 0.975)))
colnames(beta_ci) <- c("low", "high")

# Selected nutrients from Bayesian Lasso
selected_lasso_ci <- rownames(beta_ci)[beta_ci[, "low"] > 0 | beta_ci[, "high"] < 0]
selected_lasso_ci # The Bayesian lasso thinks those predictors have credible, non-zero effects on "Caloric Value".

idx_sel <- as.integer(sub("b\\.", "", selected_lasso_ci))
idx_sel 

selected_nutrients <- colnames(X_scaled)[idx_sel]
selected_nutrients



# =========================================================
#                        BAYESIAN LASSO
# ---------------------------------------------------------
#                          With FAT
# =========================================================
X_without_fat <- X %>% dplyr::select(-Fat)

X_scaled_without_fat <- scale(X_without_fat)

# Bayesian Lasso without fat
set.seed(123)
bl_fit_without_fat <- blasso(
  X = X_scaled_without_fat,
  y = y,
  T = 10000,
  thin = 5,
  verb = 0
)

# Posterior mean of each coefficient
beta_mean_without_fat <- apply(bl_fit_without_fat$beta, 2, mean)
names(beta_mean_without_fat) <- colnames(X_scaled_without_fat)

# Rank nutrients by absolute effect size
ranked_lasso_without_fat <- sort(abs(beta_mean_without_fat), decreasing = TRUE)
head(ranked_lasso_without_fat, 20)

beta_ci_without_fat <- t(apply(bl_fit_without_fat$beta, 2, quantile, probs = c(0.025, 0.975)))
colnames(beta_ci_without_fat) <- c("low", "high")

# Selected nutrients from Bayesian Lasso
selected_lasso_ci_without_fat <- rownames(beta_ci_without_fat)[beta_ci_without_fat[, "low"] > 0 | beta_ci_without_fat[, "high"] < 0]
selected_lasso_ci_without_fat # The Bayesian lasso thinks those predictors have credible, non-zero effects on "Caloric Value".

idx_sel_without_fat <- as.integer(sub("b\\.", "", selected_lasso_ci_without_fat))
idx_sel_without_fat 

selected_nutrients_without_fat <- colnames(X_scaled_without_fat)[idx_sel_without_fat]
selected_nutrients_without_fat


# ===============================================
#  Comparison of Lasso: With FAT vs. Without FAT
# ===============================================
lasso_with_fat <- lm(Caloric_Value ~ Fat + Saturated_Fats + Monounsaturated_Fats + Polyunsaturated_Fats +
                       Carbohydrates + Protein + Water + Vitamin_A + Magnesium + Phosphorus,
                     data = FOOD_DATA)

lasso_without_fat <- lm(Caloric_Value ~ Vitamin_E + Saturated_Fats + Monounsaturated_Fats + Polyunsaturated_Fats +
                          Carbohydrates + Protein + Cholesterol + Vitamin_A + Magnesium + Phosphorus,
                        data = FOOD_DATA)

vif(lasso_with_fat)
vif(lasso_without_fat)

summary(lasso_with_fat)
summary(lasso_without_fat)

# Investigating foods where sum of fats greater than fat
fat_check_df <- FOOD_DATA %>% 
  filter(Saturated_Fats + Monounsaturated_Fats + Polyunsaturated_Fats > Fat)

# ===============================================
#                 SPIKE AND SLAB PRIOR
# ===============================================
set.seed(123)
df_model <- data.frame(
  y = as.numeric(y),
  as.data.frame(X_scaled)
)

ss_fit <- lm.spike(
  y ~ .,
  data  = df_model,
  niter = 6000,                 
  expected.model.size = 10    # prior belief of how many predictors matter
)

burn_in <- 1000

ss_sum <- summary(ss_fit, burn = burn_in)

# matrix with columns:
# 1 post mean, 2 post sd,
# 3 mean|nonzero, 4 sd|nonzero,
# 5 Pr(nonzero)
coef_table <- ss_sum$coefficients

inc_prob <- coef_table[, 5]
inc_prob <- inc_prob[names(inc_prob) != "(Intercept)"]

inc_prob_sorted <- sort(inc_prob, decreasing = TRUE)
selected_vars <- names(inc_prob)[inc_prob > 0.5] # keep any nutrient whose coefficient is nonzero in more than half the posterior samples.
selected_vars


# ===============================================
#                 BAYESIAN LINEAR REGRESSION
# ===============================================
# Doing Bayesian linear regression with all predictors
# Using a conjugate Normal-Inv-Gamma prior posterior
# prior hyperparameters (noninformative)

# Model 1: Full model with all variables as regressors
n_full <- nrow(X_scaled)
p_full <- ncol(X_scaled)

xi_full <- rep(0, p_full)
Omega_full <- 10^3 * diag(p_full)
alpha_full <- 0.01
b_full <- 0.01  # using b instead of beta for rate parameter

Omega_inverse_full <- solve(Omega_full)
Q_beta_full <- t(as.matrix(X_scaled)) %*% X_scaled + Omega_inverse_full
Q_beta_inverse_full <- solve(Q_beta_full)
l_beta_full <- t(as.matrix(X_scaled)) %*% y + Omega_inverse_full %*% xi_full

# Sampling posterior
samples_full <- 10^4
sigma_2_samples_full <- rinvgamma(
  samples_full,
  alpha_full + n_full / 2,
  b_full + t(y) %*% y / 2 +
    t(xi_full) %*% Omega_inverse_full %*% xi_full / 2 -
    t(l_beta_full) %*% Q_beta_inverse_full %*% l_beta_full / 2
)

beta_samples_full <- matrix(NA, nrow = samples_full, ncol = p_full)
for (s in 1:samples_full) {
  beta_samples_full[s, ] <- mvrnorm(
    1,
    Q_beta_inverse_full %*% l_beta_full,
    Sigma = sigma_2_samples_full[s] * Q_beta_inverse_full
  )
}

hist(sigma_2_samples_full, xlab = expression(sigma^2))

# Credible intervals of betas
beta_summary_full <- data.frame(
  term = c("intercept", colnames(X[, -1])),
  mean = colMeans(beta_samples_full),
  low = apply(beta_samples_full, 2, quantile, 0.025),
  high = apply(beta_samples_full, 2, quantile, 0.975)
)
beta_summary_full

selected_betas_full <- beta_summary_full[
  beta_summary_full[, "low"] > 0 | beta_summary_full[, "high"] < 0,
]
selected_betas_full


# Model 2: Bayesian Lasso Model with Fat included initially
X_scaled_BL_F <- X_scaled[, selected_nutrients, drop = FALSE]

n_BL_F <- nrow(X_scaled_BL_F)
p_BL_F <- ncol(X_scaled_BL_F)

# Prior settings
xi_BL_F <- rep(0, p_BL_F)
Omega_BL_F <- 10^3 * diag(p_BL_F)
alpha_BL_F <- 0.01
b_BL_F <- 0.01  # rate parameter

Omega_inverse_BL_F <- solve(Omega_BL_F)

Q_beta_BL_F <- t(as.matrix(X_scaled_BL_F)) %*% X_scaled_BL_F + Omega_inverse_BL_F
Q_beta_inverse_BL_F <- solve(Q_beta_BL_F)

l_beta_BL_F <- t(as.matrix(X_scaled_BL_F)) %*% y + Omega_inverse_BL_F %*% xi_BL_F

# Posterior sampling
samples_BL_F <- 10^4

sigma_2_samples_BL_F <- rinvgamma(
  samples_BL_F,
  alpha_BL_F + n_BL_F / 2,
  b_BL_F + t(y) %*% y / 2 +
    t(xi_BL_F) %*% Omega_inverse_BL_F %*% xi_BL_F / 2 -
    t(l_beta_BL_F) %*% Q_beta_inverse_BL_F %*% l_beta_BL_F / 2
)

beta_samples_BL_F <- matrix(NA, nrow = samples_BL_F, ncol = p_BL_F)

for (s in 1:samples_BL_F) {
  beta_samples_BL_F[s, ] <- mvrnorm(
    1,
    Q_beta_inverse_BL_F %*% l_beta_BL_F,
    Sigma = sigma_2_samples_BL_F[s] * Q_beta_inverse_BL_F
  )
}

hist(sigma_2_samples_BL_F, xlab = expression(sigma^2))

# Credible intervals
beta_summary_BL_F <- data.frame(
  term = colnames(X_scaled_BL_F),
  mean = colMeans(beta_samples_BL_F),
  low = apply(beta_samples_BL_F, 2, quantile, 0.025),
  high = apply(beta_samples_BL_F, 2, quantile, 0.975)
)

beta_summary_BL_F

selected_betas_BL_F <- beta_summary_BL_F[
  beta_summary_BL_F$low > 0 | beta_summary_BL_F$high < 0,
]

selected_betas_BL_F

# Model 3: Bayesian Lasso Model without Fat included initially
X_scaled_BL <- X_scaled[, selected_nutrients_without_fat, drop = FALSE]

n_BL <- nrow(X_scaled_BL)
p_BL <- ncol(X_scaled_BL)

xi_BL <- rep(0, p_BL)
Omega_BL <- 10^3 * diag(p_BL)
alpha_BL <- 0.01
b_BL <- 0.01  # rate parameter

Omega_inverse_BL <- solve(Omega_BL)

Q_beta_BL <- t(as.matrix(X_scaled_BL)) %*% X_scaled_BL + Omega_inverse_BL
Q_beta_inverse_BL <- solve(Q_beta_BL)

l_beta_BL <- t(as.matrix(X_scaled_BL)) %*% y + Omega_inverse_BL %*% xi_BL

# Posterior sampling
samples_BL <- 10^4

sigma_2_samples_BL <- rinvgamma(
  samples_BL,
  alpha_BL + n_BL / 2,
  b_BL + t(y) %*% y / 2 +
    t(xi_BL) %*% Omega_inverse_BL %*% xi_BL / 2 -
    t(l_beta_BL) %*% Q_beta_inverse_BL %*% l_beta_BL / 2
)

beta_samples_BL <- matrix(NA, nrow = samples_BL, ncol = p_BL)

for (s in 1:samples_BL) {
  beta_samples_BL[s, ] <- mvrnorm(
    1,
    Q_beta_inverse_BL %*% l_beta_BL,
    Sigma = sigma_2_samples_BL[s] * Q_beta_inverse_BL
  )
}

# Plot sigma^2 samples
hist(sigma_2_samples_BL, xlab = expression(sigma^2))

# Credible intervals
beta_summary_BL <- data.frame(
  term = colnames(X_scaled_BL),
  mean = colMeans(beta_samples_BL),
  low = apply(beta_samples_BL, 2, quantile, 0.025),
  high = apply(beta_samples_BL, 2, quantile, 0.975)
)

beta_summary_BL

# Select coefficients with credible intervals not covering 0
selected_betas_BL <- beta_summary_BL[
  beta_summary_BL$low > 0 | beta_summary_BL$high < 0,
]

selected_betas_BL


# Model 4: Spike and Slab Prior
X_scaled_SS <- X_scaled[, selected_vars, drop = FALSE]
# If spike & slab chose nothing, stop early
if (length(selected_vars) == 0) {
  stop("Spike&Slab selected no variables (inc_prob > 0.5). Try lowering threshold.")
}

n_SS <- nrow(X_scaled_SS)
p_SS <- ncol(X_scaled_SS)

# Prior hyperparameters (same weakly-informative N-Inv-Gamma setup)
xi_SS <- rep(0, p_SS)
Omega_SS <- 10^3 * diag(p_SS)
alpha_SS <- 0.01
b_SS <- 0.01

Omega_inverse_SS <- solve(Omega_SS)

# Posterior precision and mean terms
Q_beta_SS <- t(as.matrix(X_scaled_SS)) %*% X_scaled_SS + Omega_inverse_SS
Q_beta_inverse_SS <- solve(Q_beta_SS)

l_beta_SS <- t(as.matrix(X_scaled_SS)) %*% y + Omega_inverse_SS %*% xi_SS

# Posterior sampling
samples_SS <- 10^4

sigma_2_samples_SS <- rinvgamma(
  samples_SS,
  alpha_SS + n_SS / 2,
  b_SS + t(y) %*% y / 2 +
    t(xi_SS) %*% Omega_inverse_SS %*% xi_SS / 2 -
    t(l_beta_SS) %*% Q_beta_inverse_SS %*% l_beta_SS / 2
)

beta_samples_SS <- matrix(NA, nrow = samples_SS, ncol = p_SS)

for (s in 1:samples_SS) {
  beta_samples_SS[s, ] <- mvrnorm(
    1,
    Q_beta_inverse_SS %*% l_beta_SS,
    Sigma = sigma_2_samples_SS[s] * Q_beta_inverse_SS
  )
}

# Plot sigma^2 posterior samples
hist(sigma_2_samples_SS, xlab = expression(sigma^2),
     main = "Posterior draws of sigma^2 (Spike&Slab Model)")

# Credible intervals for betas
beta_summary_SS <- data.frame(
  term = colnames(X_scaled_SS),
  mean = colMeans(beta_samples_SS),
  low = apply(beta_samples_SS, 2, quantile, 0.025),
  high = apply(beta_samples_SS, 2, quantile, 0.975)
)

beta_summary_SS

# Select coefficients whose 95% CI excludes 0
selected_betas_SS <- beta_summary_SS[
  beta_summary_SS$low > 0 | beta_summary_SS$high < 0,
]

selected_betas_SS


# ===============================================
#   TRAIN / TEST SPLIT (shared by all models)
# ===============================================
set.seed(123)

n <- nrow(df)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- X[train_idx, , drop = FALSE]
X_test  <- X[test_idx,  , drop = FALSE]
y_train <- y[train_idx]
y_test  <- y[test_idx]

# ---- scale using TRAIN stats only
train_means <- apply(X_train, 2, mean)
train_sds   <- apply(X_train, 2, sd)

scale_with_train <- function(Xmat, mu, sdv) {
  sweep(sweep(as.matrix(Xmat), 2, mu, "-"), 2, sdv, "/")
}

X_train_scaled <- scale_with_train(X_train, train_means, train_sds)
X_test_scaled  <- scale_with_train(X_test,  train_means, train_sds)

# ===============================================
#  Helper: Conjugate Bayesian LR + Test R^2
# ===============================================
fit_conjugate_blr_r2 <- function(Xtr_sc, ytr, Xte_sc, yte,
                                 samples = 1e4,
                                 alpha0 = 0.01, b0 = 0.01,
                                 Omega_scale = 1e3) {
  
  n_tr <- nrow(Xtr_sc)
  p_tr <- ncol(Xtr_sc)
  
  xi <- rep(0, p_tr)
  Omega <- Omega_scale * diag(p_tr)
  Omega_inv <- solve(Omega)
  
  Q_beta <- t(Xtr_sc) %*% Xtr_sc + Omega_inv
  Q_beta_inv <- solve(Q_beta)
  l_beta <- t(Xtr_sc) %*% ytr + Omega_inv %*% xi
  
  # posterior draws
  sigma2_draws <- rinvgamma(
    samples,
    alpha0 + n_tr / 2,
    b0 + t(ytr) %*% ytr / 2 +
      t(xi) %*% Omega_inv %*% xi / 2 -
      t(l_beta) %*% Q_beta_inv %*% l_beta / 2
  )
  
  beta_draws <- matrix(NA, nrow = samples, ncol = p_tr)
  beta_mean_post <- as.vector(Q_beta_inv %*% l_beta)
  
  for (s in 1:samples) {
    beta_draws[s, ] <- mvrnorm(
      1, mu = beta_mean_post,
      Sigma = sigma2_draws[s] * Q_beta_inv
    )
  }
  
  beta_post_mean <- colMeans(beta_draws)
  
  # posterior mean predictions on test
  yhat_test <- as.vector(Xte_sc %*% beta_post_mean)
  
  # test R^2
  sse <- sum((yte - yhat_test)^2)
  sst <- sum((yte - mean(yte))^2)
  r2  <- 1 - sse / sst
  
  list(r2 = r2, beta_post_mean = beta_post_mean)
}

# ===============================================
# Model 1: Full model (all predictors)
# ===============================================
m1 <- fit_conjugate_blr_r2(X_train_scaled, y_train,
                           X_test_scaled,  y_test)

# ===============================================
# Model 2: Bayesian Lasso selection (with Fat)
#   selection on TRAIN only
# ===============================================
set.seed(123)
bl_fit_tr <- blasso(
  X = X_train_scaled,
  y = y_train,
  T = 10000,
  thin = 5,
  verb = 0
)

beta_ci_tr <- t(apply(bl_fit_tr$beta, 2, quantile, probs = c(0.025, 0.975)))
colnames(beta_ci_tr) <- c("low", "high")

selected_lasso_ci_tr <- rownames(beta_ci_tr)[
  beta_ci_tr[, "low"] > 0 | beta_ci_tr[, "high"] < 0
]

idx_sel_tr <- as.integer(sub("b\\.", "", selected_lasso_ci_tr))
selected_nutrients_tr <- colnames(X_train_scaled)[idx_sel_tr]

# subset train/test using the SAME columns
Xtr_BL_F <- X_train_scaled[, selected_nutrients_tr, drop = FALSE]
Xte_BL_F <- X_test_scaled[,  selected_nutrients_tr, drop = FALSE]

m2 <- fit_conjugate_blr_r2(Xtr_BL_F, y_train,
                           Xte_BL_F, y_test)

# ===============================================
# Model 3: Bayesian Lasso selection (without Fat)
#   selection on TRAIN only
# ===============================================
# remove Fat BEFORE scaling to keep correct alignment
X_train_wofat <- X_train %>% dplyr::select(-Fat)
X_test_wofat  <- X_test  %>% dplyr::select(-Fat)

mu_wofat <- apply(X_train_wofat, 2, mean)
sd_wofat <- apply(X_train_wofat, 2, sd)

Xtr_wofat_sc <- scale_with_train(X_train_wofat, mu_wofat, sd_wofat)
Xte_wofat_sc <- scale_with_train(X_test_wofat,  mu_wofat, sd_wofat)

set.seed(123)
bl_fit_wofat_tr <- blasso(
  X = Xtr_wofat_sc,
  y = y_train,
  T = 10000,
  thin = 5,
  verb = 0
)

beta_ci_wofat_tr <- t(apply(bl_fit_wofat_tr$beta, 2, quantile, probs = c(0.025, 0.975)))
colnames(beta_ci_wofat_tr) <- c("low", "high")

selected_lasso_ci_wofat_tr <- rownames(beta_ci_wofat_tr)[
  beta_ci_wofat_tr[, "low"] > 0 | beta_ci_wofat_tr[, "high"] < 0
]

idx_sel_wofat_tr <- as.integer(sub("b\\.", "", selected_lasso_ci_wofat_tr))
selected_nutrients_wofat_tr <- colnames(Xtr_wofat_sc)[idx_sel_wofat_tr]

Xtr_BL <- Xtr_wofat_sc[, selected_nutrients_wofat_tr, drop = FALSE]
Xte_BL <- Xte_wofat_sc[, selected_nutrients_wofat_tr, drop = FALSE]

m3 <- fit_conjugate_blr_r2(Xtr_BL, y_train,
                           Xte_BL, y_test)

# ===============================================
# Model 4: Spike & Slab selection
#   selection on TRAIN only
# ===============================================
set.seed(123)
df_model_tr <- data.frame(
  y = as.numeric(y_train),
  as.data.frame(X_train_scaled)
)

ss_fit_tr <- lm.spike(
  y ~ .,
  data  = df_model_tr,
  niter = 6000,
  expected.model.size = 10
)

burn_in <- 1000
ss_sum_tr <- summary(ss_fit_tr, burn = burn_in)
coef_table_tr <- ss_sum_tr$coefficients

inc_prob_tr <- coef_table_tr[, 5]
inc_prob_tr <- inc_prob_tr[names(inc_prob_tr) != "(Intercept)"]

selected_vars_tr <- names(inc_prob_tr)[inc_prob_tr > 0.5]

if (length(selected_vars_tr) == 0) {
  stop("Spike&Slab selected no variables on TRAIN (inc_prob > 0.5). Lower threshold.")
}

Xtr_SS <- X_train_scaled[, selected_vars_tr, drop = FALSE]
Xte_SS <- X_test_scaled[,  selected_vars_tr, drop = FALSE]

m4 <- fit_conjugate_blr_r2(Xtr_SS, y_train,
                           Xte_SS, y_test)

# ===============================================
# Compare Test R^2 Across Models
# ===============================================
r2_table <- tibble(
  Model = c("Model 1: Full", 
            "Model 2: BLasso (with Fat)", 
            "Model 3: BLasso (no Fat)", 
            "Model 4: Spike&Slab"),
  Test_R2 = c(m1$r2, m2$r2, m3$r2, m4$r2),
  Num_Predictors = c(
    ncol(X_train_scaled),
    ncol(Xtr_BL_F),
    ncol(Xtr_BL),
    ncol(Xtr_SS)
  )
)

print(r2_table)


# ===============================================
# Model 5: Only Fat + Carbohydrates + Protein
# ===============================================

vars_m5 <- c("Fat", "Carbohydrates", "Protein")

# make sure they exist
stopifnot(all(vars_m5 %in% colnames(X_train)))

# subset original (unscaled) train/test
X_train_m5 <- X_train[, vars_m5, drop = FALSE]
X_test_m5  <- X_test[,  vars_m5, drop = FALSE]

# scale using TRAIN stats of these 3 variables
mu_m5 <- apply(X_train_m5, 2, mean)
sd_m5 <- apply(X_train_m5, 2, sd)

Xtr_m5_sc <- scale_with_train(X_train_m5, mu_m5, sd_m5)
Xte_m5_sc <- scale_with_train(X_test_m5,  mu_m5, sd_m5)

m5 <- fit_conjugate_blr_r2(Xtr_m5_sc, y_train,
                           Xte_m5_sc, y_test)

# add to comparison table
r2_table <- r2_table %>%
  add_row(
    Model = "Model 5: Fat+Carbs+Protein only",
    Test_R2 = m5$r2,
    Num_Predictors = ncol(Xtr_m5_sc)
  )

print(r2_table)

