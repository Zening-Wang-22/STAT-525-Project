library(dplyr)
library(readr)
library(tidyverse)
library(monomvn)
library(BoomSpikeSlab)
library(car)
library(MASS)


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


# Preprocess
df <- FOOD_DATA %>%
  dplyr::select(-food, -group, -Nutrition_Density) %>%
  drop_na()

y <- df$Caloric_Value
X <- df %>% dplyr::select(-Caloric_Value)

X_scaled <- scale(X)



# Bayesian Lasso 
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

# Fat multicollinearity
plot(unlist(FOOD_DATA[, "Saturated_Fats"] + FOOD_DATA[, "Monounsaturated_Fats"] + FOOD_DATA[, "Polyunsaturated_Fats"]),
    unlist(FOOD_DATA[, "Fat"]))
abline(0, 1)

# Lasso without fat
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

# Comparison of lasso with and without fat
lasso_with_fat <- lm(`Caloric_Value` ~ `Fat` + `Saturated_Fats` + `Monounsaturated_Fats` + `Polyunsaturated_Fats` + 
          `Carbohydrates` + `Protein` + `Water` + `Vitamin_A` + `Magnesium` + `Phosphorus`, data = FOOD_DATA)

lasso_without_fat <- lm(`Caloric_Value` ~ `Vitamin_E` + `Saturated_Fats` + `Monounsaturated_Fats` + `Polyunsaturated_Fats` + 
          `Carbohydrates` + `Protein` + `Cholesterol` + `Vitamin_A` + `Magnesium` + `Phosphorus`, data = FOOD_DATA)
vif(lasso_with_fat)
vif(lasso_without_fat)
summary(lasso_with_fat)
summary(lasso_without_fat)

# Investigating foods where sum of fats greater than fat
fat_check_df <- FOOD_DATA %>% filter(`Saturated_Fats` + `Monounsaturated_Fats` + `Polyunsaturated_Fats` > `Fat`)


# Spike Slab Prior
set.seed(123)

df_model <- data.frame(
  y = as.numeric(y),
  as.data.frame(X_scaled)
)

ss_fit <- lm.spike(
  y ~ .,
  data  = df_model,
  niter = 6000,                 # <-- only niter here
  expected.model.size = 10    # your prior belief of how many predictors matter
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

# Doing Bayesian linear regression with all predictors
# Using a conjugate Normal-Inv-Gamma prior posterior

# prior hyperparameters (noninformative)

n <- nrow(X_scaled)
p <- ncol(X_scaled)

xi <- rep(0, p)
Omega <- 10^3*diag(p)
alpha <- 0.01
b <- 0.01 # using b instead of beta for rate parameter for inv gamma to prevent confusion with betas coeffient in regression

Omega_inverse <- solve(Omega)
Q_beta <- t(as.matrix(X_scaled))%*%X_scaled + Omega_inverse
Q_beta_inverse <- solve(Q_beta)
l_beta <- t(as.matrix(X_scaled))%*%y + Omega_inverse%*%xi

# sampling posterior
samples <- 10^4
sigma_2_samples <- rinvgamma(samples, alpha + n/2, b + t(y)%*%y/2 + 
                              t(xi)%*%Omega_inverse%*%xi/2 - 
                              t(l_beta)%*%Q_beta_inverse%*%l_beta/2)

beta_samples <- matrix(NA, nrow = samples, ncol = p)
for (s in 1:samples) {
  beta_samples[s, ] <- mvrnorm(1, Q_beta_inverse%*%l_beta, 
                               Sigma = sigma_2_samples[s]*Q_beta_inverse)
}

hist(sigma_2_samples, xlab = expression(sigma^2))

# looking at credible interval of betas
beta_summary <- data.frame(
  term = c("intercept", colnames(X[, -1])),
  mean = colMeans(beta_samples),
  low = apply(beta_samples, 2, quantile, 0.025),
  high = apply(beta_samples, 2, quantile, 0.975)
)
beta_summary

selected_betas <- beta_summary[beta_summary[, "low"] > 0 | beta_summary[, "high"] < 0, ]
selected_betas


# EDA Plots

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
