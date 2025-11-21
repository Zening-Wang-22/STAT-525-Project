library(dplyr)
library(readr)
library(tidyverse)
library(monomvn)

FOOD_DATA_GROUP1 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP1.csv")
FOOD_DATA_GROUP2 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP2.csv")
FOOD_DATA_GROUP3 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP3.csv")
FOOD_DATA_GROUP4 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP4.csv")
FOOD_DATA_GROUP5 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP5.csv")

FOOD_DATA <- bind_rows(
  mutate(FOOD_DATA_GROUP1, group = "group1"),
  mutate(FOOD_DATA_GROUP2, group = "group2"),
  mutate(FOOD_DATA_GROUP3, group = "group3"),
  mutate(FOOD_DATA_GROUP4, group = "group4"),
  mutate(FOOD_DATA_GROUP5, group = "group5")
) %>%
  dplyr::select(-1, -2) 

# Preprocess
df <- FOOD_DATA %>%
  dplyr::select(-food, -group) %>%   # not predictors for now
  drop_na()

y <- df$`Caloric Value`
X <- df %>% dplyr::select(-`Caloric Value`)

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





