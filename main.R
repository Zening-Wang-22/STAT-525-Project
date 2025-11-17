library(dplyr)
library(readr)

FOOD_DATA_GROUP1 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP1.csv")
FOOD_DATA_GROUP2 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP2.csv")
FOOD_DATA_GROUP3 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP3.csv")
FOOD_DATA_GROUP4 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP4.csv")
FOOD_DATA_GROUP5 <- read_csv("FINAL FOOD DATASET/FOOD-DATA-GROUP5.csv")

FOOD_DATA_ALL <- bind_rows(
  mutate(FOOD_DATA_GROUP1, group = "group1"),
  mutate(FOOD_DATA_GROUP2, group = "group2"),
  mutate(FOOD_DATA_GROUP3, group = "group3"),
  mutate(FOOD_DATA_GROUP4, group = "group4"),
  mutate(FOOD_DATA_GROUP5, group = "group5")
) %>%
  select(-1, -2) 
