library(tidyverse)
library(janitor)
library(recapr)

# loading some helper functions
source("FDS_2025/R/recapr_prep.R")

# reading raw data
event1_raw <- read_csv("FDS_2025/flat_data/event1.csv", skip=1) %>%
  remove_empty() %>%
  filter(!is.na(Unique)) %>%
  mutate(`Tag Number` = as.character(`Tag Number`)) %>%
  rename(FL = `Fork Length (mm)`)
event2_raw <- read_csv("FDS_2025/flat_data/event2.csv", skip=1) %>%
  remove_empty() %>%
  filter(!is.na(Unique)) %>%
  mutate(`Tag Number` = as.character(`Tag Number`)) %>%
  rename(FL = `Fork Length (mm)`) %>%
  rename(FLc = `Fork Length Corrected (mm)`)

dim(event1_raw)
dim(event2_raw)

summary(event1_raw)
summary(event2_raw)


# creating truncated data tables
# c for corrected (length correction)
event1_550 <- event1_raw %>%
  filter(FL >= 550)
event2_550 <- event2_raw %>%
  filter(FL >= 550)
event2_550c <- event2_raw %>%
  filter(FLc >= 550)
event1_600 <- event1_raw %>%
  filter(FL >= 600)
event2_600 <- event2_raw %>%
  filter(FL >= 600)
event2_600c <- event2_raw %>%
  filter(FLc >= 600)

# creating MR objects from recapr_prep
mr550 <- recapr_prep(ID="Tag Number", event1=event1_550, event2=event2_550)
mr600 <- recapr_prep(ID="Tag Number", event1=event1_600, event2=event2_600)
mr550c <- recapr_prep(ID="Tag Number", event1=event1_550, event2=event2_550c)
mr600c <- recapr_prep(ID="Tag Number", event1=event1_600, event2=event2_600c)

# checking sample sizes against spreadsheet
sapply(mr550, \(x) sapply(x, nrow))
sapply(mr550c, \(x) sapply(x, nrow))
sapply(mr600, \(x) sapply(x, nrow))
sapply(mr600c, \(x) sapply(x, nrow))


# KS tests from non-stratified samples (all possible ways)
with(mr550, ks.test(input_data$event1$FL, recaps$matched$FL_event1))  # FTR
with(mr550, ks.test(input_data$event2$FL, recaps$matched$FL_event2))  # reject
with(mr550, ks.test(input_data$event2$FL, recaps$matched$FL_event1))  # reject

with(mr550c, ks.test(input_data$event1$FL, recaps$matched$FL_event1))  # FTR
with(mr550c, ks.test(input_data$event2$FLc, recaps$matched$FLc_event2))  # reject
# with(mr550c, ks.test(input_data$event2$FLc, recaps$matched$FL_event1))  # reject

with(mr600, ks.test(input_data$event1$FL, recaps$matched$FL_event1))  # FTR
with(mr600, ks.test(input_data$event2$FL, recaps$matched$FL_event2))  # FTR
with(mr600, ks.test(input_data$event2$FL, recaps$matched$FL_event1))  # FTR

with(mr600c, ks.test(input_data$event1$FL, recaps$matched$FL_event1))  # FTR
with(mr600c, ks.test(input_data$event2$FLc, recaps$matched$FLc_event2))  # FTR
# with(mr600c, ks.test(input_data$event2$FLc, recaps$matched$FL_event1))  # FTR


# creating stratified MR objects (all possible ways)
mr550_strat <- stratify(mr550,
                        event_names = c("event1", "event2"),
                        column_names = c("FL","FL"),
                        breaks = c(550, 650, 1200))
mr600_strat <- stratify(mr600,
                        event_names = c("event1", "event2"),
                        column_names = c("FL","FL"),
                        breaks = c(600, 650, 1200))
mr550c_strat <- stratify(mr550c,
                        event_names = c("event1", "event2"),
                        column_names = c("FL","FLc"),
                        breaks = c(550, 650, 1200))
mr600c_strat <- stratify(mr600c,
                        event_names = c("event1", "event2"),
                        column_names = c("FL","FLc"),
                        breaks = c(600, 650, 1200))

# KS tests for each stratum - all possible ways.  This is tedious!!
with(mr550_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr550_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))
with(mr550_strat,  # p-value = 0.05149
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[1]],
             recaps$matched$FL_event2[recaps$matched$FL_strat_event2==levels(recaps$matched$FL_strat_event2)[1]]))
with(mr550_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[2]],
             recaps$matched$FL_event2[recaps$matched$FL_strat_event2==levels(recaps$matched$FL_strat_event2)[2]]))
with(mr550_strat,  # p-value = 0.02559
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr550_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))

with(mr600_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr600_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))
with(mr600_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[1]],
             recaps$matched$FL_event2[recaps$matched$FL_strat_event2==levels(recaps$matched$FL_strat_event2)[1]]))
with(mr600_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[2]],
             recaps$matched$FL_event2[recaps$matched$FL_strat_event2==levels(recaps$matched$FL_strat_event2)[2]]))
with(mr600_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr600_strat,  # FTR
     ks.test(input_data$event2$FL[input_data$event2$FL_strat==levels(input_data$event2$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))

with(mr550c_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr550c_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))
with(mr550c_strat,  # p-value = 0.0628
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[1]],
             recaps$matched$FLc_event2[recaps$matched$FLc_strat_event2==levels(recaps$matched$FLc_strat_event2)[1]]))
with(mr550c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[2]],
             recaps$matched$FLc_event2[recaps$matched$FLc_strat_event2==levels(recaps$matched$FLc_strat_event2)[2]]))
with(mr550c_strat,  # p-value = 0.0628
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr550c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))

with(mr600c_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr600c_strat,  # FTR
     ks.test(input_data$event1$FL[input_data$event1$FL_strat==levels(input_data$event1$FL_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))
with(mr600c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[1]],
             recaps$matched$FLc_event2[recaps$matched$FLc_strat_event2==levels(recaps$matched$FLc_strat_event2)[1]]))
with(mr600c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[2]],
             recaps$matched$FLc_event2[recaps$matched$FLc_strat_event2==levels(recaps$matched$FLc_strat_event2)[2]]))
with(mr600c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[1]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[1]]))
with(mr600c_strat,  # FTR
     ks.test(input_data$event2$FLc[input_data$event2$FLc_strat==levels(input_data$event2$FLc_strat)[2]],
             recaps$matched$FL_event1[recaps$matched$FL_strat_event1==levels(recaps$matched$FL_strat_event1)[2]]))


# sample sizes, all possible ways.
with(mr550_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FL_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FL_strat_event2)
  print(n1)
  print(n2)
  print(m2_1)
  print(m2_2)
  print(NChapman(n1, n2, m2_1))
  print(seChapman(n1, n2, m2_1))
  print(Nstrat(n1, n2, m2_1))
  print(sestrat(n1, n2, m2_1))
})

with(mr550c_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FLc_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FLc_strat_event2)
  print(n1)
  print(n2)
  print(m2_1)
  print(m2_2)
  print(NChapman(n1, n2, m2_1))
  print(seChapman(n1, n2, m2_1))
  print(Nstrat(n1, n2, m2_1))
  print(sestrat(n1, n2, m2_1))
})

with(mr600_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FL_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FL_strat_event2)
  print(n1)
  print(n2)
  print(m2_1)
  print(m2_2)
  print(NChapman(n1, n2, m2_1))
  print(seChapman(n1, n2, m2_1))
  print(Nstrat(n1, n2, m2_1))
  print(sestrat(n1, n2, m2_1))
})

with(mr600c_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FLc_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FLc_strat_event2)
  print(n1)
  print(n2)
  print(m2_1)
  print(m2_2)
  print(NChapman(n1, n2, m2_1))
  print(seChapman(n1, n2, m2_1))
  print(Nstrat(n1, n2, m2_1))
  print(sestrat(n1, n2, m2_1))
})
