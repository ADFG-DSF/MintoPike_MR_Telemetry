library(tidyverse)
library(janitor)
library(recapr)
library(dsftools)

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
  rename(FLc = `Fork Length Corrected (mm)`) %>%
  mutate(gear = sapply(strsplit(Comments, split="Gear= "), "[", 2)) %>%
  mutate(gear = factor(ifelse(gear=="G", "GN", gear)))

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


# printing sample sizes and abundance estimates, all possible ways.
with(mr550_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FL_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FL_strat_event2)
  print("n1")
  print(n1)
  print("n2")
  print(n2)
  print("m2 v1")
  print(m2_1)
  print("m2 v2")
  print(m2_2)
  print("Nhat per strat")
  print(NChapman(n1, n2, m2_1))
  print("SE per strat")
  print(seChapman(n1, n2, m2_1))
  print("Nhat tot - stratified")
  print(Nstrat(n1, n2, m2_1))
  print("SE tot - stratified")
  print(sestrat(n1, n2, m2_1))
})

with(mr550c_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FLc_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FLc_strat_event2)
  print("n1")
  print(n1)
  print("n2")
  print(n2)
  print("m2 v1")
  print(m2_1)
  print("m2 v2")
  print(m2_2)
  print("Nhat per strat")
  print(NChapman(n1, n2, m2_1))
  print("SE per strat")
  print(seChapman(n1, n2, m2_1))
  print("Nhat tot - stratified")
  print(Nstrat(n1, n2, m2_1))
  print("SE tot - stratified")
  print(sestrat(n1, n2, m2_1))
})

with(mr600_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FL_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FL_strat_event2)
  print("n1")
  print(n1)
  print("n2")
  print(n2)
  print("m2 v1")
  print(m2_1)
  print("m2 v2")
  print(m2_2)
  print("Nhat per strat")
  print(NChapman(n1, n2, m2_1))
  print("SE per strat")
  print(seChapman(n1, n2, m2_1))
  print("Nhat tot - stratified")
  print(Nstrat(n1, n2, m2_1))
  print("SE tot - stratified")
  print(sestrat(n1, n2, m2_1))
})

with(mr600c_strat, {
  n1 <- table(input_data$event1$FL_strat)
  n2 <- table(input_data$event2$FLc_strat)
  m2_1 <- table(recaps$matched$FL_strat_event1)
  m2_2 <- table(recaps$matched$FLc_strat_event2)
  print("n1")
  print(n1)
  print("n2")
  print(n2)
  print("m2 v1")
  print(m2_1)
  print("m2 v2")
  print(m2_2)
  print("Nhat per strat")
  print(NChapman(n1, n2, m2_1))
  print("SE per strat")
  print(seChapman(n1, n2, m2_1))
  print("Nhat tot - stratified")
  print(Nstrat(n1, n2, m2_1))
  print("SE tot - stratified")
  print(sestrat(n1, n2, m2_1))
})

# to do:
# - length comps (from second event)
# - alt abundance of 600+ using length comp
# - length comp by gear
# - FIGURE OUT SPATIAL BIAS ISSUE  <--- honestly i don't know if i can do this
# - WRITE THIS ALL UP

n1 <- table(mr550_strat$input_data$event1$FL_strat)
n2 <- table(mr550_strat$input_data$event2$FL_strat)
m2_1 <- table(mr550_strat$recaps$matched$FL_strat_event1)
m2_2 <- table(mr550_strat$recaps$matched$FL_strat_event2)

Nhat_vec <- NChapman(n1, n2, m2_1)
seNhat_vec <-   seChapman(n1, n2, m2_1)
Nhat_strat <- Nstrat(n1, n2, m2_1)
seNhat_strat <- sestrat(n1, n2, m2_1)

# Length comp using 550 truncation and 650 stratification
lengthbreaks <- c(550, 600, 650, 700, 800, 900, 1000, 1200) #550 to 1146
ASL_table(age = cut(mr550_strat$input_data$event2$FL, breaks=lengthbreaks, right=FALSE, dig.lab = 4),
          stratum=as.numeric(mr550_strat$input_data$event2$FL_strat),
          Nhat=as.numeric(Nhat_vec),
          se_Nhat = as.numeric(seNhat_vec))

# this gets pop est of 600+ using the 650 strat
ASL_table(age = cut(mr550_strat$input_data$event2$FL, breaks=c(550, 600, 1200), right=FALSE, dig.lab = 4),
          stratum=as.numeric(mr550_strat$input_data$event2$FL_strat),
          Nhat=as.numeric(Nhat_vec),
          se_Nhat = as.numeric(seNhat_vec))

# this gets prop est & mn lengths for each gear type, using 650 strat?? not sure if this is relevant
ASL_table(age = mr550_strat$input_data$event2$gear,
          length = mr550_strat$input_data$event2$FL,
          stratum=as.numeric(mr550_strat$input_data$event2$FL_strat),
          Nhat=as.numeric(Nhat_vec),
          se_Nhat = as.numeric(seNhat_vec))

# plotting FL by gear
boxplot(FL~gear, data=event2_550)

# non-stratified summary  - think about which is more appropriate
event2_550 %>%
  group_by(gear) %>%
  summarise(n=length(FL),
            mean_FL=mean(FL),
            min_FL=min(FL),
            max_FL=max(FL),
            sd_FL=sd(FL)) %>%
  mutate(se_FL = sd_FL/sqrt(n)) %>%
  mutate(p_hat = n/sum(n)) %>%
  mutate(se_p_hat = sqrt(p_hat*(1-p_hat)/(sum(n)-1)))

# same thing with ASL_table
ASL_table(age = mr550_strat$input_data$event2$gear,
          length = mr550_strat$input_data$event2$FL)
