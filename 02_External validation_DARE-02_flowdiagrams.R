#:--------------------------------------------------------------------------------
#DARE CLEANING 

#Applying study exclusion criteria for external validation of C-peptide models
#in DARE

#ESM Figure 2
#:---------------------------------------------------------------------------------
#load libraries ------------------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(naniar)
library(Hmisc)

#load data ----------------------------------------------------------------------
load("~/Other/2024/YeunYi Dundee iDiabetes DARE/DARE_cpep.RData")

#flow diagram ---------------------------------------------------------------------
#keep only those with a cpep measurement
DARE_cpep <- DARE_cpep %>%
  filter(!is.na(GoldStndCpepPmolL))
#keep only those on insulin
DARE_cpep <- DARE_cpep %>%
  filter(insulinTreated == 1)
#keep only those with clinical diagnosis of type 2
DARE_cpep <- DARE_cpep %>%
  filter(TypeofDiabetes == "Type 2")
#keep only those without a positive GAD
DARE_cpep <- DARE_cpep %>%
  filter(GADPositive975 == 0 | is.na(GADPositive975))
#REMOVE THOSE WITH AGE AT TIME OF CPEPTIDE TESTING < 35
DARE_cpep <- DARE_cpep %>%
  filter(age_cpep >= 35 | (is.na(age_cpep) & Age >= 35))
#remove those with BMI outside 3 years
DARE_cpep <- DARE_cpep %>%
  filter(dur_bmi_cpep < 5 | (is.na(dur_bmi_cpep) & dur_bmi_rec < 5))
#remove those with creatinine >=300
DARE_cpep <- DARE_cpep %>%
  filter(creatinine < 300)

#remove those started insulin < 12 months
table(DARE_cpep$TimeToInsulin < 1, DARE_cpep$MonthsbetweenDmdiagnosisand, useNA = "ifany")
DARE_cpep <- DARE_cpep %>%
  filter(TimeToInsulin >= 1 | (MonthsbetweenDmdiagnosisand == ">12 months" & is.na(TimeToInsulin)))

# Remove participants with missing model variables
DARE_cpep_cc <- DARE_cpep%>%
  drop_na(sex, Age, BMI, DiabDuration, OnSU, NewHDL, NewALT, BasalBolus, NewHBA1C, 
          creatinine, NewTG, Cpep600, Cpep900)
# -------------------------------------------------------------------------------------------------
#Model all complete case
#Look at overlapping missingness in dataset
missing_vars <- DARE_cpep%>%
  select(sex, Age, BMI, DiabDuration, OnSU, NewHDL, NewALT, BasalBolus, NewHBA1C, 
          creatinine, NewTG, Cpep600, Cpep900)
vis_miss(missing_vars)
# n=3 missing because of missing basal bolus
# n=1 missing becuase of missing BMI

#Time to insulin imputation  ---------------------------------------------------------------
median(DARE_cpep_cc$TimeToInsulin, na.rm = TRUE)
#median = 7
summary(DARE_cpep_cc$TimeToInsulin)
summary(DARE_cpep$TimeToInsulin)
#with mising GAD
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.000   3.000   7.000   7.928  11.000  45.000      28 

#Perform histogram on those (n=31) that have YearContinousInsulin + time to insulin "<12months" (n=12)
ggplot(DARE_cpep_cc, aes(x=TimeToInsulin)) + 
  geom_histogram() +
  theme_classic() +
  xlab("Time to insulin (years) for those with YearContinousInsulin and \"<12months\"")

#impute or replace with median -------------------------------------------------------------
DARE_cpep_cc <- DARE_cpep_cc %>%
  mutate(TimeToInsulin = ifelse(is.na(TimeToInsulin), 7, TimeToInsulin))

#Save DARE dataset for type 2s with c-peptide and complete case model variables
#this has BMI within 5 years
save(DARE_cpep_cc, file = "DARE_cpep_cc.RData")

