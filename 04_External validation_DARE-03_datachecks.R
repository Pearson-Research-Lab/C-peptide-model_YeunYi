#:-------------------------------------------------------------------------------
#DARE CLEANING 

#Examining DARE data for external validation of C-peptide models
#to maximise data inclusion for validation 

#:-------------------------------------------------------------------------------

#load libraries -----------------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(naniar)
library(Hmisc)

#load data ----------------------------------------------------------------------
load("DARE_cpep_cc.RData")
insulin_types <- read_xlsx("sorted insulin.xlsx")

#Cpep_year ----------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(is.na(cpep_year)) %>%
  count()
#93/250 participants have DiabDuration coded via difference between recruitment 
#and diabetes diagnosis, not Cpeptide

#Basal bolus --------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(BasalBolus == 0) %>%
  mutate(
    insulin_combo = case_when(
      (Insulin1Type %in% insulin_types$Long 
       & Insulin2Type %in% insulin_types$Rapid) | 
        (Insulin1Type %in% insulin_types$Rapid 
         & Insulin2Type %in% insulin_types$Long) | 
        (Insulin1Type %in% insulin_types$Rapid 
         & Insulin3Type %in% insulin_types$Long)  ~ "Rapid and Long",
      (Insulin1Type %in% insulin_types$Long 
       & Insulin2Type %in% insulin_types$Mixed) | 
        (Insulin1Type %in% insulin_types$Mixed 
         & Insulin2Type %in% insulin_types$Long) ~ "Mixed and Long",
      (Insulin1Type %in% insulin_types$Mixed 
       & Insulin2Type %in% insulin_types$Rapid 
       & Insulin3Type %in% insulin_types$Long) ~ "Rapid, Mixed and Long"
    )
  ) %>%
  group_by(insulin_combo) %>%
  count()

#Rapid and Long           35
#Rapid, Mixed and Long     1

#TRIGS -------------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(TGDate_min > 365) %>%
  count()
DARE_cpep_cc %>%
  filter(is.na(trigvalue)) %>%
  count()
DARE_cpep_cc %>%
  filter(TGDate_min > 365 | is.na(trigvalue)) %>%
  count()
#170/250 has a replaced triglyceride
#ALT --------------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(ALTDate_min > 365) %>%
  count()
DARE_cpep_cc %>%
  filter(is.na(ALTvalue)) %>%
  count()
DARE_cpep_cc %>%
  filter(ALTDate_min > 365 | is.na(ALTvalue)) %>%
  count()
#45/250 has a replaced ALT
#HDL --------------------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(HDLDate_min > 365) %>%
  count()
DARE_cpep_cc %>%
  filter(is.na(HDLvalue)) %>%
  count()
DARE_cpep_cc %>%
  filter(HDLDate_min > 365 | is.na(HDLvalue)) %>%
  count()
#106/250 has a replaced HDL
#HBA1C -------------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(HBA1CDate_min > 182.5) %>%
  count()
DARE_cpep_cc %>%
  filter(is.na(HBA1Cvalue)) %>%
  count()
DARE_cpep_cc %>%
  filter(HBA1CDate_min > 182.5 | is.na(HBA1Cvalue)) %>%
  count()
#26/250 has a replaced HBA1C

#Age ------------------------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(is.na(age_cpep)) %>%
  count()
#93/250 has age = age at study recruitment

#cREATININE DATE -----------------------------------------------------------------------
DARE_cpep_cc %>%
  filter(dur_creatinine > 1) %>%
  count()
#66/250 outside of 1 year
DARE_cpep_cc %>%
  filter(dur_creatinine > 1.5) %>%
  count()
#59/250 outside of 1.5 year
DARE_cpep_cc %>%
  filter(dur_creatinine > 2) %>%
  count()
#52/250 outside of 2 year