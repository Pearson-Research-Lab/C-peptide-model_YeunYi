#:----------------------------------------------------------------------------------
#DARE CLEANING

#DARE used as external validation dataset for C-peptide models 
#:--------------------------------------------------------------------------------
#load libraries ------------------------------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(naniar)
library(Hmisc)

#Datset preparation and loading ------------------------------------------------------
##load data -------------------------------------------------------------------------------------
DARE <- read_dta("anita uptodate goldstandard.dta")
DARE_lipids <- read_dta("DareLipidDataSurvivalModel.dta")
insulin_types <- read_xlsx("sorted insulin.xlsx")

##join data ------------------------------------------------------------------------
DARE <- left_join(DARE, DARE_lipids, by = "ID")

###########################################################################################
#DEFINING VARIABLES 
##############################################################################################
##DiabDuration ---------------------------------------------------------------------------------
#calculate diabetes duration from C-peptide date  to YearofDIganosis
DARE <- DARE %>%
  mutate(
    #extract year of cpeptide measurement from date
    cpep_year  = format(as.POSIXct(GoldStndCpepDte), "%Y"), 
    #Calculate DiabDuration
    #if they have a year at C-peptide measurement, minus year of diagnosis, 
    #other wise use Duration of diabetes (years) at recruitment
    DiabDuration = ifelse(
      !is.na(cpep_year), 
      as.numeric(cpep_year) - as.numeric(YearofDiagnosis), 
      DurationYrsdiabetes)
  )
summary(DARE$DiabDuration)

#TIMETOINSULIN  --------------------------------------------------------------------------
#(DURATION BETWEEN DATE OF DIAGNOSIS AND FIRST INSULIN PRESCRIPTION)
#calculate Time to insulin
DARE <- DARE %>%
  mutate(
    TimeToInsulin = Timetoinsulinmonths/12
  )
summary(DARE$TimeToInsulin)

#BasalBolus -----------------------------------------------------------------------
#Investigate combinations of relevant variables
#check number of people who have two types of insulin in the below two variables
DARE %>%
  filter(Insulin1Type != "" & Insulin2Type != "") %>%
  count()
#check number of people who have two types of insulin in the below three 
#variables
DARE %>%
  filter(Insulin1Type != "" & Insulin2Type != ""& Insulin3Type != "") %>%
  count()
#check number of people who have 2 types of insulin, in Insulin1Type and 
#Insulin3Type, and empty Insulin2Type
DARE %>%
  filter(Insulin1Type != "" & Insulin2Type == ""& Insulin3Type != "") %>%
  select(Insulin1Type, Insulin3Type)
#check if there are any missing insulin1type but have in the other two
DARE %>%
  filter(Insulin1Type == "" & Insulin2Type != ""& Insulin3Type != "") %>%
  count()
#none
#count the number who have a missing insulin1type
DARE %>%
  filter(Insulin1Type == "") %>%
  count()

#Calculate BasalBolus
DARE <- DARE %>%
  mutate(
    BasalBolus = ifelse(
      (Insulin1Type %in% insulin_types$Long 
       & Insulin2Type %in% insulin_types$Rapid) | 
      (Insulin1Type %in% insulin_types$Rapid 
       & Insulin2Type %in% insulin_types$Long) | 
      (Insulin1Type %in% insulin_types$Rapid 
       & Insulin3Type %in% insulin_types$Long) | 
      (Insulin1Type %in% insulin_types$Long 
       & Insulin2Type %in% insulin_types$Mixed) | 
      (Insulin1Type %in% insulin_types$Mixed 
       & Insulin2Type %in% insulin_types$Long) | 
      (Insulin1Type %in% insulin_types$Mixed 
       & Insulin2Type %in% insulin_types$Rapid 
       & Insulin3Type %in% insulin_types$Long), 
      "Yes", 
      ifelse(Insulin1Type == "", NA, "No")))
table(DARE$BasalBolus, useNA = "ifany")

## Biomarkers -------------------------------------------------------------------
### Calulating TG (within 12 months of cpeptide measurement) --------------------
#build dataset to pivot longer
#Select all variables that start with TRIG and corresponding SpecimenDate
#There should be 75 of each
Trig_clean <- DARE %>%
  select(ID, DateSeen, GoldStndCpepDte, 
         starts_with("TRIG"), 
         starts_with("SpecimenDate"))
#Add suffixes to TRIGs
Trig_clean <- Trig_clean %>%
  rename_with(~paste0(.,"_trigvalue"), TRIG1:TRIG75)
#Add suffixes to DATES
Trig_clean <- Trig_clean %>%
  rename_with(~paste0(.,"_date"), SpecimenDate1:SpecimenDate75)
#Replace all original variables starting with "SpecimenDate" with "TRIG"
names(Trig_clean) <- Trig_clean %>%
  names() %>%
  str_replace("SpecimenDate", "TRIG")

#pivot longer
### Trig_clean is dataframe with all C-peptide levels, C-peptide dates, 
#Triglycerides levels, Triglycerides dates for each individual
Trig_clean <- pivot_longer(
  Trig_clean, 
  cols = starts_with("TRIG"), 
  names_to = c("measurment", ".value"), 
  names_sep = "_")

###Cpep_Date is C-peptide date
Trig_clean$GoldStndCpepDte <- as.Date(Trig_clean$GoldStndCpepDte)  
### date is Triglycerides date
Trig_clean$date <- as.Date(Trig_clean$date)  
Trig_clean$DateSeen <- as.Date(Trig_clean$DateSeen)

#Identify number of individuals without any trig values (across 75 measurements)
trig_missingness <- Trig_clean %>%
  group_by(ID) %>%
  summarise(no_missing = sum(is.na(trigvalue)), 
            have_value = sum(!is.na(trigvalue))) %>%
  filter(have_value == 0) %>%
  count()
#n=2708/5993 have no trig values in any visit

### TG3 is filtered dataframe with Triglycerides within 365 days of 
#C-peptide measurement for each individual
#this is filtering for those measurements that have a trig value, then slicing 
#the minimum date of those per individual
#this is to maximise obtaining trig values over minimum date
#if they don't have a trig value then just slicing on minimum date
TG3<- Trig_clean %>% 
  group_by(ID) %>% 
  mutate(
    TGDate_min= ifelse(
      !is.na(GoldStndCpepDte), 
      abs(GoldStndCpepDte - date), 
      abs(DateSeen - date)), 
    no_missing = sum(is.na(trigvalue)), 
    have_value = sum(!is.na(trigvalue)), 
    Datemin_have = ifelse(sum(!is.na(trigvalue)) > 0, "Yes", "No"),
    Datemin_have_value = ifelse(!is.na(trigvalue), TGDate_min, NA),
    rank = order(TGDate_min)) %>% 
  filter(Datemin_have == "Yes" 
         & Datemin_have_value == min(Datemin_have_value, na.rm = TRUE) 
         | Datemin_have == "No" & rank == 1) %>% 
  slice_min(TGDate_min, with_ties = F)

### check for missing value
TG3 %>% 
  ungroup() %>% 
  filter(is.na(trigvalue)) %>% 
  count() 
#with just slicing on minimum date end up with n=2708/5993 without a trig value; 
#and therefore n=3285/5993 with trig values

TG3 %>% 
  ungroup() %>% 
  filter(is.na(trigvalue) | TGDate_min > 365) %>% 
  count()
#with just slicing on minimum date end up with n=4115 without a trig value or date > 365
#therefore n=4115/5993 would have a replaced trig value 
#and n=1878/5993 with trig values in the right time frame

#identify median value that should replace missing values with
TG3 %>%
  filter(!is.na(trigvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(trigvalue, na.rm = TRUE),
    n = n()
  )
#median is 1.5 from n=3285
median(TG3$trigvalue, na.rm = TRUE)
#median is 1.5

#this identifies median value in those within 1 year of cpep
#use this value
TG3 %>%
  filter(TGDate_min <= 365 & !is.na(trigvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(trigvalue, na.rm = TRUE),
    n = n()
  )
#median is 1.46 from n=1878

### I replaced missing values and those >=365 days with median TG 1.46 
#(median of all TG values within 12 months of C-pep)
### 28 individuals
TG3 %>%
  ungroup() %>%
  filter(TG3$TGDate_min >365) %>%
  count()
#Replace 1992/5993 trig values with dates outside 365 range
TG3$NewTG <-ifelse((TG3$TGDate_min >365), 1.46, TG3$trigvalue)
TG3 %>%
  ungroup() %>%
  filter(is.na(TG3$NewTG)) %>%
  count()
#Replace further 2123/5993 missing trig values 
TG3$NewTG <- ifelse(is.na(TG3$NewTG), 1.46, TG3$NewTG)

### The filtered dataframe (TG3) is then merged with the main dataset 
#(to add a Triglycerides column) - later

### ALT (WITHIN 12 MONTHS OF C-PEPTIDE MEASUREMENT) ------------------------------
#build dataset to pivot longer
#Repeat process of TRIGs for ALT
ALT_clean <- DARE %>%
  select(ID, DateSeen,GoldStndCpepDte, 
         starts_with("ALT"), 
         starts_with("SpecimenDate"))
ALT_clean <- ALT_clean %>%
  rename_with(~paste0(.,"_ALTvalue"), ALT1:ALT75)
ALT_clean <- ALT_clean %>%
  rename_with(~paste0(.,"_date"), SpecimenDate1:SpecimenDate75)
names(ALT_clean) <- ALT_clean %>%
  names() %>%
  str_replace("SpecimenDate", "ALT")

#pivotlonger
### ALT_clean is dataframe with all C-peptide levels, C-peptide dates, 
#Triglycerides levels, Triglycerides dates for each individual
ALT_clean <- pivot_longer(
  ALT_clean, 
  cols = starts_with("ALT"), 
  names_to = c("measurment", ".value"), 
  names_sep = "_")
###GoldStndCpepDte is C-peptide date
ALT_clean$GoldStndCpepDte <- as.Date(ALT_clean$GoldStndCpepDte)  
###date is ALT date
ALT_clean$date <- as.Date(ALT_clean$date)  
ALT_clean$DateSeen <- as.Date(ALT_clean$DateSeen)

#Identify number of individuals without any ALT values (across 75 measurements)
ALT_missingness <- ALT_clean %>%
  group_by(ID) %>%
  summarise(no_missing = sum(is.na(ALTvalue)), 
            have_value = sum(!is.na(ALTvalue))) %>%
  filter(have_value == 0) %>%
  count()
#n=1103/5993 have no ALT values in any visit

### ALT3 is filtered dataframe with ALT within 365 days of C-peptide measurement 
#for each individual
#this is filtering for those measurements that have a ALT value, then slicing 
#the minimum date of those per individual
#this is to maximise obtaining trig values over minimum date
#if they don't have a trig value then just slicing on minimum date
ALT3<- ALT_clean %>% 
  group_by(ID) %>% 
  mutate(
    ALTDate_min= ifelse(
      !is.na(GoldStndCpepDte), 
      abs(GoldStndCpepDte - date), 
      abs(DateSeen - date)), 
    no_missing = sum(is.na(ALTvalue)), 
    have_value = sum(!is.na(ALTvalue)), 
    Datemin_have = ifelse(sum(!is.na(ALTvalue)) > 0, "Yes", "No"),
    Datemin_have_value = ifelse(!is.na(ALTvalue), ALTDate_min, NA),
    rank = order(ALTDate_min)) %>% 
  filter(Datemin_have == "Yes" 
         & Datemin_have_value == min(Datemin_have_value, na.rm = TRUE) 
         | Datemin_have == "No" & rank == 1) %>% 
  slice_min(ALTDate_min, with_ties = F)

### check for missing value
ALT3 %>% 
  ungroup() %>% 
  filter(is.na(ALTvalue)) %>% 
  count() 
#with just slicing on minimum date end up with n=1103/5993 without a ALT value; 
#and therefore n=4890/5993 with ALT values
ALT3 %>% 
  ungroup() %>% 
  filter(is.na(ALTvalue) | ALTDate_min > 365) %>% 
  count()
#with just slicing on minimum date end up with n=4115 without a ALT value or date > 365
#therefore n=1854/5993 would have a replaced ALT value 
#and n=4139/5993 with ALT values in the right time frame

#identify median value that should replace missing values with
ALT3 %>%
  filter(!is.na(ALTvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(ALTvalue, na.rm = TRUE),
    n = n()
  )
median(ALT3$ALTvalue, na.rm = TRUE)
#median is 23
#this identifies median value in those within 1 year of cpep
#use this value
ALT3 %>%
  filter(ALTDate_min <= 365 & !is.na(ALTvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(ALTvalue, na.rm = TRUE),
    n = n()
  )
#Median is 23 from 4139 values

### I replaced missing values and those >=365 days with median alt 23 
#(median of all TG values within 12 months of C-pep)
### 28 individuals
ALT3 %>%
  filter(ALTDate_min >365) %>%
  ungroup() %>%
  count()
#replaced 894/5993 with dates outside 1 year
ALT3$NewALT <- ifelse((ALT3$ALTDate_min >365), 23, ALT3$ALTvalue)
ALT3 %>%
  filter(is.na(NewALT)) %>%
  ungroup() %>%
  count()
#replaced further 960/5993 with missing ALT
ALT3$NewALT <- ifelse(is.na(ALT3$NewALT), 23, ALT3$NewALT)

### HDL (WITHIN 12 MONTHS OF C-PEPTIDE MEASUREMTNE) -----------------------------
#build dataset to pivot longer
#Repeat process for TRIGs and ALT for HDL
HDL_clean <- DARE %>%
  select(ID, DateSeen, GoldStndCpepDte, 
         starts_with("HDL"), 
         starts_with("SpecimenDate"))
HDL_clean <- HDL_clean %>%
  rename_with(~paste0(.,"_HDLvalue"), HDL1:HDL75)
HDL_clean <- HDL_clean %>%
  rename_with(~paste0(.,"_date"), SpecimenDate1:SpecimenDate75)
names(HDL_clean) <- HDL_clean %>%
  names() %>%
  str_replace("SpecimenDate", "HDL")
#PIVOT LONGER
### HDL_clean is dataframe with all C-peptide levels, C-peptide dates, 
#HDL levels, HDL dates for each individual
HDL_clean <- pivot_longer(
  HDL_clean, 
  cols = starts_with("HDL"), 
  names_to = c("measurment", ".value"), 
  names_sep = "_")
###GoldStndCpepDte is C-peptide date
HDL_clean$GoldStndCpepDte <- as.Date(HDL_clean$GoldStndCpepDte)  
### date is HDL date
HDL_clean$date <- as.Date(HDL_clean$date)  
HDL_clean$DateSeen <- as.Date(HDL_clean$DateSeen)  

#Identify number of individuals without any HDL values (across 75 measurements)
HDL_missingness <- HDL_clean %>%
  group_by(ID) %>%
  summarise(no_missing = sum(is.na(HDLvalue)), 
            have_value = sum(!is.na(HDLvalue))) %>%
  filter(have_value == 0) %>%
  count()
#n=1502/5993 have no HDL values in any visit

#HDL3 is filtered dataframe with HDL within 365 days of C-peptide measurement 
#for each individual
#this is filtering for those measurements that have a HDL value, then slicing 
#the minimum date of those per individual
#this is to maximise obtaining HDL values over minimum date
#if they don't have a HDL value then just slicing on minimum date
HDL3<- HDL_clean %>% 
  group_by(ID) %>% 
  mutate(
    HDLDate_min= ifelse(
      !is.na(GoldStndCpepDte), 
      abs(GoldStndCpepDte - date), 
      abs(DateSeen - date)), 
    no_missing = sum(is.na(HDLvalue)), 
    have_value = sum(!is.na(HDLvalue)), 
    Datemin_have = ifelse(sum(!is.na(HDLvalue)) > 0, "Yes", "No"),
    Datemin_have_value = ifelse(!is.na(HDLvalue), HDLDate_min, NA),
    rank = order(HDLDate_min)) %>% 
  filter(Datemin_have == "Yes" 
         & Datemin_have_value == min(Datemin_have_value, na.rm = TRUE) 
         | Datemin_have == "No" & rank == 1) %>% 
  slice_min(HDLDate_min, with_ties = F)

### check for missing value
HDL3 %>% 
  ungroup() %>% 
  filter(is.na(HDLvalue)) %>% 
  count() 
#with just slicing on minimum date end up with n=1502/5993 without a HDL value; 
#and therefore n=4491/5993 with HDL values
HDL3 %>% 
  ungroup() %>% 
  filter(is.na(HDLvalue) | HDLDate_min > 365) %>% 
  count()
#with just slicing on minimum date end up with n=4115 without a HDL value 
#or date > 365
#therefore n=2863/5993 would have a replaced HDL value 
#and n=3130/5993 with HDL values in the right time frame

#identify median value that should replace missing values with
HDL3 %>%
  ungroup() %>%
  filter(!is.na(HDLvalue)) %>%
  summarise(
    median = median(HDLvalue, na.rm = TRUE),
    n = n()
  )
median(HDL3$HDLvalue, na.rm = TRUE)
#this identifies median value in those within 1 year of cpep
#use this value
HDL3 %>%
  ungroup() %>%
  filter(HDLDate_min <= 365 & !is.na(HDLvalue))   %>%
  summarise(
    median = median(HDLvalue, na.rm = TRUE),
    n = n()
  )
#Median is 1.29 OF N=3130

### I replaced missing values and those >=365 days with median HDL 1.29 
#(median of all TG values within 12 months of C-pep)
HDL3 %>%
  ungroup() %>%
  filter(HDLDate_min >365) %>%
  count()
#Replace 1652/5993 with values outside a year with 1.29
HDL3$NewHDL <- ifelse((HDL3$HDLDate_min >365), 1.29, HDL3$HDLvalue)
HDL3 %>%
  ungroup() %>%
  filter(is.na(NewHDL)) %>%
  count()
#replace fruther 1211/5993 with missing hdl values with 1.29
HDL3$NewHDL<- ifelse(is.na(HDL3$NewHDL), 1.29, HDL3$NewHDL)

### HBA1C (WITHIN 6 MONTHS OF CPEPTIDE MEAUSRE; mmol/mol) ------------------------
#build dataset to pivot longer
#Repeat same process as above with HbA1c
HBA1C_clean <- DARE %>%
  select(ID, DateSeen, GoldStndCpepDte, 
         starts_with("HbA1c"), starts_with("SpecimenDate")) %>%
  select(-starts_with("HbA1c_")) %>%
  select(-HbA1cAtRecrDate, 
         -HBA1cNearestToGAD, 
         -HBA1cNearestToGADDate, 
         -HBA1cNearestToGADPerc, 
         -HBA1CClosestToRec, 
         -HbA1cAtDiagnosis)
HBA1C_clean <- HBA1C_clean %>%
  rename_with(~paste0(.,"_HBA1Cvalue"), HbA1c1:HbA1c75)
HBA1C_clean <- HBA1C_clean %>%
  rename_with(~paste0(.,"_date"), SpecimenDate1:SpecimenDate75)
names(HBA1C_clean) <- HBA1C_clean %>%
  names() %>%
  str_replace("SpecimenDate", "HbA1c")

#pivot longer
### HBA1C_clean is dataframe with all C-peptide levels, C-peptide dates, 
#HBA1C levels, HBA1C dates for each individual
HBA1C_clean <- pivot_longer(
  HBA1C_clean, 
  cols = starts_with("HbA1c"), 
  names_to = c("measurment", ".value"), 
  names_sep = "_")
###GoldStndCpepDte is C-peptide date
HBA1C_clean$GoldStndCpepDte <- as.Date(HBA1C_clean$GoldStndCpepDte)  
### date is HBA1C date
HBA1C_clean$date <- as.Date(HBA1C_clean$date)  
HBA1C_clean$DateSeen <- as.Date(HBA1C_clean$DateSeen)  

#Identify number of individuals without any HBA1C values (across 75 measurements)
HBA1C_missingness <- HBA1C_clean %>%
  group_by(ID) %>%
  summarise(no_missing = sum(is.na(HBA1Cvalue)), 
            have_value = sum(!is.na(HBA1Cvalue))) %>%
  filter(have_value == 0) %>%
  count()
#n=614/5993 have no HbA1c values in any visit

### HBA1C3 is filtered dataframe with Triglycerides within 365 days of 
#C-peptide measurement for each individual
#this is filtering for those measurements that have a HBA1C value, 
#then slicing the minimum date of those per individual
#this is to maximise obtaining HBA1C values over minimum date
#if they don't have a HBA1C value then just slicing on minimum date
HBA1C3<- HBA1C_clean %>% 
  group_by(ID) %>% 
  mutate(
    HBA1CDate_min= ifelse(
      !is.na(GoldStndCpepDte), 
      abs(GoldStndCpepDte - date), 
      abs(DateSeen - date)), 
    no_missing = sum(is.na(HBA1Cvalue)), 
    have_value = sum(!is.na(HBA1Cvalue)), 
    Datemin_have = ifelse(sum(!is.na(HBA1Cvalue)) > 0, "Yes", "No"),
    Datemin_have_value = ifelse(!is.na(HBA1Cvalue), HBA1CDate_min, NA),
    rank = order(HBA1CDate_min)) %>% 
  filter(Datemin_have == "Yes" 
         & Datemin_have_value == min(Datemin_have_value, na.rm = TRUE)
         | Datemin_have == "No" & rank == 1) %>% 
  slice_min(HBA1CDate_min, with_ties = F)

### check for missing value
HBA1C3 %>% 
  ungroup() %>% 
  filter(is.na(HBA1Cvalue)) %>% 
  count() 
#with just slicing on minimum date end up with n=614/5993 without a HBA1C value; 
#and therefore n=5379/5993 with HBA1C values
HBA1C3 %>% 
  ungroup() %>% 
  filter(is.na(HBA1Cvalue) | HBA1CDate_min > 182.5) %>% 
  count()
#with just slicing on minimum date end up with n=1344 without a HBA1C value or 
#date > 182.5
#therefore n=1344/5993 would have a replaced HBA1C value 
#and n=4649/5993 with HBA1C values in the right time frame

#identify median value that should replace missing values with
HBA1C3 %>%
  filter(!is.na(HBA1Cvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(HBA1Cvalue, na.rm = TRUE),
    n = n()
  )
median(HBA1C3$HBA1Cvalue, na.rm = TRUE)
#Median is 54 of n=5390
#this identifies median value in those within 1/2 year of cpep
#use this value
HBA1C3 %>%
  filter(HBA1CDate_min <= 182.5 & !is.na(HBA1Cvalue)) %>%
  ungroup() %>%
  summarise(
    median = median(HBA1Cvalue, na.rm = TRUE),
    n = n()
  )
#Median is 54 of n=4649

### I replaced missing values and those >=182.5 days with median HBA1C 54 
#(median of all HBA1C values within 6 months of C-pep)
HBA1C3 %>%
  filter(HBA1CDate_min >182.5) %>%
  ungroup() %>%
  count()
#replace 730/5993 values outside 6 months
HBA1C3$NewHBA1C <- ifelse((HBA1C3$HBA1CDate_min >182.5), 54, HBA1C3$HBA1Cvalue)
HBA1C3 %>%
  filter(is.na(NewHBA1C)) %>%
  ungroup() %>%
  count()
#replace further 614/5993 missing values
HBA1C3$NewHBA1C <- ifelse(is.na(HBA1C3$NewHBA1C), 54, HBA1C3$NewHBA1C)
summary(HBA1C3$NewHBA1C)
summary(HBA1C3$HBA1Cvalue)


## make remaining variables -----------------------------------------------------
DARE <- DARE %>%
  mutate(
    #the duration between the Goldc-peptide and the recruitment date
    dur_rec_goldcpep = as.numeric(
      difftime(GoldStndCpepDte, DateSeen, units = "weeks")), 
    #age at cpep measuremtns
    age_cpep = as.numeric(
      difftime(GoldStndCpepDte, DateofBirth, units = "weeks"))/52.14, 
    #the difference between age at recruitment and age at cpep
    diff_age_cpep_age = age_cpep - Age, 
    #duration in years between bmi and cpeptide measurement
    dur_bmi_cpep = as.numeric(
      difftime(GoldStndCpepDte, bmi1Date, units = "weeks")/52.14), 
    #duration in years between bmi and recruitment/visit date
    dur_bmi_rec = as.numeric(
      difftime(DateSeen, bmi1Date, units = "weeks")/52.14), 
    #combined duration between bmi date and c-peptide date
    dur_BMI = ifelse(
      !is.na(GoldStndCpepDte), 
      as.numeric(difftime(GoldStndCpepDte, bmi1Date, units = "weeks")/52.14), 
      as.numeric(difftime(DateSeen, bmi1Date, units = "weeks")/52.14)), 
    #duration between creatinine date and c-peptide date
    dur_creatinine = ifelse(
      !is.na(GoldStndCpepDte), 
      as.numeric(difftime(GoldStndCpepDte, DateSerumCreat, units = "weeks")/52.14), 
      as.numeric(difftime(DateSeen, DateSerumCreat, units = "weeks")/52.14)), 
    #combined c-peptide reference date
    ref_date = if_else(!is.na(GoldStndCpepDte),
                       as.Date(GoldStndCpepDte),
                       as.Date(DateSeen)),
    sex = ifelse(sexMale == 1, 0, 1), 
    BMI = bmi, 
    Age_recruitment = Age, 
    #composite age at c-peptide
    Age = ifelse(!is.na(age_cpep), 
                 age_cpep, 
                 Age_recruitment),
    BasalBolus = ifelse(BasalBolus == "Yes", 0, 1),
    OnSU = ifelse(sulphonylurea == 1, 0, 1),
    creatinine = SerumCreat, 
    Cpep900 = ifelse(GoldStndCpepPmolL >= 900, "High", "Low"), 
    Cpep600 = ifelse(GoldStndCpepPmolL >= 600, "High", "Low")
  )

save(DARE, file = "DARE.RData")

# Make dataset with only relevant variables --------------------------------------
DARE_cpep <- DARE %>%
  select(ID, sex, Age, age_cpep, Age_recruitment, cpep_year,DateSeen, dur_rec_goldcpep, 
         BMI, bmi1Date, dur_BMI, dur_bmi_cpep, dur_bmi_rec, OnSU, creatinine, 
         DateSerumCreat, dur_creatinine, GoldStndCpepDte, GoldStndCpepPmolL, 
         Cpep600, Cpep900, GADPositive975, TypeofDiabetes, insulinTreated, 
         BasalBolus,Insulin1Type, Insulin2Type, Insulin3Type, TimeToInsulin, 
         MonthsbetweenDmdiagnosisand, DiabDuration, 
         DurationYrsdiabetes, Durationyearsofdiagnosisat)
#Join in HBA1C
DARE_cpep <- left_join(
  DARE_cpep, 
  HBA1C3[, c("ID", "date", "HBA1CDate_min", "NewHBA1C", "HBA1Cvalue")], 
  by = "ID") %>%
  rename(HBA1C_DATE = date)
#Join in HDL
DARE_cpep <- left_join(
  DARE_cpep, 
  HDL3[, c("ID", "date", "HDLDate_min", "NewHDL", "HDLvalue")], 
  by = "ID") %>%
  rename(HDL_DATE = date)
#Join in triglycerides
DARE_cpep <- left_join(
  DARE_cpep, 
  TG3[, c("ID", "date", "TGDate_min", "NewTG", "trigvalue")], 
  by = "ID") %>%
  rename(TG_DATE = date)
#Join in ALT
DARE_cpep <- left_join(
  DARE_cpep, 
  ALT3[, c("ID", "date", "ALTDate_min", "NewALT", "ALTvalue")], 
  by = "ID") %>%
  rename(ALT_DATE = date)
save(DARE_cpep, file = "DARE_cpep.RData")
