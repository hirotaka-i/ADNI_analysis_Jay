# load library
library(tidyverse)
library(skimr)
library(ADNIMERGE)
library(haven)

# load adni
data(adnimerge)

dbl <- adnimerge %>% 
  filter(VISCODE == "bl") %>% 
  select(RID, AGE, PTGENDER, PTEDUCAT, PTMARRY, PTETHCAT, PTRACCAT, DX.bl, APOE4)

# check if dbl$RID is unique
length(unique(dbl$RID)) == nrow(dbl)

# check variable types
sapply(dbl,class)
attributes(dbl$DX.bl)
attributes(dbl$AGE)
attributes(dbl$PTGENDER)
attributes(dbl$PTEDUCAT)
attributes(dbl$PTMARRY)
attributes(dbl$PTETHCAT)
attributes(dbl$PTRACCAT)
attributes(dbl$APOE4)

# remove the label and make them flat
dbl <- dbl %>%
  mutate(
    RID = as.numeric(RID),
    AGE = as.numeric(AGE),
    PTGENDER = as.factor(as.character(PTGENDER)),
    PTEDUCAT = as.numeric(PTEDUCAT),
    PTETHCAT = as.factor(as.character(PTETHCAT)),    
    PTRACCAT = as.factor(as.character(PTRACCAT)),
    DX.bl = as.factor(as.character(DX.bl)),
    APOE4 = as.factor(as.numeric(APOE4))
  )

skim(dbl) # some missing in APOE4

# make a complete dataset for baseline variables
dbl_set <- dbl[complete.cases(dbl), ]
# check the ID uniquness
length(unique(dbl_set$RID)) == nrow(dbl_set)
skim(dbl_set) # n rows = 2393

# PET status (should be defined at baseline)
pet <- read.csv("data/Image_Analyses/ADSP_PHC_PET_Amyloid_Simple_18Sep2025.csv") %>%
  select(RID, VISCODE, VISCODE2, PTID, PHC_CENTILOIDS, PHC_AMYLOID_STATUS, PHC_CL_FAIL)
pet %>% with(table(VISCODE, VISCODE2)) # VISCODE2 looks better

pet_bl = pet %>% filter(VISCODE2 == "bl") %>% 
  filter(!is.na(PHC_AMYLOID_STATUS)) %>% 
  select(RID, PHC_AMYLOID_STATUS)
table(pet_bl$PHC_AMYLOID_STATUS)

# combine bl pet with dbl_set
dbl_pet <- inner_join(dbl_set, pet_bl, by = "RID")
skim(dbl_pet) # 1536
dbl_pet %>% with(table(DX.bl, PHC_AMYLOID_STATUS))

# picc outcome - longitudinal
dd <- adnimerge %>%
  select(ADASQ4, LDELTOTAL, DIGITSCOR, MMSE, TRABSCOR, DX.bl, VISCODE, RID, PTID)

adni_pacc <- pacc(dd, keepComponents = FALSE) %>%
  select(RID, VISCODE, mPACCdigit, mPACCtrailsB) %>%
  mutate(
    RID = as.numeric(RID),
    month = if_else(VISCODE == "bl", 0, as.numeric(str_extract(VISCODE, "\\d+")))
  )

# check month and viscode
adni_pacc %>% with(table(VISCODE, month))
# check the uniqueness of RID and month combination
{adni_pacc %>% group_by(RID, month) %>% filter(n() > 1) %>% nrow()} == 0


skim(adni_pacc)

adni_pacc_complete = adni_pacc[complete.cases(adni_pacc), ]

skim(adni_pacc_complete) # n rows = 11508

# combine with baseline variables (dbl_pet)
df = inner_join(adni_pacc_complete, dbl_pet, by = "RID")
df %>% with(table(month, PHC_AMYLOID_STATUS))

skim(df) # nrows 6318

write.csv(df, 'temp/df.csv')
