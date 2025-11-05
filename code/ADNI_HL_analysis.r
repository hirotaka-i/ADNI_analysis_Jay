# load library
library(tidyverse)
library(skimr)
library(ADNIMERGE)
library(haven)
library(lme4)
library(lmerTest)
library(emmeans)

df = read.csv('temp/df.csv')

# Baseline DX.bl vs PHC_AMYLOID_STATUS
df %>% filter(month==0) %>%
  with(table(DX.bl, PHC_AMYLOID_STATUS))

### EXPOSURES (BASELINE STATUS) ####
## Use test data
data(adni1lipidomicsrader)
skim(adni1lipidomicsrader) # 792
table(adni1lipidomicsrader$VISCODE, useNA = "ifany") # all bl
id_test_all =  unique(adni1lipidomicsrader$RID)
length(id_test_all) # 773
id_test_200 = adni1lipidomicsrader %>% filter(CHOL > 200) %>%
  pull(RID) %>% unique()
length(id_test_200) # 249

# Use medication and history data
medications <- read.csv("data/Medication/RECCMEDS_21Oct2025.csv")
medical_hist <- read.csv("data/Medical_History/INITHEALTH_21Oct2025.csv")



