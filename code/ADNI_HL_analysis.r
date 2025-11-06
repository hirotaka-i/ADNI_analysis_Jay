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

# medication
med_bl <- medications %>% filter(VISCODE2 %in% c('bl', 'sc'))
id_med_bl_all = unique(med_bl$RID)
length(id_med_bl_all) ##3948
id_med_bl_chol <- med_bl %>%
  mutate(CHOLMED = case_when(
      str_detect(tolower(CMMED), "lipitor|lescol|altoprev|livalo|pravachol|
crestor|zocor|zetia|praluent|repatha|nexletol|nexlizet|prevalite|welchol|colestid|
vytorin|caduet|antara|lipofen|lopid|niacor|niaspan|lovaza|omacor|vascepa|atorvastatin|
fluvastatin|lovastatin|pitavastatin|pravastatin|rosuvastatin|simvastatin|ezetimibe|alirocumab|
evolocumab|bempedoic|cholestyramine|colesevelam|colestipol| 
statin|fenofibrate|gemfibrozil|niacint|cycloset|nesina|tradjenta|glyxambi|onglyza|januvia|trulicity|byetta|
      Bydureon BCise|Saxenda|Victoza") ~ "High",
      TRUE ~ "Normal"
    )
  ) %>% filter(CHOLMED == "High") %>% pull(RID) %>% unique()
length(id_med_bl_chol) ##0 rows with high cholesterol

# medical hist
hist_sc = medical_hist %>% filter(VISCODE2 =='sc')
id_hist_sc_all = unique(hist_sc$RID)
length(id_hist_sc_all) ##1802
id_hist_sc_chol <- medical_hist %>% 
  filter(VISCODE2 =='sc') %>%
  mutate(
    CHOLHIST = case_when(
      str_detect(
        str_to_lower(IHDESC),
        "cholesterol|hyperlipidemia"
      ) ~ "High",
      TRUE ~ "Normal"
    )
  ) %>% filter(CHOLHIST == "High") %>% pull(RID) %>% unique()

length(id_hist_sc_chol) ##872 rows with high cholesterol

# union ids

# combined diabetes exposure at baseline
id_hl = union(id_med_bl_chol, id_hist_sc_chol) %>% union(., id_test_200)
id_all = union(id_med_bl_all, id_hist_sc_all) %>% union(., id_test_all)

df_hl = df %>% filter(RID %in% id_all) %>%
  mutate(HL = if_else(RID %in% id_hl, 1, 0)) %>%
  mutate(HL_G4 = case_when(
    (PHC_AMYLOID_STATUS==1)&(HL == 1) ~ "Abeta yes HL yes",
    (PHC_AMYLOID_STATUS==1)&(HL == 0) ~ "Abeta yes HL no",
    (PHC_AMYLOID_STATUS==0)&(HL == 1) ~ "Abeta no HL yes",
    (PHC_AMYLOID_STATUS==0)&(HL == 0) ~ "Abeta no HL no",
  )) %>% 
  filter(!is.na(HL_G4)) %>%
  mutate(
    HL_G4 = factor(HL_G4, levels = c("Abeta no HL no", "Abeta no HL yes", "Abeta yes HL no", "Abeta yes HL yes")),
    APOEe4_carrier = if_else(APOE4 %in% c("1", "2"), 1, 0),
    year = month / 12
  )

skim(df_hl)

# separate by baseline diagnosis
df_hl %>% filter(month==0) %>%
  with(table(HL_G4, DX.bl,  useNA = "ifany"))

df_hl_ad = df_hl %>% filter(DX.bl == "AD")
df_hl_mci = df_hl %>% filter(DX.bl %in% c("EMCI", "LMCI", 'SMC'))
df_hl_cn = df_hl %>% filter(DX.bl %in% c("CN"))

PACC_HL <- lmer(mPACCdigit ~ year * HL_G4 + AGE + PTGENDER + APOEe4_carrier + PTRACCAT + PTEDUCAT + PTMARRY + (1 | RID), data = df_hl_cn)

summary(PACC_HL)


# Group 4 level comparisons (intecept)
em_intercept <- emmeans(
  PACC_HL,
  ~ HL_G4,
  at = list(year = 0),
  lmerTest.limit = 1e5,
  pbkrtest.limit = 1e5
)
summary(em_intercept)

# test for interaction worse than additive
contrast(em_intercept, list("Interaction worse than additive" = c(1, -1, -1, 1)))
# 1-to-1 comparisons
pairs(em_intercept, adjust = "tukey")  # 或 adjust = "tukey" 做多重比較校正

# Group 4 level comparisons (slope)
em_slope <- emtrends(
  PACC_HL,
  ~ HL_G4,
  var = "year"
)
summary(em_slope)
contrast(em_slope, list("Interaction worse than additive" = c(1, -1, -1, 1)))
pairs(em_slope, adjust = "tukey")
