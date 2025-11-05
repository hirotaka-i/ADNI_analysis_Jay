# load library
library(tidyverse)
library(lme4)
library(skimr)
library(lmerTest)
library(emmeans)

df = read.csv('temp/df.csv')

# Baseline DX.bl vs PHC_AMYLOID_STATUS
df %>% filter(month==0) %>%
  with(table(DX.bl, PHC_AMYLOID_STATUS))

### EXPOSURES (BASELINE STATUS) ####
## Use medication data and initial health evaluation data to define diabetes exposure at baseline
medications <- read.csv("data/Medication/RECCMEDS_21Oct2025.csv")
medical_hist <- read.csv("data/Medical_History/INITHEALTH_21Oct2025.csv")

# diabetes exposure dataset
med_bl <- medications %>% filter(VISCODE2 %in% c('bl', 'sc'))
id_med_bl_all = unique(med_bl$RID)
length(id_med_bl_all)

dm_stems <- c(
  # Insulins (brand/generic)
  "insul",      # insulin / insuline / insulin glargine etc.
  "humulin",
  "novolin",
  "novolog",
  "humalog",
  "fiasp",
  "lyumjev",
  "basaglar",
  "tresiba",
  "lantus",
  "toujeo",
  "semglee",
  "afrezza",
  
  # Biguanide
  "metfor",     # metformin, metfornin, etc.
  
  # SUs and meglitinides
  "glipizid",
  "glybur",     # glyburide
  "glimepir",   # glimepiride
  "gliclazid",
  "repaglin",   # repaglinide
  "nateglin",   # nateglinide / starlix
  
  # DPP-4 inhibitors
  "sitaglipt",  # sitagliptin
  "aloglipt",
  "linaglipt",  # tradjenta
  "saxaglipt",  # ongliza
  "nesin",      # nesina
  
  # GLP-1 agonists
  "liraglut",   # victoza / saxenda
  "semaglut",   # ozempic, rybelsus
  "dulaglut",   # trulicity
  "exenat",     # byetta / bydureon
  "lixisenat",  # adlyxin
  "tirzepat",   # mounjaro
  
  # SGLT2 inhibitors
  "invokan",    # invokana / invokamet
  "canaglif",   # canagliflozin
  "dapaglif",   # farxiga / xigduo
  "empaglif",   # jardiance / synjardy / trijardy
  "ertuglif",   # steglatro
  
  # Combo products (stems are often enough)
  "xigduo",
  "synjardy",
  "trijardy",
  "invokamet",
  "kazano",
  "actoplus",   # actoplus met
  "glucovance",
  "jentadueto",
  "duetact",
  "oseni",
  
  # TZDs etc.
  "pioglitaz",  # actos
  "rosiglitaz"  # avandia
)
# Build regex like: "\\b(insul[a-z0-9]*|basaglar[a-z0-9]*|...)\\b"
dm_pattern <- paste0(
  "\\b(",
  paste0(dm_stems, "[a-z0-9]*", collapse = "|"),
  ")\\b"
)

id_med_bl_dm <- med_bl %>%
  mutate(
    CMMED_low = str_to_lower(CMMED),
    DMMED = if_else(
      str_detect(CMMED_low, dm_pattern),
      "Yes",
      "No"
    )
  ) %>%
  filter(DMMED == "Yes") %>%
  pull(RID) %>%
  unique()
length(id_med_bl_dm)

# Medical history diabetes exposure
table(medical_hist$VISCODE2) # no bl only sc
hist_sc = medical_hist %>% filter(VISCODE2 =='sc')
id_hist_sc_all = unique(hist_sc$RID)
length(id_hist_sc_all)
id_hist_sc_dm <- medical_hist %>% 
  filter(VISCODE2 =='sc') %>%
  mutate(
    DMHIST = case_when(
      str_detect(
        str_to_lower(IHDESC),
        "diabetes|t2dm|type 2 dm|type ii dm|type 2 diabet|dm 2|dm ii|^dm$"
      ) ~ "Yes",
      TRUE ~ "No"
    )
  ) %>% filter(DMHIST == "Yes") %>% pull(RID) %>% unique()

length(id_hist_sc_dm)


# combined diabetes exposure at baseline
id_dm = union(id_med_bl_dm, id_hist_sc_dm)
id_all = union(id_med_bl_all, id_hist_sc_all)

df_dm = df %>% filter(RID %in% id_all) %>%
  mutate(DM = if_else(RID %in% id_dm, 1, 0)) %>%
  mutate(DM_G4 = case_when(
    (PHC_AMYLOID_STATUS==1)&(DM == 1) ~ "Abeta yes DM yes",
    (PHC_AMYLOID_STATUS==1)&(DM == 0) ~ "Abeta yes DM no",
    (PHC_AMYLOID_STATUS==0)&(DM == 1) ~ "Abeta no DM yes",
    (PHC_AMYLOID_STATUS==0)&(DM == 0) ~ "Abeta no DM no",
  )) %>% 
  filter(!is.na(DM_G4)) %>%
  mutate(
    DM_G4 = factor(DM_G4, levels = c("Abeta no DM no", "Abeta no DM yes", "Abeta yes DM no", "Abeta yes DM yes")),
    APOEe4_carrier = if_else(APOE4 %in% c("1", "2"), 1, 0),
    year = month / 12
  )

skim(df_dm) # nrows 6315 will be in the analysis

# separate by baseline diagnosis
df_dm %>% filter(month==0) %>%
  with(table(DM_G4, DX.bl,  useNA = "ifany"))

df_dm_ad = df_dm %>% filter(DX.bl == "AD")
df_dm_mci = df_dm %>% filter(DX.bl %in% c("EMCI", "LMCI", 'SMC'))
df_dm_cn = df_dm %>% filter(DX.bl %in% c("CN", "SMC"))
#### ANALYSIS ####

PACC_DM <- lmer(mPACCdigit ~ year * DM_G4 + AGE + PTGENDER + APOEe4_carrier + PTRACCAT + PTEDUCAT + PTMARRY + (1 | RID), data = df_dm_cn)

summary(PACC_DM)

# Group 4 level comparisons (intecept)
em_intercept <- emmeans(
  PACC_DM,
  ~ DM_G4,
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
  PACC_DM,
  ~ DM_G4,
  var = "year"
)
summary(em_slope)
contrast(em_slope, list("Interaction worse than additive" = c(1, -1, -1, 1)))
pairs(em_slope, adjust = "tukey")