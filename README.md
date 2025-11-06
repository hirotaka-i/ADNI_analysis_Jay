# Introduction

Repository to repliate Alice's study on A4. 

# Requirements
adnimerge R package

# set-up
```{bash}
mkdir data
mkdir temp
```
Place the adnimerge package in the data folder

# Data Preparation

1. `dbl`: Basic variables from adnimerge baseline data:
   
   complete dataset n= 2,392 (~40 missing APOE4 status)

   Included variables `RID, AGE, PTGENDER, PTEDUCAT, PTMARRY, PTETHCAT, PTRACCAT, DX.bl, APOE4`

2. `dbl_pet`: PET results from `ADSP_PHC_PET_Amyloid_Simple_18Sep2025.csv` n=1,535
```
> dbl_pet %>% with(table(DX.bl, PHC_AMYLOID_STATUS))
      PHC_AMYLOID_STATUS
DX.bl    0   1
  AD    39 178
  CN   229  73
  EMCI 234 171
  LMCI 114 169
  SMC  228 100
```

3. PACC results generated using pacc algorithm: n_obs = 11508

4. `df`: all above results combined: n_obs = 6,317
```
> df %>% filter(month==0) %>%
+   with(table(DX.bl, PHC_AMYLOID_STATUS))
      PHC_AMYLOID_STATUS
DX.bl    0   1
  AD    39 178
  CN   229  73
  EMCI 234 171
  LMCI 114 169
  SMC  228 100
```
So it seems half of the PACC data are dropped because they don't have PET data. Those with PET data seem to have PACC data available. 

# Testing DM
After DM exposure defined by medication and initial health evaluation data, the combined diabetes exposure at baseline was created.
```
> df_dm %>% filter(month==0) %>%
+   with(table(DM_G4, DX.bl,  useNA = "ifany"))
                  DX.bl
DM_G4               AD  CN EMCI LMCI SMC
  Abeta no DM no    33 208  204   98 196
  Abeta no DM yes    6  20   30   16  31
  Abeta yes DM no  159  68  154  155  91
  Abeta yes DM yes  19   5   17   14   9

```
Made thress strata and tested for df_dm_cn
```
df_dm_ad = df_dm %>% filter(DX.bl == "AD")
df_dm_mci = df_dm %>% filter(DX.bl %in% c("EMCI", "LMCI", 'SMC'))
df_dm_cn = df_dm %>% filter(DX.bl %in% c("CN", "SMC"))
```

```
> summary(em_slope)
 DM_G4            year.trend     SE   df lower.CL upper.CL
 Abeta no DM no      -0.0662 0.0262 1849   -0.118  -0.0148
 Abeta no DM yes      0.0186 0.0998 1901   -0.177   0.2143
 Abeta yes DM no     -0.5745 0.0419 1856   -0.657  -0.4924
 Abeta yes DM yes    -0.3260 0.2020 1898   -0.722   0.0696
```
And there is no interaction effect.

# Testing HL
```
 with(table(HL_G4, DX.bl,  useNA = "ifany"))
                  DX.bl
HL_G4               AD  CN EMCI LMCI SMC
  Abeta no HL no    17 118  128   44 117
  Abeta no HL yes   22 110  106   70 110
  Abeta yes HL no   78  40   89   81  45
  Abeta yes HL yes 100  33   82   88  55
```
```
> summary(em_slope)
 HL_G4            year.trend     SE   df lower.CL upper.CL
 Abeta no HL no      -0.0664 0.0394 1105   -0.144   0.0110
 Abeta no HL yes     -0.0640 0.0465 1118   -0.155   0.0272
 Abeta yes HL no     -0.6430 0.0666 1104   -0.774  -0.5123
 Abeta yes HL yes    -0.6784 0.0830 1133   -0.841  -0.5155

Results are averaged over the levels of: PTGENDER, APOEe4_carrier, PTRACCAT, PTMARRY 
Degrees-of-freedom method: kenward-roger 
Confidence level used: 0.95 


> contrast(em_slope, list("Interaction worse than additive" = c(1, -1, -1, 1)))
 contrast                        estimate    SE   df t.ratio p.value
 Interaction worse than additive  -0.0378 0.123 1121  -0.309  0.7577

Results are averaged over the levels of: PTGENDER, APOEe4_carrier, PTRACCAT, PTMARRY 
Degrees-of-freedom method: kenward-roger 
```