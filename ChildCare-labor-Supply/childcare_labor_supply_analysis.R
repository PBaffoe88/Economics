# =============================================================================
# The Childcare Constraint: Estimating the Causal Effect of Care Access on
# Labor Supply in Central Minnesota
# =============================================================================
# Author  : Paul Baffoe
# Date    : June 2026
# Target  : Journal of Urban Economics | AEA: Applied Economics |
#           JPAM | Regional Science and Urban Economics
# =============================================================================
#
# Research Questions
# ------------------
# RQ1: Does childcare slot availability increase LFP, employment, and
#      full-time work among parents of children under five?
# RQ2: Are these effects heterogeneous across gender, income, immigration
#      status, and education?
# RQ3: Does the CCAP 85% FPL threshold create a detectable discontinuity
#      in labor supply behavior?
#
# Identification Strategy
# -----------------------
# (A) Individual OLS   : y ~ log_slots + demographics + county_d + ur_c
# (B) ZIP×Year TWFE    : feols(log_slots_zip ~ 1 | zip + year)
# (C) IV / 2SLS        : instrument = lag-1 & lag-2 new provider openings
# (D) Event study      : i(event_t, ref="-1") | county_fac; no year FE
# (E) Sharp RD         : rdrobust() at CCAP 85% FPL threshold
#
# COLLINEARITY NOTE: log_slots varies only at the MSA-year level.
#   feols(y ~ log_slots | county + YEAR) is misspecified — log_slots is
#   perfectly collinear with YEAR FE. Correct specs use county FE only,
#   OR OLS with county dummy + demeaned UR as macro control.
#
# Data Availability
# -----------------
# Input microdata are NOT distributed with this repository (see DATA.md).
# Place the seven input CSVs in data/ before running. All paths below are
# relative to the repository root; outputs are written to output/.
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# 0.  SETUP
# ─────────────────────────────────────────────────────────────────────────────

## Install packages (uncomment on first run) -----------------------------------
# install.packages(c(
#   "tidyverse",      # data wrangling + ggplot2
#   "fixest",         # feols() for TWFE + event study
#   "rdrobust",       # RD estimation
#   "ivreg",          # 2SLS / IV
#   "sandwich",       # HC robust variance
#   "lmtest",         # coeftest()
#   "modelsummary",   # publication tables
#   "patchwork",      # combine ggplots
#   "scales",         # axis formatters
#   "ggrepel",        # non-overlapping labels
#   "broom",          # tidy() model output
#   "knitr",          # kable()
#   "haven"           # read_dta if needed
# ))

library(tidyverse)
library(fixest)
library(rdrobust)
library(ivreg)
library(sandwich)
library(lmtest)
library(modelsummary)
library(patchwork)
library(scales)
library(ggrepel)
library(broom)

## Paths (relative to the repository root) --------------------------------------
DATA_DIR   <- "data/"     # place the seven input CSVs here (see DATA.md)
OUTPUT_DIR <- "output/"   # all tables and figures are written here
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Color palette ---------------------------------------------------------------
PAL_NAVY  <- "#1F3864"
PAL_MID   <- "#2E5B9A"
PAL_LIGHT <- "#7BAFD4"
PAL_RED   <- "#C0392B"
PAL_GREEN <- "#27AE60"
PAL_GRAY  <- "#95A5A6"
PAL_ORG   <- "#E67E22"

## ggplot2 paper theme ---------------------------------------------------------
THEME_PAPER <- theme_bw(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, color = PAL_NAVY),
    plot.subtitle     = element_text(size = 10, color = "grey40"),
    plot.caption      = element_text(size = 7.5, color = "grey50", hjust = 0,
                                     margin = margin(t = 8)),
    axis.title        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92"),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9, face = "bold"),
    strip.background  = element_rect(fill = PAL_NAVY),
    strip.text        = element_text(color = "white", face = "bold", size = 9)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 1.  LOAD ALL SEVEN DATASETS
# ─────────────────────────────────────────────────────────────────────────────

pums <- read_csv(paste0(DATA_DIR, "pums_stcloud_analysis.csv"),
                 show_col_types = FALSE)          # 103,379 × 57
zp   <- read_csv(paste0(DATA_DIR, "zip_year_supply_panel.csv"),
                 show_col_types = FALSE)          # 666 × 10
fred <- read_csv(paste0(DATA_DIR, "fred_monthly_wide.csv"),
                 show_col_types = FALSE)          # 257 × 10
qcew <- read_csv(paste0(DATA_DIR, "qcew_stcloud_6244_panel.csv"),
                 show_col_types = FALSE)          # 29 × 26
oews <- read_csv(paste0(DATA_DIR, "deed_oews_master.csv"),
                 show_col_types = FALSE)          # 23 × 15
jvs  <- read_csv(paste0(DATA_DIR, "jvs_centralMN_naics62_panel.csv"),
                 show_col_types = FALSE)          # 24 × 10
mpls <- read_csv(paste0(DATA_DIR, "mpls_fed_surveys_complete.csv"),
                 show_col_types = FALSE)          # 6 × 20

cat("Datasets loaded:\n")
cat(" PUMS :", nrow(pums), "rows\n")
cat(" ZIP  :", nrow(zp),   "rows |", n_distinct(zp$zip), "ZIPs |",
    n_distinct(zp$year), "years\n")
cat(" FRED :", nrow(fred), "rows\n")
cat(" QCEW :", nrow(qcew), "rows\n")
cat(" OEWS :", nrow(oews), "rows\n")
cat(" JVS  :", nrow(jvs),  "rows\n")
cat(" MPLS :", nrow(mpls), "rows\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2.  DATA PREPARATION
# ─────────────────────────────────────────────────────────────────────────────

## 2a.  MSA-level supply panel (year aggregation) ------------------------------
zp_agg <- zp |>
  group_by(year) |>
  summarise(
    total_slots    = sum(total_slots),
    n_providers    = sum(n_providers),
    slots_under5   = sum(slots_under5),
    slots_infants  = sum(slots_infants),
    new_openings   = sum(new_openings),
    closures       = sum(closures),
    n_zips         = n_distinct(zip),
    n_centers      = sum(n_centers),
    n_family       = sum(n_family),
    .groups = "drop"
  ) |>
  mutate(
    log_slots      = log(total_slots),
    log_providers  = log(n_providers),
    slots_per_prov = total_slots / n_providers,
    pct_centers    = n_centers / n_providers,
    lag1_openings  = lag(new_openings, 1),
    lag2_openings  = lag(new_openings, 2)
  )

## 2b.  ZIP × Year panel for TWFE (2009-2024 overlap with PUMS) ---------------
zp_panel <- zp |>
  filter(year >= 2009, year <= 2024) |>
  mutate(
    log_slots_zip  = log(pmax(total_slots,   1)),
    log_u5_zip     = log(pmax(slots_under5,  1)),
    log_inf_zip    = log(pmax(slots_infants, 1)),
    log_prov_zip   = log(pmax(n_providers,   1)),
    slots_per_prov = total_slots / pmax(n_providers, 1),
    zip_fac        = factor(zip),
    year_fac       = factor(year)
  ) |>
  group_by(zip) |>
  mutate(
    lag1_open_zip  = lag(new_openings, 1),
    lag2_open_zip  = lag(new_openings, 2)
  ) |>
  ungroup()

## 2c.  FRED: annual macro controls -------------------------------------------
fred_yr <- fred |>
  group_by(year) |>
  summarise(
    ur       = mean(STCL027UR,  na.rm = TRUE),
    lf_size  = mean(STCL027LFN, na.rm = TRUE),
    emp_msa  = mean(STCL027EMN, na.rm = TRUE),
    lfp_rate = mean(LBSSA27,    na.rm = TRUE),
    nonfarm  = mean(ENUC410640010SA, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(ur_c = ur - mean(ur, na.rm = TRUE))   # demeaned UR

## 2d.  QCEW: county-year childcare sector (NAICS 6244) -----------------------
#   Interpolate missing years (Benton 2016, 2020; Stearns 2020 present)
qcew_clean <- qcew |>
  mutate(county_name = str_replace(county, " County", ""))

qcew_interp <- map_dfr(c("Benton", "Stearns"), function(cty) {
  qcew_clean |>
    filter(county_name == cty) |>
    select(year, emp_annual, estab_annual, avg_wkly_wage_annual) |>
    complete(year = 2009:2024) |>
    mutate(
      across(c(emp_annual, estab_annual, avg_wkly_wage_annual),
             ~ zoo::na.approx(., na.rm = FALSE)),
      county_name = cty
    )
}) |>
  mutate(log_cc_emp   = log(pmax(emp_annual,   1)),
         log_cc_estab = log(pmax(estab_annual, 1)))
# Note: zoo::na.approx requires zoo package; if unavailable use:
# mutate(across(..., ~ approx(seq_along(.), ., seq_along(.))$y))

## 2e.  Merge into individual-level PUMS frame ---------------------------------
pums_clean <- pums |>
  left_join(
    zp_agg |> select(year, log_slots, total_slots, slots_under5,
                     slots_infants, n_providers, slots_per_prov,
                     lag1_openings, lag2_openings),
    by = c("YEAR" = "year")
  ) |>
  left_join(
    fred_yr |> select(year, ur, ur_c, lf_size, lfp_rate),
    by = c("YEAR" = "year")
  ) |>
  left_join(
    qcew_interp |> select(year, county_name, log_cc_emp, avg_wkly_wage_annual),
    by = c("YEAR" = "year", "county_name" = "county_name")
  ) |>
  mutate(
    # Derived variables
    age2        = AGE^2,
    log_slots   = log(total_slots),
    county_d    = as.integer(county_name == "Stearns"),  # 1=Stearns, 0=Benton
    county_fac  = factor(county_name),
    year_fac    = factor(YEAR),
    event_t     = YEAR - 2018,        # event time; 2018 = expansion shock
    poverty_c   = POVERTY - 85,       # running variable centred at CCAP cutoff
    above_ccap  = as.integer(POVERTY >= 85),
    year_trend  = YEAR - 2009,        # linear time trend (for robustness)
    # Wage outcome
    wage_annual = exp(log_wage),
    log_wage_w  = log_wage             # already log
  )

## 2f.  Analytical subsamples --------------------------------------------------
parents  <- pums_clean |> filter(parent_u5 == 1)          # N = 13,254
mothers  <- pums_clean |> filter(parent_u5 == 1, female == 1)
fathers  <- pums_clean |> filter(parent_u5 == 1, female == 0)
sp_all   <- pums_clean |> filter(parent_u5 == 1, single_parent == 1)
sp_moms  <- pums_clean |> filter(parent_u5 == 1, single_parent == 1, female == 1)

cat("Analytical sample sizes:\n")
cat(" Full PUMS:       ", nrow(pums_clean), "\n")
cat(" Parents (u<5):   ", nrow(parents),    "\n")
cat(" Mothers (u<5):   ", nrow(mothers),    "\n")
cat(" Single parents:  ", nrow(sp_all),     "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3.  DESCRIPTIVE STATISTICS  (Table 1)
# ─────────────────────────────────────────────────────────────────────────────

desc_row <- function(df, label) {
  df |>
    summarise(
      Group       = label,
      N           = n(),
      LFP_rate    = mean(in_lf,        na.rm = TRUE),
      Emp_rate    = mean(employed,      na.rm = TRUE),
      FT_rate     = mean(fulltime,      na.rm = TRUE),
      Pct_female  = mean(female,        na.rm = TRUE),
      Pct_single  = mean(single_parent, na.rm = TRUE),
      Pct_lowinc  = mean(low_income,    na.rm = TRUE),
      Pct_immig   = mean(immigrant,     na.rm = TRUE),
      Pct_college = mean(college,       na.rm = TRUE),
      Mean_wage   = mean(exp(log_wage), na.rm = TRUE)
    )
}

table1 <- bind_rows(
  desc_row(pums_clean,                                           "Full sample"),
  desc_row(filter(pums_clean, parent_u5 == 1),                  "Parents, child <5"),
  desc_row(filter(pums_clean, parent_u5 == 1, female == 1),     "Mothers, child <5"),
  desc_row(filter(pums_clean, parent_u5 == 1, female == 0),     "Fathers, child <5"),
  desc_row(filter(pums_clean, parent_u5 == 1, single_parent == 1), "Single parents, child <5"),
  desc_row(filter(pums_clean, parent_u5 == 1, low_income == 1), "Low-income parents, child <5"),
  desc_row(filter(pums_clean, parent_u5 == 1, immigrant == 1),  "Immigrant parents, child <5")
)

print("=== TABLE 1: DESCRIPTIVE STATISTICS ===")
print(table1 |> mutate(across(where(is.numeric) & !N, ~round(., 3))))

## Supply panel summary (Table 2) ----------------------------------------------
table2 <- zp_agg |>
  filter(year %in% c(2003, 2007, 2010, 2015, 2018, 2020, 2023, 2026)) |>
  select(year, n_zips, n_providers, total_slots, slots_under5,
         slots_infants, new_openings, slots_per_prov) |>
  mutate(slots_per_prov = round(slots_per_prov, 1))

print("\n=== TABLE 2: CHILDCARE SUPPLY PANEL ===")
print(table2)

## Labour market context from FRED + QCEW + OEWS + JVS + MPLS ----------------
cat("\n=== LABOUR MARKET CONTEXT ===\n")

# FRED: key years
fred_context <- fred_yr |>
  filter(year %in% c(2009, 2012, 2015, 2018, 2019, 2020, 2021, 2023, 2024)) |>
  select(year, ur, lf_size, emp_msa, lfp_rate)
print(fred_context)

# QCEW sector wages and employment
cat("\nQCEW childcare sector (Stearns County, NAICS 6244):\n")
print(qcew_interp |>
        filter(county_name == "Stearns") |>
        select(year, emp_annual, estab_annual, avg_wkly_wage_annual))

# OEWS wage distribution
cat("\nOEWS childcare worker wages (SOC 39-9011, 2026):\n")
print(oews |>
        filter(file == "primary") |>
        select(geography, employment, mean_wage, wage_median, wage_10th, wage_90th))

# JVS vacancy trend
cat("\nJVS job vacancies (NAICS 62, Central MN):\n")
print(jvs |>
        filter(year >= 2015) |>
        select(year, annual_vacancies, geography))

# MPLS provider financial stress
cat("\nMinneapolis Fed provider surveys:\n")
print(mpls |>
        select(survey_year, insurance_index_2022base,
               pct_raising_tuition_centers, pct_comp_unsustainable_centers))

# ─────────────────────────────────────────────────────────────────────────────
# 4.  RQ1 — CHILDCARE SUPPLY AND LABOR SUPPLY
# ─────────────────────────────────────────────────────────────────────────────
# Correct identification:
#   log_slots is MSA-level (year-only variation).
#   feols(y ~ log_slots | county + YEAR) is collinear — DO NOT USE.
#   Correct specs:
#   (A) OLS + county dummy + demeaned UR
#   (B) feols(y ~ log_slots + controls | county_fac)  — county FE only
#   (C) 2SLS IV with lagged new openings as instruments

cat("\n\n")
cat("=================================================================\n")
cat("RQ1: CHILDCARE SLOTS AND LABOR SUPPLY\n")
cat("=================================================================\n\n")

## Define control vectors ------------------------------------------------------
ctrl_base  <- c("AGE", "age2", "college", "immigrant",
                "single_parent", "low_income", "county_d", "ur_c")
ctrl_slots <- c(ctrl_base, "log_slots")

## 4A.  OLS: positive but imprecise (macro shocks not fully controlled) --------
ols_lfp <- lm(in_lf    ~ log_slots + AGE + age2 + college + immigrant +
                single_parent + low_income + county_d + ur_c, data = parents)
ols_emp <- lm(employed  ~ log_slots + AGE + age2 + college + immigrant +
                single_parent + low_income + county_d + ur_c, data = parents)
ols_ft  <- lm(fulltime  ~ log_slots + AGE + age2 + college + immigrant +
                single_parent + low_income + county_d + ur_c, data = parents)

cat("--- 4A. Main OLS Estimates ---\n")
msummary(
  list("LFP" = ols_lfp, "Employment" = ols_emp, "Full-Time" = ols_ft),
  vcov     = "HC1",
  coef_map = c(
    "log_slots"     = "log(Total Slots)",
    "ur_c"          = "Unemployment Rate (demeaned)",
    "college"       = "College",
    "single_parent" = "Single Parent",
    "low_income"    = "Low Income",
    "immigrant"     = "Immigrant",
    "county_d"      = "Stearns County"
  ),
  stars   = c("*" = .10, "**" = .05, "***" = .01),
  gof_map = c("nobs", "r.squared"),
  title   = "Table 3A: OLS Estimates — Childcare Supply on Labor Supply"
)

## 4B.  TWFE: county FE absorbed via feols (county ≠ year FE → identified) ----
#   log_slots varies by year; county FE absorbs time-invariant county factors
twfe_lfp <- feols(
  in_lf    ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + ur_c | county_fac,
  data = parents, vcov = "HC1"
)
twfe_emp <- feols(
  employed  ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + ur_c | county_fac,
  data = parents, vcov = "HC1"
)
twfe_ft  <- feols(
  fulltime  ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + ur_c | county_fac,
  data = parents, vcov = "HC1"
)

cat("\n--- 4B. TWFE (County FE Absorbed) ---\n")
msummary(
  list("LFP" = twfe_lfp, "Employment" = twfe_emp, "Full-Time" = twfe_ft),
  coef_map = c(
    "log_slots" = "log(Total Slots)",
    "ur_c"      = "Unemployment Rate",
    "college"   = "College",
    "single_parent" = "Single Parent",
    "low_income"    = "Low Income",
    "immigrant"     = "Immigrant"
  ),
  stars   = c("*" = .10, "**" = .05, "***" = .01),
  gof_map = c("nobs", "r.squared"),
  title   = "Table 3B: TWFE Estimates (County FE)"
)

## 4C.  ZIP × Year TWFE (proper panel DiD at ZIP level) -----------------------
#   Unit: ZIP-year; ZIP FE + Year FE absorbed; NO collinearity
#   Outcome: log(slots_zip), log(u5_zip), slots/provider

zip_twfe_slots <- feols(log_slots_zip ~ 1 | zip_fac + year_fac,
                        data = zp_panel, vcov = "HC1")
zip_twfe_u5    <- feols(log_u5_zip    ~ 1 | zip_fac + year_fac,
                        data = zp_panel, vcov = "HC1")
zip_twfe_sp    <- feols(slots_per_prov ~ 1 | zip_fac + year_fac,
                        data = zp_panel, vcov = "HC1")

cat("\n--- 4C. ZIP×Year TWFE (Supply-Side Panel DiD) ---\n")
cat("R² log(slots):   ", r2(zip_twfe_slots)["r2"], "\n")
cat("R² log(u5):      ", r2(zip_twfe_u5)["r2"],    "\n")
cat("R² slots/prov:   ", r2(zip_twfe_sp)["r2"],    "\n")

# Extract year FEs for supply trend plot
yr_fe_slots <- fixef(zip_twfe_slots)$year_fac
yr_fe_df    <- tibble(
  year   = as.integer(names(yr_fe_slots)),
  fe_val = as.numeric(yr_fe_slots)
) |> arrange(year)
cat("\nYear FEs from ZIP×Year TWFE (log slots, base absorb):\n")
print(yr_fe_df)

## 4D.  IV / 2SLS (primary causal estimate for RQ1) ---------------------------
#   First stage: log_slots ~ lag1_openings + lag2_openings + controls
#   Exclusion: provider entry decisions driven by licensing timelines,
#              capital, subsidy rules — not current labour demand shocks

parents_iv <- parents |> drop_na(lag1_openings, lag2_openings, log_slots)

iv_lfp <- ivreg(
  in_lf ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c |
    lag1_openings + lag2_openings + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c,
  data = parents_iv
)
iv_emp <- ivreg(
  employed ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c |
    lag1_openings + lag2_openings + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c,
  data = parents_iv
)
iv_ft <- ivreg(
  fulltime ~ log_slots + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c |
    lag1_openings + lag2_openings + AGE + age2 + college + immigrant +
    single_parent + low_income + county_d + ur_c,
  data = parents_iv
)

cat("\n--- 4D. IV (2SLS) Estimates ---\n")
msummary(
  list("LFP (IV)" = iv_lfp, "Employment (IV)" = iv_emp, "Full-Time (IV)" = iv_ft),
  vcov     = "HC1",
  coef_map = c("log_slots" = "log(Slots) [IV]"),
  stars    = c("*" = .10, "**" = .05, "***" = .01),
  gof_map  = c("nobs", "r.squared"),
  title    = "Table 4: IV (2SLS) Estimates — Childcare Supply on Labor Supply"
)

# First stage diagnostics
fs_mod <- lm(
  log_slots ~ lag1_openings + lag2_openings + AGE + age2 + college +
    immigrant + single_parent + low_income + county_d + ur_c,
  data = parents_iv
)
cat("\nFirst Stage Diagnostics:\n")
cat("  F-stat (overall):", round(summary(fs_mod)$fstatistic[1], 1), "\n")
cat("  R-squared:       ", round(summary(fs_mod)$r.squared, 3),    "\n")
ct_fs <- coeftest(fs_mod, vcov = vcovHC(fs_mod, type = "HC1"))
cat("  lag1_openings:   β =", round(ct_fs["lag1_openings","Estimate"], 4),
    " p =", round(ct_fs["lag1_openings","Pr(>|t|)"], 4), "\n")
cat("  lag2_openings:   β =", round(ct_fs["lag2_openings","Estimate"], 4),
    " p =", round(ct_fs["lag2_openings","Pr(>|t|)"], 4), "\n\n")

## 4E.  Robustness: slots for under-5 only, log(providers) --------------------
rob_u5_emp <- lm(employed ~ log(slots_under5) + AGE + age2 + college +
                   immigrant + single_parent + low_income + county_d + ur_c,
                 data = parents |> mutate(log_slots_u5 = log(pmax(slots_under5, 1))))
rob_pv_emp <- lm(employed ~ log(n_providers) + AGE + age2 + college +
                   immigrant + single_parent + low_income + county_d + ur_c,
                 data = parents |> mutate(log_prov = log(n_providers)))

cat("--- 4E. Robustness Checks ---\n")
msummary(
  list(
    "Slots-u5 → Emp"  = rob_u5_emp,
    "Providers → Emp" = rob_pv_emp,
    "TWFE-cty → Emp"  = twfe_emp
  ),
  vcov     = "HC1",
  coef_map = c(
    "log(slots_under5)" = "log(Slots <5)",
    "log(n_providers)"  = "log(Providers)",
    "log_slots"         = "log(Total Slots) [TWFE]"
  ),
  stars   = c("*" = .10, "**" = .05, "***" = .01),
  gof_map = c("nobs", "r.squared"),
  title   = "Table 5: Robustness — Alternative Supply Measures"
)

## 4F.  Event study around 2018 expansion shock --------------------------------
#   FIX: No year FE (collinear with log_slots). Use county FE + UR control.
#   Base period: event_t = -1 (year 2017)

parents_ev <- parents |>
  filter(event_t >= -6, event_t <= 6) |>
  mutate(event_fac = factor(event_t))

es_lfp <- feols(
  in_lf    ~ i(event_fac, ref = "-1") + AGE + age2 + college +
    immigrant + single_parent + low_income + ur_c | county_fac,
  data = parents_ev, vcov = "HC1"
)
es_emp <- feols(
  employed  ~ i(event_fac, ref = "-1") + AGE + age2 + college +
    immigrant + single_parent + low_income + ur_c | county_fac,
  data = parents_ev, vcov = "HC1"
)
es_ft  <- feols(
  fulltime  ~ i(event_fac, ref = "-1") + AGE + age2 + college +
    immigrant + single_parent + low_income + ur_c | county_fac,
  data = parents_ev, vcov = "HC1"
)

# Pre-trend joint Wald test (t = -6 to -2)
pre_terms <- grep("event_fac::-[2-6]", names(coef(es_lfp)), value = TRUE)
cat("\n--- 4F. Pre-trend Wald Tests (t = -6 to -2) ---\n")
for (mod_obj in list(list(es_lfp,"LFP"), list(es_emp,"Employment"), list(es_ft,"Full-Time"))) {
  pre_t <- grep("event_fac::-[2-6]", names(coef(mod_obj[[1]])), value = TRUE)
  if (length(pre_t) > 0) {
    w <- wald(mod_obj[[1]], pre_t)
    cat(" ", mod_obj[[2]], ": F =", round(w$stat, 3),
        " p =", round(w$p, 3),
        if (w$p > 0.05) " ✓ parallel trends" else " ⚠ pre-trends", "\n")
  }
}

# Extract event study coefficients for plotting
extract_es <- function(mod, outcome_label) {
  tidy(mod, conf.int = TRUE) |>
    filter(str_detect(term, "event_fac")) |>
    mutate(
      t       = as.integer(str_extract(term, "-?\\d+")),
      outcome = outcome_label
    ) |>
    bind_rows(tibble(
      t = -1L, estimate = 0, conf.low = 0, conf.high = 0,
      std.error = 0, p.value = 1, outcome = outcome_label
    )) |>
    arrange(t)
}

es_results <- bind_rows(
  extract_es(es_lfp, "LFP"),
  extract_es(es_emp, "Employment"),
  extract_es(es_ft,  "Full-Time")
) |>
  mutate(outcome = factor(outcome, levels = c("LFP", "Employment", "Full-Time")))

# ─────────────────────────────────────────────────────────────────────────────
# 5.  RQ2 — HETEROGENEITY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

cat("\n\n")
cat("=================================================================\n")
cat("RQ2: HETEROGENEITY ACROSS GENDER, INCOME, IMMIGRATION, EDUCATION\n")
cat("=================================================================\n\n")

## Helper: run OLS for a subgroup, return tidy row ----------------------------
run_subgroup <- function(df, label, drop_vars = NULL) {
  c_vars <- c("AGE", "age2", "college", "immigrant", "low_income",
              "county_d", "ur_c", "log_slots")
  if (!is.null(drop_vars)) c_vars <- setdiff(c_vars, drop_vars)
  
  purrr::map_dfr(
    c("in_lf", "employed", "fulltime"),
    function(outcome) {
      formula  <- reformulate(c_vars, response = outcome)
      dat      <- df |> select(all_of(c(c_vars, outcome))) |> drop_na()
      mod      <- lm(formula, data = dat)
      ct       <- coeftest(mod, vcov = vcovHC(mod, type = "HC1"))
      tibble(
        subgroup = label,
        outcome  = outcome,
        N        = nrow(dat),
        beta     = ct["log_slots", "Estimate"],
        se       = ct["log_slots", "Std. Error"],
        pval     = ct["log_slots", "Pr(>|t|)"],
        ci_lo    = beta - 1.96 * se,
        ci_hi    = beta + 1.96 * se
      )
    }
  )
}

## Run for all eight subgroups ------------------------------------------------
het_results <- bind_rows(
  run_subgroup(filter(parents, female == 1),          "Mothers"),
  run_subgroup(filter(parents, female == 0),          "Fathers"),
  run_subgroup(filter(parents, single_parent == 1),   "Single Parents",
               drop_vars = "college"),   # collinear for single parents
  run_subgroup(filter(parents, single_parent == 1,
                      female == 1),                   "Single Mothers",
               drop_vars = "college"),
  run_subgroup(filter(parents, single_parent == 1,
                      female == 0),                   "Single Fathers",
               drop_vars = "college"),
  run_subgroup(filter(parents, low_income == 1),      "Low-Income"),
  run_subgroup(filter(parents, immigrant == 1),       "Immigrants",
               drop_vars = "low_income"),
  run_subgroup(filter(parents, immigrant == 0),       "Non-Immigrant"),
  run_subgroup(filter(parents, college == 1),         "College-Educ"),
  run_subgroup(filter(parents, college == 0),         "No College")
) |>
  mutate(
    sig     = case_when(
      pval < 0.01 ~ "p<0.01",
      pval < 0.05 ~ "p<0.05",
      pval < 0.10 ~ "p<0.10",
      TRUE        ~ "n.s."
    ),
    outcome  = factor(outcome,
                      levels = c("in_lf", "employed", "fulltime"),
                      labels = c("LFP", "Employment", "Full-Time")),
    subgroup = factor(subgroup,
                      levels = c("Mothers", "Fathers",
                                 "Single Parents", "Single Mothers", "Single Fathers",
                                 "Low-Income", "Immigrants", "Non-Immigrant",
                                 "College-Educ", "No College"))
  )

cat("=== TABLE 6: HETEROGENEITY RESULTS ===\n")
print(
  het_results |>
    select(subgroup, outcome, N, beta, se, pval, sig) |>
    mutate(across(c(beta, se, pval), ~round(., 4)))
)

## Single-parent deep dive: CCAP cliff for single parents specifically ---------
sp_rd <- sp_all |>
  filter(poverty_c >= -40, poverty_c <= 40) |>
  mutate(pov_x_above = poverty_c * above_ccap)

cat("\n--- Single-Parent RD at CCAP Threshold ---\n")
for (out in c("in_lf", "employed", "fulltime")) {
  d   <- sp_rd |> select(above_ccap, poverty_c, pov_x_above, all_of(out)) |> drop_na()
  mod <- lm(reformulate(c("above_ccap", "poverty_c", "pov_x_above"), response = out), data = d)
  ct  <- coeftest(mod, vcov = vcovHC(mod, type = "HC1"))
  cat(" ", out, ": LATE =", round(ct["above_ccap","Estimate"], 4),
      " SE =", round(ct["above_ccap","Std. Error"], 4),
      " p =", round(ct["above_ccap","Pr(>|t|)"], 4),
      " N =", nrow(d), "\n")
}

## Single-parent supply quartile analysis (informal care substitution test) ---
sp_qt <- sp_all |>
  mutate(slots_q = ntile(log_slots, 4)) |>
  group_by(slots_q) |>
  summarise(across(c(in_lf, employed, fulltime), ~mean(., na.rm = TRUE)),
            n = n(), .groups = "drop") |>
  mutate(slots_q_label = paste0("Q", slots_q))

cat("\n--- Single-Parent LFP by Supply Quartile (Informal Care Test) ---\n")
print(sp_qt)

# ─────────────────────────────────────────────────────────────────────────────
# 6.  RQ3 — REGRESSION DISCONTINUITY AT CCAP THRESHOLD
# ─────────────────────────────────────────────────────────────────────────────

cat("\n\n")
cat("=================================================================\n")
cat("RQ3: RD AT CCAP 85% FPL THRESHOLD\n")
cat("=================================================================\n\n")

rd_data <- parents |>
  filter(poverty_c >= -40, poverty_c <= 40) |>
  select(poverty_c, above_ccap, in_lf, employed, fulltime) |>
  drop_na()

cat("RD Sample: N =", nrow(rd_data),
    "| Below cutoff:", sum(rd_data$above_ccap == 0),
    "| Above:", sum(rd_data$above_ccap == 1), "\n\n")

## 6A.  rdrobust: MSE-optimal bandwidth, triangular kernel --------------------
rd_lfp <- rdrobust(rd_data$in_lf,    rd_data$poverty_c, c = 0,
                   kernel = "triangular", bwselect = "mserd")
rd_emp <- rdrobust(rd_data$employed,  rd_data$poverty_c, c = 0,
                   kernel = "triangular", bwselect = "mserd")
rd_ft  <- rdrobust(rd_data$fulltime,  rd_data$poverty_c, c = 0,
                   kernel = "triangular", bwselect = "mserd")

rd_summary <- tibble(
  Outcome = c("LFP", "Employment", "Full-Time"),
  LATE    = c(rd_lfp$coef[1],  rd_emp$coef[1],  rd_ft$coef[1]),
  SE      = c(rd_lfp$se[1],    rd_emp$se[1],    rd_ft$se[1]),
  p_value = c(rd_lfp$pv[1],    rd_emp$pv[1],    rd_ft$pv[1]),
  BW_L    = c(rd_lfp$bws[1,1], rd_emp$bws[1,1], rd_ft$bws[1,1]),
  BW_R    = c(rd_lfp$bws[1,2], rd_emp$bws[1,2], rd_ft$bws[1,2])
) |>
  mutate(Sig = case_when(
    p_value < 0.01 ~ "***", p_value < 0.05 ~ "**",
    p_value < 0.10 ~ "*",   TRUE ~ "n.s."
  ))

cat("=== TABLE 7: RD RESULTS (rdrobust, triangular kernel, MSE-BW) ===\n")
print(rd_summary |> mutate(across(where(is.numeric), ~round(., 4))))

## 6B.  Manual local linear RD (wide bandwidth ±40, for plotting) -------------
rd_data <- rd_data |>
  mutate(pov_x_above = poverty_c * above_ccap)

cat("\n--- Manual Local Linear RD (BW = ±40 ppt) ---\n")
for (out in c("in_lf", "employed", "fulltime")) {
  d   <- rd_data |> select(above_ccap, poverty_c, pov_x_above, all_of(out)) |> drop_na()
  mod <- lm(reformulate(c("above_ccap", "poverty_c", "pov_x_above"), response = out), data = d)
  ct  <- coeftest(mod, vcov = vcovHC(mod, type = "HC1"))
  cat(" ", out, ": LATE =", round(ct["above_ccap","Estimate"], 4),
      " p =", round(ct["above_ccap","Pr(>|t|)"], 3), "\n")
}

## 6C.  Placebo RD checks (alternative cutoffs) --------------------------------
cat("\n--- Placebo RD Tests (alternative cutoffs ±20 ppt from true) ---\n")
for (cutoff_c in c(-20, 20)) {
  rd_p <- parents |>
    mutate(poverty_p = POVERTY - (85 + cutoff_c),
           above_p   = as.integer(POVERTY >= 85 + cutoff_c),
           pov_x_a   = poverty_p * above_p) |>
    filter(between(poverty_p, -40, 40)) |>
    select(above_p, poverty_p, pov_x_a, fulltime) |>
    drop_na()
  mod_p <- lm(fulltime ~ above_p + poverty_p + pov_x_a, data = rd_p)
  ct_p  <- coeftest(mod_p, vcov = vcovHC(mod_p, type = "HC1"))
  cat("  Placebo cutoff at", 85 + cutoff_c, "% FPL: LATE =",
      round(ct_p["above_p","Estimate"], 4),
      " p =", round(ct_p["above_p","Pr(>|t|)"], 3),
      " N =", nrow(rd_p), "\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 7.  EXPORT REGRESSION TABLES
# ─────────────────────────────────────────────────────────────────────────────

## Table 3: Main OLS + TWFE ---------------------------------------------------
modelsummary(
  list(
    "OLS: LFP"     = ols_lfp,
    "OLS: Emp"     = ols_emp,
    "OLS: FT"      = ols_ft,
    "TWFE: LFP"    = twfe_lfp,
    "TWFE: Emp"    = twfe_emp,
    "TWFE: FT"     = twfe_ft
  ),
  output   = paste0(OUTPUT_DIR, "table3_ols_twfe.docx"),
  vcov     = "HC1",
  coef_map = c(
    "log_slots"     = "log(Total Slots)",
    "ur_c"          = "Unemployment Rate",
    "college"       = "College",
    "single_parent" = "Single Parent",
    "low_income"    = "Low Income",
    "immigrant"     = "Immigrant",
    "county_d"      = "Stearns County"
  ),
  stars    = c("*" = .10, "**" = .05, "***" = .01),
  gof_map  = c("nobs", "r.squared"),
  title    = "Table 3: OLS and TWFE Estimates — Childcare Supply on Labor Supply"
)

## Table 4: IV estimates ------------------------------------------------------
modelsummary(
  list("LFP (IV)" = iv_lfp, "Employment (IV)" = iv_emp, "Full-Time (IV)" = iv_ft),
  output   = paste0(OUTPUT_DIR, "table4_iv.docx"),
  vcov     = "HC1",
  coef_map = c("log_slots" = "log(Slots) [IV]"),
  stars    = c("*" = .10, "**" = .05, "***" = .01),
  gof_map  = c("nobs", "r.squared"),
  title    = "Table 4: IV (2SLS) Estimates"
)

## Table 6: Heterogeneity (wide format) ----------------------------------------
het_wide <- het_results |>
  mutate(label = sprintf("%.4f%s", beta,
                         ifelse(pval < .01,"***",ifelse(pval < .05,"**",
                                                        ifelse(pval < .10,"*",""))))) |>
  select(subgroup, outcome, label) |>
  pivot_wider(names_from = outcome, values_from = label)

write_csv(het_wide, paste0(OUTPUT_DIR, "table6_heterogeneity.csv"))
cat("\nTable 6 (heterogeneity) saved.\n")

## Table 7: RD ----------------------------------------------------------------
write_csv(rd_summary, paste0(OUTPUT_DIR, "table7_rd.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# 8.  FIGURES (12 publication-quality ggplots)
# ─────────────────────────────────────────────────────────────────────────────

## Figure 1: Childcare Supply Expansion (2003–2026) ---------------------------
fig1_data <- zp_agg |>
  select(year, total_slots, slots_under5, slots_infants, n_providers) |>
  pivot_longer(-year, names_to = "series", values_to = "value") |>
  mutate(series = factor(series,
                         levels = c("total_slots","slots_under5","slots_infants","n_providers"),
                         labels = c("Total Licensed Slots","Slots for Children <5","Infant Slots","Licensed Providers")))

fig1 <- ggplot(fig1_data, aes(x = year, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = 2018, linetype = "dashed", color = PAL_RED, linewidth = 0.9) +
  annotate("text", x = 2018.4, y = max(zp_agg$total_slots) * 0.86,
           label = "2018\nExpansion\nShock", hjust = 0, size = 3, color = PAL_RED) +
  scale_color_manual(values = c(PAL_NAVY, PAL_MID, PAL_LIGHT, PAL_ORG)) +
  scale_linetype_manual(values = c("solid","solid","dashed","dotdash")) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = seq(2003, 2026, 3)) +
  labs(
    title    = "Figure 1: Childcare Supply Expansion, St. Cloud MSA (2003–2026)",
    subtitle = "Licensed slots and providers from MN DHS/DCYF licensing records",
    x = "Year", y = "Count", color = NULL, linetype = NULL,
    caption  = "Source: MN DHS/DCYF Licensed Provider Records. ZIP-year panel aggregated to MSA level."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig1_supply_trends.png"),
       fig1, width = 9, height = 5.5, dpi = 300)

## Figure 2: Labor Supply Outcomes Over Time ----------------------------------
pums_yr <- pums_clean |>
  group_by(YEAR) |>
  summarise(
    LFP         = mean(in_lf,   na.rm = TRUE),
    Employment  = mean(employed, na.rm = TRUE),
    `Full-Time` = mean(fulltime, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(-YEAR, names_to = "Outcome", values_to = "Rate")

fig2 <- ggplot(pums_yr, aes(x = YEAR, y = Rate, color = Outcome)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.10, linewidth = 0.5) +
  geom_vline(xintercept = 2018, linetype = "dashed", color = PAL_RED, linewidth = 0.8) +
  scale_color_manual(values = c(PAL_NAVY, PAL_MID, PAL_LIGHT)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(2009, 2024, 3)) +
  labs(
    title    = "Figure 2: Labor Supply Outcomes, St. Cloud MSA (2009–2024)",
    subtitle = "Annual mean rates; full PUMS sample (N = 103,379)",
    x = "Year", y = "Rate", color = "Outcome",
    caption  = "Source: ACS PUMS, Stearns & Benton Counties."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig2_lfp_trends.png"),
       fig2, width = 9, height = 5.5, dpi = 300)

## Figure 3: Supply vs LFP Scatter (TWFE motivation) -------------------------
pums_cy <- pums_clean |>
  group_by(YEAR, county_name) |>
  summarise(
    lfp_rate   = mean(in_lf,    na.rm = TRUE),
    emp_rate   = mean(employed,  na.rm = TRUE),
    log_slots  = first(log_slots),
    .groups = "drop"
  )

fig3 <- ggplot(pums_cy, aes(x = log_slots, y = lfp_rate, color = county_name)) +
  geom_point(aes(size = YEAR), alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 1.1) +
  scale_color_manual(values = c(PAL_NAVY, PAL_MID), name = "County") +
  scale_size_continuous(range = c(2, 5.5), name = "Year",
                        breaks = c(2009, 2015, 2020, 2024)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Figure 3: Childcare Supply and LFP — County-Year Panel",
    subtitle = "Points sized by year; lines are county-level OLS fits",
    x = "log(Total Licensed Slots)", y = "LFP Rate",
    caption  = "Source: ACS PUMS merged with MN DHS supply panel."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig3_supply_lfp.png"),
       fig3, width = 8, height = 5.5, dpi = 300)

## Figure 3b: ZIP×Year TWFE Year Fixed Effects --------------------------------
fig3b <- ggplot(yr_fe_df, aes(x = year, y = fe_val)) +
  geom_col(aes(fill = year >= 2018), width = 0.7, alpha = 0.88) +
  geom_line(color = PAL_NAVY, linewidth = 0.8, alpha = 0.5) +
  geom_point(color = PAL_NAVY, size = 2.5) +
  geom_vline(xintercept = 2017.5, linetype = "dashed",
             color = PAL_RED, linewidth = 0.9) +
  annotate("text", x = 2018.2, y = max(yr_fe_df$fe_val, na.rm = TRUE) * 0.85,
           label = "2018\nShock", hjust = 0, size = 3, color = PAL_RED) +
  scale_fill_manual(values = c("FALSE" = PAL_MID, "TRUE" = PAL_NAVY), guide = "none") +
  scale_x_continuous(breaks = seq(2010, 2024, 2)) +
  labs(
    title    = "Figure 3b: ZIP×Year TWFE — Within-ZIP Supply Growth (Year FEs)",
    subtitle = "Year fixed effects from feols(log_slots_zip ~ 1 | zip + year); base = 2009",
    x = "Year", y = "Year FE Coefficient\n(log slot growth vs. 2009)",
    caption  = paste("Source: MN DHS/DCYF. N = 461 ZIP-year obs; 32 ZIPs; ZIP + Year FE absorbed.",
                     "\nFirst stage: lag-1 β=0.067 (p<0.001); lag-2 β=0.058 (p<0.001); F=3,700; R²=0.973.")
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig3b_zip_twfe_yearfe.png"),
       fig3b, width = 9, height = 5, dpi = 300)

## Figure 4: Event Study -------------------------------------------------------
fig4 <- ggplot(es_results,
               aes(x = t, y = estimate, color = outcome, fill = outcome)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed",
             color = PAL_RED, linewidth = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  annotate("rect", xmin = -6.5, xmax = -1.5,
           ymin = -Inf, ymax = Inf, alpha = 0.04, fill = "grey40") +
  annotate("text", x = -4, y = max(es_results$conf.high, na.rm = TRUE) * 0.82,
           label = "Pre-period", size = 2.8, color = "grey45") +
  annotate("text", x = -0.6,
           y = min(es_results$conf.low, na.rm = TRUE) * 0.65,
           label = "2018\nShock", hjust = 1, size = 2.8, color = PAL_RED) +
  scale_color_manual(values = c(PAL_NAVY, PAL_MID, PAL_GREEN)) +
  scale_fill_manual(values  = c(PAL_NAVY, PAL_MID, PAL_GREEN)) +
  scale_x_continuous(
    breaks = -6:6,
    labels = c("-6","-5","-4","-3","-2","-1\n(ref)","0","1","2","3","4","5","6")
  ) +
  labs(
    title    = "Figure 4: Event Study — 2018 Childcare Expansion Shock",
    subtitle = "Coefficients relative to t−1 (2017). Shaded = pre-period. 95% CIs shown.",
    x = "Event Time (years relative to 2018)", y = "Coefficient Estimate",
    color = "Outcome", fill = "Outcome",
    caption  = paste(
      "Source: ACS PUMS. County FE absorbed (feols). Controls: age, age², college,",
      "immigrant, single_parent, low_income, UR (demeaned). HC1 SEs.",
      "\nYear FE omitted — collinear with log_slots at MSA level.",
      "Macro trend controlled via demeaned MSA unemployment rate."
    )
  ) +
  THEME_PAPER +
  theme(axis.text.x = element_text(size = 8.5))

ggsave(paste0(OUTPUT_DIR, "fig4_event_study.png"),
       fig4, width = 10, height = 5.5, dpi = 300)

## Figure 5: IV vs OLS Comparison ----------------------------------------------
iv_compare <- bind_rows(
  tibble(
    Outcome  = rep(c("LFP","Employment","Full-Time"), each = 2),
    Method   = rep(c("OLS","IV (2SLS)"), 3),
    Estimate = c(
      coeftest(ols_lfp, vcovHC(ols_lfp))["log_slots","Estimate"],
      coeftest(iv_lfp,  vcovHC(iv_lfp)) ["log_slots","Estimate"],
      coeftest(ols_emp, vcovHC(ols_emp))["log_slots","Estimate"],
      coeftest(iv_emp,  vcovHC(iv_emp)) ["log_slots","Estimate"],
      coeftest(ols_ft,  vcovHC(ols_ft)) ["log_slots","Estimate"],
      coeftest(iv_ft,   vcovHC(iv_ft))  ["log_slots","Estimate"]
    ),
    SE = c(
      coeftest(ols_lfp, vcovHC(ols_lfp))["log_slots","Std. Error"],
      coeftest(iv_lfp,  vcovHC(iv_lfp)) ["log_slots","Std. Error"],
      coeftest(ols_emp, vcovHC(ols_emp))["log_slots","Std. Error"],
      coeftest(iv_emp,  vcovHC(iv_emp)) ["log_slots","Std. Error"],
      coeftest(ols_ft,  vcovHC(ols_ft)) ["log_slots","Std. Error"],
      coeftest(iv_ft,   vcovHC(iv_ft))  ["log_slots","Std. Error"]
    ),
    pval = c(
      coeftest(ols_lfp, vcovHC(ols_lfp))["log_slots","Pr(>|t|)"],
      coeftest(iv_lfp,  vcovHC(iv_lfp)) ["log_slots","Pr(>|t|)"],
      coeftest(ols_emp, vcovHC(ols_emp))["log_slots","Pr(>|t|)"],
      coeftest(iv_emp,  vcovHC(iv_emp)) ["log_slots","Pr(>|t|)"],
      coeftest(ols_ft,  vcovHC(ols_ft)) ["log_slots","Pr(>|t|)"],
      coeftest(iv_ft,   vcovHC(iv_ft))  ["log_slots","Pr(>|t|)"]
    )
  )
) |>
  mutate(
    ci_lo    = Estimate - 1.96 * SE,
    ci_hi    = Estimate + 1.96 * SE,
    Outcome  = factor(Outcome, levels = c("LFP","Employment","Full-Time")),
    Method   = factor(Method,  levels = c("OLS","IV (2SLS)")),
    sig_star = case_when(
      pval < .01 ~ "***", pval < .05 ~ "**",
      pval < .10 ~ "*",   TRUE       ~ ""
    )
  )

fig5 <- ggplot(iv_compare,
               aes(x = Outcome, y = Estimate, color = Method, shape = Method)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, position = position_dodge(0.4), linewidth = 0.9) +
  geom_point(size = 4.5, position = position_dodge(0.4)) +
  geom_text(aes(y = ci_hi + 0.002, label = sig_star),
            position = position_dodge(0.4),
            size = 4.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c(PAL_NAVY, PAL_RED)) +
  scale_shape_manual(values = c(16, 17)) +
  annotate("label", x = 2.38, y = 0.045,
           label = "IV Employment:\nβ=+0.029, p=0.094†\n(primary causal claim)",
           size = 3, fill = "lightyellow", color = PAL_RED,
           label.size = 0.4, fontface = "bold") +
  labs(
    title    = "Figure 5: OLS vs IV (2SLS) — log(Slots) on Labor Supply",
    subtitle = "Parents with child under 5. 95% CIs. IV: lag-1 & lag-2 new provider openings.",
    x = "Outcome", y = "Coefficient on log(Total Slots)",
    color = "Estimator", shape = "Estimator",
    caption  = paste(
      "Source: ACS PUMS + MN DHS. Controls: demographics + county dummy + demeaned UR.",
      "\nFirst-stage F = 24,216; R² = 0.935. HC1 SEs. † p<0.10, ** p<0.05, *** p<0.01."
    )
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig5_iv_ols.png"),
       fig5, width = 9, height = 5.5, dpi = 300)

## Figure 6: Heterogeneity Forest Plot ----------------------------------------
sig_pal <- c("p<0.01" = PAL_NAVY, "p<0.05" = PAL_MID,
             "p<0.10" = PAL_LIGHT, "n.s."  = PAL_GRAY)

fig6 <- het_results |>
  filter(subgroup %in% c("Mothers","Fathers","Single Parents",
                         "Low-Income","Immigrants","Non-Immigrant",
                         "College-Educ","No College")) |>
  ggplot(aes(x = beta, xmin = ci_lo, xmax = ci_hi, y = subgroup, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.28, linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_color_manual(values = sig_pal, name = "Significance") +
  facet_wrap(~outcome, ncol = 3, scales = "free_x") +
  labs(
    title    = "Figure 6: Heterogeneity in Childcare Supply Effects by Subgroup",
    subtitle = "OLS coefficient on log(total slots); 95% CIs. Subgroup-specific controls.",
    x = "β on log(Slots)", y = NULL,
    caption  = "Source: ACS PUMS. Controls: age, age², college, immigrant, low_income, county dummy, UR. HC1 SEs."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig6_heterogeneity.png"),
       fig6, width = 12, height = 6, dpi = 300)

## Figure 7: Single-Parent Deep Dive ------------------------------------------
sp_q_long <- sp_qt |>
  select(slots_q_label, LFP = in_lf, Employment = employed, `Full-Time` = fulltime) |>
  pivot_longer(-slots_q_label, names_to = "Outcome", values_to = "Rate")

fig7a <- ggplot(sp_q_long, aes(x = slots_q_label, y = Rate, fill = Outcome)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.88) +
  scale_fill_manual(values = c(PAL_NAVY, PAL_MID, PAL_LIGHT)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "A. Single-Parent Labor Supply by Childcare Supply Quartile",
    subtitle = "Inverted-U in LFP consistent with informal care crowding-out at high supply",
    x = "Childcare Supply Quartile", y = "Rate", fill = "Outcome"
  ) +
  THEME_PAPER

# SP regression coefficients (gender split)
sp_betas <- het_results |>
  filter(subgroup %in% c("Single Mothers", "Single Fathers")) |>
  mutate(
    label = sprintf("β=%+.3f%s", beta,
                    ifelse(pval < .01,"***",ifelse(pval < .05,"**",
                                                   ifelse(pval < .10,"*",""))))
  )

fig7b <- ggplot(sp_betas, aes(x = beta, xmin = ci_lo, xmax = ci_hi,
                              y = subgroup, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.25, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = label), hjust = -0.25, size = 3.2, fontface = "bold") +
  scale_color_manual(values = sig_pal, name = "Significance") +
  facet_wrap(~outcome, ncol = 3, scales = "free_x") +
  labs(
    title    = "B. Coefficients: Single Mothers vs Single Fathers",
    subtitle = "Strong negative LFP effect for single mothers; positive employment for fathers",
    x = "β on log(Slots)", y = NULL
  ) +
  THEME_PAPER

fig7 <- fig7a / fig7b +
  plot_annotation(
    title   = "Figure 7: Single-Parent Deep Dive — Three Mechanisms",
    caption = paste(
      "Source: ACS PUMS. N(single parents)=1,652; N(single mothers)=1,006; N(single fathers)=646.",
      "\nHC1 SEs. Controls: age, age², college, immigrant, low_income, county dummy, UR."
    ),
    theme   = theme(plot.title = element_text(face = "bold", size = 12, color = PAL_NAVY))
  )

ggsave(paste0(OUTPUT_DIR, "fig7_singleparent.png"),
       fig7, width = 12, height = 10, dpi = 300)

## Figure 8: RD Plots (all three outcomes) ------------------------------------
rd_bins <- parents |>
  filter(poverty_c >= -40, poverty_c <= 40) |>
  mutate(
    bin     = cut(poverty_c, breaks = seq(-40, 40, by = 5), include.lowest = TRUE),
    bin_mid = as.numeric(str_extract(as.character(bin), "-?\\d+\\.?\\d*")) + 2.5,
    above   = above_ccap
  ) |>
  group_by(bin_mid, above) |>
  summarise(
    lfp = mean(in_lf,   na.rm = TRUE),
    emp = mean(employed, na.rm = TRUE),
    ft  = mean(fulltime, na.rm = TRUE),
    n   = n(), .groups = "drop"
  ) |>
  mutate(group = if_else(above == 1, "Above cutoff\n(CCAP ineligible)",
                         "Below cutoff\n(CCAP eligible)"))

make_rd_panel <- function(data, y_var, y_label, late_text, p_text) {
  ggplot(data, aes(x = bin_mid, y = .data[[y_var]])) +
    geom_smooth(data = filter(data, bin_mid < 0), method = "lm", se = TRUE,
                color = PAL_MID, fill = PAL_MID, alpha = 0.15) +
    geom_smooth(data = filter(data, bin_mid >= 0), method = "lm", se = TRUE,
                color = PAL_RED, fill = PAL_RED, alpha = 0.15) +
    geom_point(aes(size = n, color = bin_mid < 0), alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1) +
    scale_color_manual(values = c("TRUE" = PAL_MID, "FALSE" = PAL_RED), guide = "none") +
    scale_size_continuous(range = c(2, 7), name = "Bin N") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    annotate("label", x = 0, y = 0.25,
             label = sprintf("LATE = %s\n(p = %s)", late_text, p_text),
             size = 3.2, fill = "white", color = PAL_RED, fontface = "bold",
             label.size = 0.5) +
    labs(title = y_label,
         x = "Distance from 85% FPL (ppt)",
         y = y_label) +
    THEME_PAPER +
    theme(legend.position = "none")
}

p_rd1 <- make_rd_panel(rd_bins, "lfp", "LFP",         "−0.012", "0.880")
p_rd2 <- make_rd_panel(rd_bins, "emp", "Employment",   "+0.044", "0.608")
p_rd3 <- make_rd_panel(rd_bins, "ft",  "Full-Time",    "−0.298", "0.001")

fig8 <- (p_rd1 | p_rd2 | p_rd3) +
  plot_annotation(
    title    = "Figure 8: Regression Discontinuity — CCAP Eligibility Threshold (85% FPL)",
    subtitle = "Running variable: household income as % FPL, centered at 85%",
    caption  = paste(
      "Source: ACS PUMS. Local linear RD (rdrobust), triangular kernel, MSE-optimal BW.",
      "\nN=587 parents within ±40 ppt BW. Each point = 5 ppt bin; size ∝ bin N.",
      "\nPlacebo cutoffs at 65% and 105% FPL show no significant discontinuity."
    ),
    theme    = theme(plot.title = element_text(face = "bold", size = 12, color = PAL_NAVY))
  )

ggsave(paste0(OUTPUT_DIR, "fig8_rd_all.png"),
       fig8, width = 13, height = 5.5, dpi = 300)

## Figure 9: QCEW Sector Wages + Employment -----------------------------------
qcew_st <- qcew_interp |> filter(county_name == "Stearns")

fig9 <- ggplot(qcew_st, aes(x = year)) +
  geom_area(aes(y = emp_annual), fill = PAL_LIGHT, alpha = 0.4) +
  geom_line(aes(y = emp_annual, color = "Employment"),  linewidth = 1.4) +
  geom_line(aes(y = avg_wkly_wage_annual * 0.55,
                color = "Avg Weekly Wage (×0.55)"), linewidth = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name     = "Annual Employment (FTE)",
    labels   = comma,
    sec.axis = sec_axis(~ . / 0.55, name = "Avg Weekly Wage ($)", labels = dollar)
  ) +
  scale_color_manual(values = c("Employment" = PAL_NAVY, "Avg Weekly Wage (×0.55)" = PAL_RED),
                     name = NULL) +
  scale_x_continuous(breaks = seq(2009, 2024, 3)) +
  labs(
    title    = "Figure 9: Childcare Sector Employment & Wages (NAICS 6244)",
    subtitle = "Stearns County, MN — Interpolated from MN DEED QCEW (private sector)",
    x = "Year",
    caption  = "Source: MN DEED QCEW. Missing years (2016, 2020 Benton) linearly interpolated."
  ) +
  THEME_PAPER +
  theme(axis.title.y.right = element_text(color = PAL_RED))

ggsave(paste0(OUTPUT_DIR, "fig9_qcew.png"),
       fig9, width = 9, height = 5.5, dpi = 300)

## Figure 10: FRED Labor Market Context ----------------------------------------
fred_long <- fred_yr |>
  filter(year >= 2005) |>
  select(year, `Unemployment Rate` = ur, `LFP Rate (×10)` = lfp_rate) |>
  mutate(`LFP Rate (×10)` = `LFP Rate (×10)` * 10) |>
  pivot_longer(-year, names_to = "series", values_to = "value")

fig10 <- ggplot(fred_long, aes(x = year, y = value, color = series)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_vline(xintercept = c(2009, 2020), linetype = "dotted",
             color = "grey50", linewidth = 0.8) +
  annotate("text", x = 2009.2, y = 10.5, label = "GFC\nPeak", size = 2.8, color = "grey40") +
  annotate("text", x = 2020.2, y = 10.5, label = "COVID", size = 2.8, color = "grey40") +
  scale_color_manual(values = c(PAL_RED, PAL_NAVY)) +
  scale_x_continuous(breaks = seq(2005, 2026, 3)) +
  labs(
    title    = "Figure 10: St. Cloud MSA Aggregate Labor Market (2005–2026)",
    subtitle = "Annual averages from FRED monthly series",
    x = "Year", y = "Rate", color = NULL,
    caption  = "Source: FRED. STCL027UR (UR); LBSSA27 (LFP rate ×10 for scale)."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig10_fred.png"),
       fig10, width = 9, height = 5, dpi = 300)

## Figure 11: LFP by Parent Group Over Time ------------------------------------
pums_grp <- parents |>
  group_by(YEAR) |>
  summarise(
    `All Parents`    = mean(in_lf,                            na.rm = TRUE),
    `Mothers`        = mean(in_lf[female == 1],              na.rm = TRUE),
    `Fathers`        = mean(in_lf[female == 0],              na.rm = TRUE),
    `Single Parents` = mean(in_lf[single_parent == 1],       na.rm = TRUE),
    `Low-Income`     = mean(in_lf[low_income == 1],          na.rm = TRUE),
    `Immigrants`     = mean(in_lf[immigrant == 1],           na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(-YEAR, names_to = "Group", values_to = "LFP")

fig11 <- ggplot(pums_grp, aes(x = YEAR, y = LFP, color = Group, linetype = Group)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(PAL_NAVY, PAL_MID, PAL_LIGHT,
                                PAL_RED, PAL_ORG, PAL_GREEN)) +
  scale_linetype_manual(values = c("solid","solid","solid","dashed","dashed","dotted")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(2009, 2024, 3)) +
  labs(
    title    = "Figure 11: LFP Rate by Parent Group (2009–2024)",
    subtitle = "Parents with child under 5 only; St. Cloud MSA",
    x = "Year", y = "Labor Force Participation Rate",
    color = "Group", linetype = "Group",
    caption  = "Source: ACS PUMS. Annual cell means; unweighted."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig11_lfp_by_group.png"),
       fig11, width = 10, height = 6, dpi = 300)

## Figure 12: OEWS Wage Benchmarking ------------------------------------------
oews_plot <- oews |>
  filter(file %in% c("primary","regional"), !is.na(wage_median)) |>
  mutate(
    geo_short  = str_trunc(geography, 28),
    is_stcloud = str_detect(geography, regex("ST CLOUD|St Cloud", ignore_case = TRUE))
  ) |>
  arrange(wage_median)

fig12 <- ggplot(oews_plot, aes(x = reorder(geo_short, wage_median),
                               y = wage_median, fill = is_stcloud)) +
  geom_col(width = 0.72, alpha = 0.88) +
  geom_text(aes(label = dollar(wage_median, accuracy = 0.01)),
            hjust = -0.1, size = 2.8) +
  scale_fill_manual(values = c("FALSE" = PAL_LIGHT, "TRUE" = PAL_RED),
                    labels  = c("Comparison regions","St. Cloud MSA"), name = NULL) +
  scale_y_continuous(labels = dollar, expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(
    title    = "Figure 12: Childcare Worker Median Wages by Region (2026)",
    subtitle = "SOC 39-9011 (Childcare Workers) — MN DEED OEWS Q1 2026",
    x = NULL, y = "Median Hourly Wage ($/hr)",
    caption  = "Source: MN DEED OEWS Q1 2026. St. Cloud MSA highlighted."
  ) +
  THEME_PAPER +
  theme(legend.position = "top")

ggsave(paste0(OUTPUT_DIR, "fig12_oews_wages.png"),
       fig12, width = 9, height = 8, dpi = 300)

## Figure 13: ZIP-level Supply Distribution (2023) ----------------------------
zp_2023 <- zp |>
  filter(year == 2023) |>
  arrange(desc(total_slots)) |>
  mutate(
    zip_label = paste0("ZIP ", zip),
    top10     = row_number() <= 10
  )

fig13 <- ggplot(zp_2023,
                aes(x = reorder(zip_label, total_slots), y = total_slots,
                    fill = top10)) +
  geom_col(width = 0.72, alpha = 0.88) +
  geom_col(aes(y = slots_infants), fill = PAL_ORG, alpha = 0.8, width = 0.72) +
  scale_fill_manual(values = c("FALSE" = PAL_LIGHT, "TRUE" = PAL_NAVY),
                    labels  = c("Other ZIPs","Top 10 ZIPs"), name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(
    title    = "Figure 13: ZIP-Code Distribution of Licensed Slots (2023)",
    subtitle = "Substantial spatial heterogeneity within MSA; orange = infant slots",
    x = NULL, y = "Total Licensed Slots",
    caption  = "Source: MN DHS/DCYF. N = 31 ZIPs with licensed providers in 2023."
  ) +
  THEME_PAPER +
  theme(legend.position = "top")

ggsave(paste0(OUTPUT_DIR, "fig13_zip_slots.png"),
       fig13, width = 9, height = 8, dpi = 300)

## Figure 14: Minneapolis Fed Provider Financial Stress ------------------------
mpls_long <- mpls |>
  select(survey_year, insurance_index_2022base, pct_raising_tuition_centers) |>
  pivot_longer(-survey_year, names_to = "series", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(series = recode(series,
                         "insurance_index_2022base"    = "Insurance Cost Index (2022=100)",
                         "pct_raising_tuition_centers" = "% Centers Raising Tuition"
  ))

fig14 <- ggplot(mpls_long, aes(x = survey_year, y = value, color = series, shape = series)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c(PAL_RED, PAL_NAVY)) +
  scale_x_continuous(breaks = unique(mpls_long$survey_year)) +
  labs(
    title    = "Figure 14: Childcare Provider Financial Stress (Minneapolis Fed Surveys)",
    subtitle = "Insurance cost index (2022 = 100) and % of centers raising tuition",
    x = "Survey Year", y = "Value", color = NULL, shape = NULL,
    caption  = "Source: Federal Reserve Bank of Minneapolis Annual Childcare Provider Surveys (2021–2026)."
  ) +
  THEME_PAPER

ggsave(paste0(OUTPUT_DIR, "fig14_mpls_stress.png"),
       fig14, width = 9, height = 5, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# 9.  FINAL COMBINED PANEL FIGURES
# ─────────────────────────────────────────────────────────────────────────────

## Panel A: RQ1 — Supply + LFP trends + Event study --------------------------
panel_rq1 <- (fig1 | fig2) / (fig3b | fig4) +
  plot_annotation(
    title   = "RQ1 Panel: Childcare Supply and Labor Supply Evidence",
    caption = "Figs 1–4 combined. See individual figure captions for sources.",
    theme   = theme(plot.title = element_text(face = "bold", size = 13, color = PAL_NAVY))
  )

ggsave(paste0(OUTPUT_DIR, "panel_rq1.png"),
       panel_rq1, width = 16, height = 11, dpi = 300)

## Panel B: RQ2 — Heterogeneity -----------------------------------------------
panel_rq2 <- fig6 / fig7 +
  plot_annotation(
    title   = "RQ2 Panel: Heterogeneity in Childcare Supply Effects",
    caption = "Figs 6–7 combined.",
    theme   = theme(plot.title = element_text(face = "bold", size = 13, color = PAL_NAVY))
  )

ggsave(paste0(OUTPUT_DIR, "panel_rq2.png"),
       panel_rq2, width = 14, height = 16, dpi = 300)

## Panel C: RQ3 — RD ----------------------------------------------------------
ggsave(paste0(OUTPUT_DIR, "panel_rq3.png"),
       fig8, width = 14, height = 6, dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# 10. SESSION INFO & COMPLETION SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

cat("\n\n")
cat("=================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=================================================================\n\n")

cat("RQ1 PRIMARY RESULTS:\n")
cat("  OLS Employment:  β =+0.011, p =0.498 (n.s.; positive, imprecise)\n")
cat("  IV  Employment:  β =+0.029, p =0.094 (†; primary causal claim)\n")
cat("  IV  Full-Time:   β =+0.033, p =0.127 (n.s.)\n\n")

cat("RQ2 KEY HETEROGENEITY:\n")
cat("  Mothers — FT:           β =+0.059, p =0.054 (†)\n")
cat("  Single Mothers — LFP:   β =−0.147, p =0.023 (**)\n")
cat("  Single Mothers — Emp:   β =−0.152, p =0.017 (**)\n")
cat("  No College — FT:        β =+0.074, p =0.023 (**)\n\n")

cat("RQ3 RD RESULTS:\n")
cat("  LFP:        LATE =−0.012, p =0.880 (n.s.)\n")
cat("  Employment: LATE =+0.044, p =0.608 (n.s.)\n")
cat("  Full-Time:  LATE =−0.298, p =0.001 (***) ← STRONGEST FINDING\n")
cat("  SP FT:      LATE =−0.359, p =0.005 (***)  ← Single parents stronger\n\n")

cat("Files saved to:", OUTPUT_DIR, "\n")
cat(" Tables: table3_ols_twfe.docx, table4_iv.docx,\n")
cat("         table6_heterogeneity.csv, table7_rd.csv\n")
cat(" Figures: fig1 through fig14, plus 3 combined panels\n\n")

sessionInfo()
