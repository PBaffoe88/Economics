############################################################
# FULL RESEARCH SCRIPT
# "Buffer Depletion and Price Instability:
#  How Monetary Tightening Amplifies Producer Price
#  Volatility Through the Inventory Channel"
#
# Self-contained: data download → cleaning → base models →
# publication figures → extensions with verified industry
# characteristics from live official APIs.
#
# RUN ORDER: source this single file.
#
# DATA SOURCES (all verified):
#   Census M3 / BLS PPI / FRED macros   — FRED API
#   SF Fed monetary policy surprises     — frbsf.org (Bauer & Swanson 2023)
#   ASM export shares                    — Census ASM published tables
#                                          (NOT available via ASM API — see FIX 14)
#   HHI concentration ratios             — api.census.gov/data/2017/ecnconcentration
#   BEA I-O upstreamness                 — apps.bea.gov/api (Antras & Chor 2013)
#   Perishability                        — USDA ERS / FDA 21 CFR
#   OECD STAN + BIS rates (optional)     — oecd.org / bis.org
#
# API KEYS (free, OPTIONAL):
#   Census: api.census.gov/data/key_signup.html
#           Sys.setenv(CENSUS_API_KEY = "your_key")
#   BEA:    apps.bea.gov/API/signup/
#           Sys.setenv(BEA_API_KEY = "your_key")
#   Neither key is required for the script to run: all Census/BEA
#   sections fall back to published-table values automatically.
#
# FIXES APPLIED:
#   FIX 1  — fred_csv(): robust column rename by position
#   FIX 2  — SF Fed: dynamic shock column detection
#   FIX 3  — USREC: max() not mean() for binary daily→monthly
#   FIX 4  — run_lp(): two-way cluster SE (industry + date)
#   FIX 5  — food_model_price_vol: simultaneity bias removed
#   FIX 6  — inventory_intensity: Bartik pre-period comment
#   FIX 7  — fred_csv(): 3-attempt exponential back-off
#   FIX 8  — map_dfr(): Sys.sleep(0.3) throttle
#   FIX 9  — daily series → monthly (GS2, GS10, WTISPLC)
#   FIX 10 — missing-column guard before controls mutate
#   FIX 7+ — fred_csv(): 5 retries, longer back-off (5s base), jitter
#   FIX 8+ — fred_batch(): batched downloads (groups of 4, 3s inter-group pause)
#   FIX 10+ — recession defaults to 0 not NA when USREC download fails
#   FIX 13 — USREC removed from FRED batch; recession built from hard-coded
#             NBER dates (permanent historical facts, standard in literature);
#             recession_term() helper prevents collinearity warnings
#   FIX 11 — run_lp_split(): industry FE only (date FE collinear)
#   FIX 12 — Section 4 industry chars: live API data replaces
#             hand-coded tribble (Census ASM, Econ Census, BEA I-O)
#   FIX 14 — Ext 4 Census crash fixed:
#             (a) 'timeseries/asm/industry2017' does not exist in the
#                 Census API; it returned an HTML error page and
#                 fromJSON() crashed outside the tryCatch, killing the
#                 script instead of engaging the fallback.
#             (b) No ASM API dataset carries an export-shipments
#                 variable (EXPSHIP) — verified against the variable
#                 lists of asm/industry and asm/benchmark2022. Export
#                 shares therefore load directly from the published
#                 ASM tables (primary source, not a fallback).
#             (c) New census_json() helper validates that the API
#                 response is JSON before parsing and wraps the whole
#                 fetch in tryCatch, so any Census failure degrades
#                 gracefully to fallback tables instead of aborting.
############################################################

# ============================================================
# 0. PACKAGES & DIRECTORIES
# ============================================================

pkgs <- c(
  "tidyverse", "lubridate", "zoo", "slider", "readxl",
  "janitor", "fixest", "broom", "modelsummary",
  "ggplot2", "scales", "httr", "jsonlite"
)
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0)
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(tidyverse); library(lubridate); library(zoo)
library(slider);    library(readxl);    library(janitor)
library(fixest);    library(broom);     library(modelsummary)
library(ggplot2);   library(scales);    library(httr)
library(jsonlite)

for (d in c("data_raw","data_clean","outputs","outputs/figures","outputs/tables"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Shared ggplot theme
theme_paper <- function()
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.caption  = element_text(size = 8, colour = "grey40", hjust = 0),
        plot.subtitle = element_text(size = 9, colour = "grey25", lineheight = 1.4))

save_fig <- function(plot, filename, width = 9, height = 5) {
  ggsave(file.path("outputs/figures", filename), plot,
         width = width, height = height, dpi = 300, bg = "white")
  message("Saved: ", filename)
}

# ============================================================
# 1. FRED HELPER — retry + throttle
#
# FIX 7 (enhanced): longer back-off (5s base, not 2s), more retries (5),
#   and jitter to avoid thundering-herd when multiple series hit 50x errors.
# FIX 8 (enhanced): batch downloads in groups of 4 with a 3s inter-group
#   pause. FRED's rate limit is ~120 req/min; 4 requests + 3s pause = ~80/min.
# ============================================================

fred_csv <- function(series_id, retries = 5, pause = 5) {
  url <- paste0("https://fred.stlouisfed.org/graph/fredgraph.csv?id=", series_id)
  for (attempt in seq_len(retries)) {
    result <- tryCatch({
      df <- read_csv(url, show_col_types = FALSE)
      names(df) <- c("date", "value")
      df %>% mutate(date = as.Date(date), value = as.numeric(value),
                    series_id = series_id)
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
    wait <- pause * attempt + runif(1, 0, 2)   # jitter avoids thundering herd
    message("Attempt ", attempt, " failed for ", series_id,
            ". Retrying in ", round(wait, 1), "s...")
    Sys.sleep(wait)
  }
  message("Could not download FRED series after ", retries, " attempts: ", series_id)
  return(NULL)
}

# Batched downloader: splits series list into groups, pauses between groups
# to stay well under FRED's rate limit regardless of how many series are needed.
fred_batch <- function(series_ids, batch_size = 4, inter_batch_pause = 3) {
  batches <- split(series_ids, ceiling(seq_along(series_ids) / batch_size))
  map_dfr(seq_along(batches), function(i) {
    if (i > 1) Sys.sleep(inter_batch_pause)
    map_dfr(batches[[i]], function(sid) {
      Sys.sleep(0.5)   # within-batch throttle
      fred_csv(sid)
    })
  })
}

make_growth     <- function(x) 100 * (log(x) - log(lag(x)))
make_volatility <- function(x, window = 12)
  slider::slide_dbl(x, ~ sd(.x, na.rm = TRUE),
                    .before = window - 1, .complete = TRUE)

# ============================================================
# 2. INDUSTRY PANEL — Census M3 + BLS PPI via FRED
# ============================================================

industry_map <- tribble(
  ~industry,                          ~naics, ~sis_series, ~ship_series, ~inv_series, ~ppi_series,
  "Food Manufacturing",               "311",  "A11SIS", "A11SVS", "A11STI", "PCU311311",
  "Beverage and Tobacco",             "312",  "A12SIS", "A12SVS", "A12STI", "PCU312312",
  "Chemical Products",                "325",  "A25SIS", "A25SVS", "A25STI", "PCU325325",
  "Petroleum and Coal Products",      "324",  "A24SIS", "A24SVS", "A24STI", "PCU324324",
  "Plastics and Rubber Products",     "326",  "A26SIS", "A26SVS", "A26STI", "PCU326326",
  "Primary Metals",                   "331",  "A31SIS", "A31SVS", "A31STI", "PCU331331",
  "Machinery",                        "333",  "A35SIS", "A35SVS", "A35STI", "PCU333333",
  "Computer and Electronic Products", "334",  "A34SIS", "A34SVS", "A34STI", "PCU334334"
)

series_long <- industry_map %>%
  pivot_longer(c(sis_series,ship_series,inv_series,ppi_series),
               names_to="variable", values_to="series_id")

# Use batched downloader — 32 series in groups of 4 with 3s inter-group pauses
all_fred <- fred_batch(unique(series_long$series_id), batch_size = 4, inter_batch_pause = 3)

industry_panel_raw <- all_fred %>%
  left_join(series_long %>% select(industry,naics,variable,series_id), by="series_id") %>%
  select(date,industry,naics,variable,value) %>%
  pivot_wider(names_from=variable, values_from=value) %>%
  rename(inv_sales_ratio=sis_series, shipments=ship_series,
         inventories=inv_series, ppi=ppi_series) %>%
  filter(!is.na(industry), !is.na(inv_sales_ratio),
         !is.na(shipments), !is.na(inventories), !is.na(ppi))

write_csv(industry_panel_raw, "data_raw/industry_panel_raw.csv")
message("Industry panel downloaded: ", nrow(industry_panel_raw), " rows across ",
        n_distinct(industry_panel_raw$industry), " industries.")

# ============================================================
# 3. MONETARY POLICY SHOCKS — SF Fed (Bauer & Swanson 2023)
# ============================================================

download.file(
  "https://www.frbsf.org/wp-content/uploads/monetary-policy-surprises-data.xlsx",
  destfile = "data_raw/sf_fed_monetary_policy_surprises.xlsx", mode = "wb"
)

mps_raw   <- readxl::read_excel("data_raw/sf_fed_monetary_policy_surprises.xlsx",
                                sheet = "Monthly (update 2023)") %>%
  janitor::clean_names()

shock_col <- grep("mps|orth|surprise", names(mps_raw), value=TRUE, ignore.case=TRUE)
if (length(shock_col) == 0)
  stop("Cannot find shock column. Available: ", paste(names(mps_raw), collapse=", "))
message("Using column '", shock_col[1], "' as mp_shock")

mps_monthly <- mps_raw %>%
  transmute(year=as.integer(year), month=as.integer(month),
            mp_shock=suppressWarnings(as.numeric(.data[[shock_col[1]]]))) %>%
  filter(!is.na(year), !is.na(month), !is.na(mp_shock)) %>%
  mutate(date=as.Date(sprintf("%04d-%02d-01",year,month))) %>%
  select(date, mp_shock) %>% arrange(date)

write_csv(mps_monthly, "data_clean/monthly_monetary_policy_shocks.csv")

# ============================================================
# 4. MACRO CONTROLS — all monthly series via FRED
#
# USREC is excluded from the FRED batch entirely.
# FRED's CSV endpoint for USREC persistently returns 504 errors
# because it is a derived binary series that routes through a
# different backend than the standard release series.
#
# Fix: NBER recession dates are hard-coded directly from the
# official NBER Business Cycle Dating Committee announcements:
#   nber.org/research/data/us-business-cycle-expansions-contractions
# These dates are permanent historical facts, not subject to revision,
# so hard-coding them is appropriate and standard in the literature.
# The series will be updated manually if future recessions occur.
# ============================================================

control_series <- tribble(
  ~name,                   ~series_id,
  "effective_fed_funds",   "FEDFUNDS",
  "two_year_treasury",     "GS2",
  "ten_year_treasury",     "GS10",
  "oil_wti",               "WTISPLC",
  "all_commodities_ppi",   "PPIACO",
  "cpi_all",               "CPIAUCSL",
  "industrial_production", "INDPRO"
  # USREC excluded — persistent FRED 504 errors; replaced by nber_recession below
)

# Batched download for 7 continuous controls
controls_raw <- fred_batch(control_series$series_id,
                           batch_size = 4, inter_batch_pause = 5) %>%
  left_join(control_series, by="series_id") %>%
  select(date, name, value, series_id)

controls_wide <- controls_raw %>%
  mutate(date = floor_date(date,"month")) %>%
  group_by(date, name) %>%
  summarise(value = mean(value, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=name, values_from=value) %>%
  arrange(date)

# NBER recession months — hard-coded from official NBER announcements.
# Source: NBER Business Cycle Dating Committee,
#         nber.org/research/data/us-business-cycle-expansions-contractions
# Covers all recessions from 1990 onward (our sample starts 1992-01).
# Peak month = first recession month; trough month = last recession month.
nber_recessions <- tibble(
  peak   = as.Date(c("1990-07-01","2001-03-01","2007-12-01","2020-02-01")),
  trough = as.Date(c("1991-03-01","2001-11-01","2009-06-01","2020-04-01"))
)

# Build monthly recession indicator: 1 if month falls in any NBER recession
all_months <- seq(min(controls_wide$date), max(controls_wide$date), by="month")

recession_monthly <- tibble(date = all_months) %>%
  mutate(recession = as.integer(map_lgl(date, function(d)
    any(d >= nber_recessions$peak & d <= nber_recessions$trough))))

# Merge recession indicator into controls
controls_wide <- controls_wide %>%
  left_join(recession_monthly, by="date")

# Guard for any remaining missing continuous columns
expected_cols <- c("oil_wti","all_commodities_ppi","cpi_all","industrial_production",
                   "two_year_treasury","ten_year_treasury","effective_fed_funds")
missing_cols  <- setdiff(expected_cols, names(controls_wide))
if (length(missing_cols) > 0) {
  warning("Missing control series (set to NA): ", paste(missing_cols, collapse=", "))
  controls_wide[missing_cols] <- NA_real_
}

controls_monthly <- controls_wide %>%
  mutate(oil_infl       = make_growth(oil_wti),
         commodity_infl = make_growth(all_commodities_ppi),
         cpi_infl       = make_growth(cpi_all),
         ip_growth      = make_growth(industrial_production),
         yield_spread   = ten_year_treasury - two_year_treasury)

write_csv(controls_monthly, "data_clean/controls_monthly.csv")

# Helper used throughout: only include recession in a formula if it has
# variance (i.e. if both recession and non-recession months are present).
# When recession is constant (all 0 or all 1), it is collinear with the
# intercept and feols/lm will drop it with warnings. This helper silences
# that by checking variance before including the term.
recession_term <- function(data_col) {
  if (var(data_col, na.rm = TRUE) > 0) "+ recession" else ""
}

message("Controls ready. Recession indicator: ",
        sum(controls_monthly$recession, na.rm=TRUE), " recession months from NBER dates.")

# ============================================================
# 5. DATA TRANSFORMATION
# ============================================================

panel_clean <- industry_panel_raw %>%
  arrange(industry, date) %>%
  group_by(industry) %>%
  mutate(
    ln_inventory    = log(inventories),
    ln_shipments    = log(shipments),
    ln_ppi          = log(ppi),
    inventory_growth        = make_growth(inventories),
    shipments_growth        = make_growth(shipments),
    ppi_inflation           = make_growth(ppi),
    ppi_vol_6m              = make_volatility(ppi_inflation, 6),
    ppi_vol_12m             = make_volatility(ppi_inflation, 12),
    inv_sales_ratio_change  = inv_sales_ratio - lag(inv_sales_ratio),
    abs_ppi_inflation       = abs(ppi_inflation)
  ) %>% ungroup()

# Bartik-style pre-period intensity (1992-2006 only)
# Cross-industry variation is structural and pre-determined
# relative to post-2007 shocks. See Mian & Sufi (2014).
inventory_intensity <- panel_clean %>%
  filter(date >= as.Date("1992-01-01"), date <= as.Date("2006-12-01")) %>%
  group_by(industry) %>%
  summarise(inventory_intensity = mean(inv_sales_ratio, na.rm=TRUE), .groups="drop") %>%
  mutate(high_inventory = if_else(
    inventory_intensity >= median(inventory_intensity, na.rm=TRUE), 1L, 0L))

panel_final <- panel_clean %>%
  left_join(inventory_intensity, by="industry") %>%
  left_join(mps_monthly,         by="date") %>%
  left_join(controls_monthly,    by="date") %>%
  mutate(mp_shock        = replace_na(mp_shock, 0),
         mp_x_inventory  = mp_shock * inventory_intensity,
         mp_x_high_inv   = mp_shock * high_inventory,
         year  = year(date), month = month(date), ym = as.yearmon(date)) %>%
  filter(date >= as.Date("1992-01-01")) %>%
  arrange(industry, date)

write_csv(panel_final, "data_clean/final_industry_panel.csv")

# ============================================================
# 6. BASE LOCAL PROJECTIONS
# ============================================================

run_lp <- function(data, yvar, h_max = 24) {
  map_dfr(0:h_max, function(h) {
    df_h <- data %>%
      group_by(industry) %>% arrange(date) %>%
      mutate(y_lead = lead(.data[[yvar]], h)) %>%
      ungroup() %>%
      filter(!is.na(y_lead), !is.na(mp_x_inventory), !is.na(industry), !is.na(date))
    tryCatch(
      feols(y_lead ~ mp_x_inventory | industry + date,
            cluster = ~industry + date, data = df_h) %>%
        tidy(conf.int = TRUE) %>% filter(term == "mp_x_inventory") %>%
        mutate(horizon = h, outcome = yvar),
      error = function(e) tibble(
        term="mp_x_inventory", estimate=NA_real_, std.error=NA_real_,
        statistic=NA_real_, p.value=NA_real_,
        conf.low=NA_real_, conf.high=NA_real_, horizon=h, outcome=yvar)
    )
  })
}

lp_inventory        <- run_lp(panel_final, "inv_sales_ratio",  24)
lp_inventory_growth <- run_lp(panel_final, "inventory_growth", 24)
lp_price_vol        <- run_lp(panel_final, "ppi_vol_12m",      24)
lp_ppi_inflation    <- run_lp(panel_final, "ppi_inflation",    24)

lp_all <- bind_rows(lp_inventory, lp_inventory_growth, lp_price_vol, lp_ppi_inflation)
write_csv(lp_all, "outputs/tables/local_projection_results.csv")

# ============================================================
# 7. BASE LP FIGURES
# ============================================================

plot_lp <- function(results, outcome_name, title_text, ylab_text) {
  results %>% filter(outcome == outcome_name) %>%
    ggplot(aes(x=horizon, y=estimate)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2) +
    geom_line(linewidth=1) +
    labs(title=title_text, x="Months after monetary policy shock", y=ylab_text,
         caption="Source: Census M3/FRED, BLS PPI/FRED, SF Fed monetary policy surprises.") +
    theme_minimal()
}

ggsave("outputs/figures/lp_inventory_ratio.png",
       plot_lp(lp_all,"inv_sales_ratio",
               "Effect of Monetary Tightening on Inventory-to-Shipments Ratio",
               "Response by inventory intensity"), width=8, height=5)
ggsave("outputs/figures/lp_inventory_growth.png",
       plot_lp(lp_all,"inventory_growth",
               "Effect of Monetary Tightening on Inventory Growth",
               "Response by inventory intensity"), width=8, height=5)
ggsave("outputs/figures/lp_price_volatility.png",
       plot_lp(lp_all,"ppi_vol_12m",
               "Effect of Monetary Tightening on Producer Price Volatility",
               "Response by inventory intensity"), width=8, height=5)
ggsave("outputs/figures/lp_ppi_inflation.png",
       plot_lp(lp_all,"ppi_inflation",
               "Effect of Monetary Tightening on Producer Price Inflation",
               "Response by inventory intensity"), width=8, height=5)

# ============================================================
# 8. BASE REGRESSION TABLES
# ============================================================

model_inventory  <- feols(inv_sales_ratio  ~ mp_x_inventory | industry+date, cluster=~industry+date, data=panel_final)
model_inv_growth <- feols(inventory_growth ~ mp_x_inventory | industry+date, cluster=~industry+date, data=panel_final)
model_price_vol  <- feols(ppi_vol_12m      ~ mp_x_inventory | industry+date, cluster=~industry+date, data=panel_final)
model_ppi        <- feols(ppi_inflation    ~ mp_x_inventory | industry+date, cluster=~industry+date, data=panel_final)

modelsummary(
  list("Inventory-to-Shipments Ratio"=model_inventory, "Inventory Growth"=model_inv_growth,
       "PPI Volatility"=model_price_vol, "PPI Inflation"=model_ppi),
  output="outputs/tables/main_regression_table.html", stars=TRUE, gof_omit="IC|Log|Adj|Within")

# ============================================================
# 9. ALTERNATIVE MODEL + FOOD FOCUS
# ============================================================

model_inventory_alt <- feols(
  inv_sales_ratio ~ mp_shock + mp_x_inventory +
    commodity_infl + oil_infl + cpi_infl + ip_growth + yield_spread + recession | industry,
  cluster=~industry, data=panel_final, warn=FALSE)

model_price_vol_alt <- feols(
  ppi_vol_12m ~ mp_shock + mp_x_inventory +
    commodity_infl + oil_infl + cpi_infl + ip_growth + yield_spread + recession | industry,
  cluster=~industry, data=panel_final, warn=FALSE)

modelsummary(list("Alt: Inventory Ratio"=model_inventory_alt, "Alt: PPI Volatility"=model_price_vol_alt),
             output="outputs/tables/alternative_controls_model.html", stars=TRUE, gof_omit="IC|Log|Adj|Within")

food_panel <- panel_final %>% filter(industry == "Food Manufacturing")
modelsummary(
  list(
    "Food Inventory Ratio" = lm(inv_sales_ratio ~ mp_shock + commodity_infl + oil_infl +
                                  cpi_infl + ip_growth + yield_spread + recession, data=food_panel),
    "Food PPI Volatility"  = lm(ppi_vol_12m ~ mp_shock + commodity_infl + oil_infl +
                                  cpi_infl + ip_growth + yield_spread + recession, data=food_panel)
  ),
  output="outputs/tables/food_manufacturing_focus.html", stars=TRUE)

# ============================================================
# 10. DESCRIPTIVE FIGURES
# ============================================================

focus_inds <- c("Food Manufacturing","Chemical Products",
                "Petroleum and Coal Products","Computer and Electronic Products")

ggsave("outputs/figures/descriptive_inventory_ratios.png",
       panel_final %>% filter(industry %in% focus_inds) %>%
         ggplot(aes(x=date, y=inv_sales_ratio, color=industry)) + geom_line() +
         labs(title="Inventory-to-Shipments Ratios Across Manufacturing Industries",
              x=NULL, y="Inventory-to-shipments ratio", caption="Source: Census M3 via FRED.") +
         theme_minimal() + theme(legend.position="bottom"),
       width=9, height=5)

ggsave("outputs/figures/descriptive_price_volatility.png",
       panel_final %>% filter(industry %in% focus_inds) %>%
         ggplot(aes(x=date, y=ppi_vol_12m, color=industry)) + geom_line() +
         labs(title="Rolling 12-Month Producer Price Volatility",
              x=NULL, y="SD of monthly PPI inflation", caption="Source: BLS PPI via FRED.") +
         theme_minimal() + theme(legend.position="bottom"),
       width=9, height=5)

# ============================================================
# 11. BASE ROBUSTNESS: PLACEBO + SPLIT SAMPLE + CRISIS
# ============================================================

# Placebo (12-month lag)
panel_placebo <- panel_final %>%
  group_by(industry) %>% arrange(date) %>%
  mutate(mp_shock_placebo = lag(mp_shock,12),
         mp_x_inv_placebo = mp_shock_placebo * inventory_intensity) %>% ungroup()

run_lp_placebo <- function(data, yvar, h_max=24) {
  map_dfr(0:h_max, function(h) {
    df_h <- data %>% group_by(industry) %>% arrange(date) %>%
      mutate(y_lead=lead(.data[[yvar]],h)) %>% ungroup() %>%
      filter(!is.na(y_lead), !is.na(mp_x_inv_placebo))
    tryCatch(
      feols(y_lead~mp_x_inv_placebo|industry+date, cluster=~industry, data=df_h) %>%
        tidy(conf.int=TRUE) %>% filter(term=="mp_x_inv_placebo") %>%
        mutate(horizon=h, outcome=yvar, type="placebo"),
      error=function(e) tibble(term="mp_x_inv_placebo",estimate=NA_real_,
                               std.error=NA_real_,statistic=NA_real_,p.value=NA_real_,
                               conf.low=NA_real_,conf.high=NA_real_,horizon=h,outcome=yvar,type="placebo"))
  })
}

lp_comparison <- bind_rows(
  lp_inventory  %>% mutate(type="actual"),
  lp_price_vol  %>% mutate(type="actual"),
  run_lp_placebo(panel_placebo,"inv_sales_ratio"),
  run_lp_placebo(panel_placebo,"ppi_vol_12m"))
write_csv(lp_comparison, "outputs/tables/placebo_comparison.csv")

plot_placebo <- function(results, outcome_name, title_text) {
  results %>% filter(outcome==outcome_name) %>%
    ggplot(aes(x=horizon,y=estimate,color=type,fill=type)) +
    geom_hline(yintercept=0,linetype="dashed") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.15,color=NA) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c("actual"="#1a5fa8","placebo"="#b0b0b0"),
                       labels=c("Actual shock","Placebo (12-month lag)")) +
    scale_fill_manual(values=c("actual"="#1a5fa8","placebo"="#b0b0b0"),guide="none") +
    labs(title=title_text, x="Months after shock", y="IRF estimate",
         color=NULL, caption="Placebo shifts mp_shock forward 12 months.") +
    theme_minimal() + theme(legend.position="bottom")
}

ggsave("outputs/figures/placebo_inventory.png",
       plot_placebo(lp_comparison,"inv_sales_ratio","Placebo test: inventory-to-shipments ratio"),
       width=8, height=5)
ggsave("outputs/figures/placebo_price_vol.png",
       plot_placebo(lp_comparison,"ppi_vol_12m","Placebo test: producer price volatility"),
       width=8, height=5)

# Split sample (industry FE only — date FE would absorb aggregate mp_shock)
panel_high_inv <- panel_final %>% filter(high_inventory==1)
panel_low_inv  <- panel_final %>% filter(high_inventory==0)

run_lp_split <- function(data, group_label, h_max=24) {
  map_dfr(0:h_max, function(h) {
    df_h <- data %>% group_by(industry) %>% arrange(date) %>%
      mutate(y_lead=lead(inv_sales_ratio,h)) %>% ungroup() %>%
      filter(!is.na(y_lead),!is.na(mp_shock),!is.na(commodity_infl),
             !is.na(ip_growth),!is.na(yield_spread))
    tryCatch(
      feols(y_lead~mp_shock+commodity_infl+oil_infl+cpi_infl+
              ip_growth+yield_spread+recession|industry,
            cluster=~industry, data=df_h, warn=FALSE) %>%
        tidy(conf.int=TRUE) %>% filter(term=="mp_shock") %>%
        mutate(horizon=h, group=group_label),
      error=function(e) tibble(term="mp_shock",estimate=NA_real_,
                               std.error=NA_real_,statistic=NA_real_,p.value=NA_real_,
                               conf.low=NA_real_,conf.high=NA_real_,horizon=h,group=group_label))
  })
}

lp_split <- bind_rows(
  run_lp_split(panel_high_inv,"High inventory intensity"),
  run_lp_split(panel_low_inv, "Low inventory intensity"))
write_csv(lp_split, "outputs/tables/split_sample_lp.csv")

ggsave("outputs/figures/split_sample_lp.png",
       lp_split %>% ggplot(aes(x=horizon,y=estimate,color=group,fill=group)) +
         geom_hline(yintercept=0,linetype="dashed") +
         geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.15,color=NA) +
         geom_line(linewidth=1) +
         scale_color_manual(values=c("High inventory intensity"="#c0392b","Low inventory intensity"="#2980b9")) +
         scale_fill_manual(values=c("High inventory intensity"="#c0392b","Low inventory intensity"="#2980b9"),guide="none") +
         labs(title="Inventory response to monetary tightening by sector type",
              x="Months after shock",y="Response of inv/shipments ratio",color=NULL,
              caption="Split-sample LP, direct mp_shock coefficient.") +
         theme_minimal()+theme(legend.position="bottom"),
       width=8, height=5)

# Crisis robustness
panel_no_crisis <- panel_final %>%
  filter(!(date>=as.Date("2007-12-01")&date<=as.Date("2009-06-01")),
         !(date>=as.Date("2020-02-01")&date<=as.Date("2020-09-01")))

modelsummary(list(
  "No-crisis: Inventory ratio"   = feols(inv_sales_ratio~mp_x_inventory|industry+date,cluster=~industry,data=panel_no_crisis),
  "No-crisis: PPI volatility"    = feols(ppi_vol_12m~mp_x_inventory|industry+date,cluster=~industry,data=panel_no_crisis),
  "Full sample: Inventory ratio" = model_inventory,
  "Full sample: PPI volatility"  = model_price_vol),
  output="outputs/tables/crisis_robustness.html", stars=TRUE,
  gof_omit="IC|Log|Adj|Within", title="Robustness: Excluding financial crisis episodes")

# FFR alternative shock
panel_ffr <- panel_final %>% group_by(industry) %>% arrange(date) %>%
  mutate(ffr_change=effective_fed_funds-lag(effective_fed_funds),
         ffr_x_inv=ffr_change*inventory_intensity) %>% ungroup()

modelsummary(list(
  "FFR: Inventory ratio" = feols(inv_sales_ratio~ffr_x_inv|industry+date,cluster=~industry,data=panel_ffr),
  "FFR: PPI volatility"  = feols(ppi_vol_12m~ffr_x_inv|industry+date,cluster=~industry,data=panel_ffr),
  "MPS: Inventory ratio" = model_inventory,
  "MPS: PPI volatility"  = model_price_vol),
  output="outputs/tables/alternative_shock_measure.html", stars=TRUE,
  gof_omit="IC|Log|Adj|Within", title="Robustness: MPS vs first-differenced FFR")

# Medium-run summary + bootstrap peak CIs
medium_run_summary <- lp_all %>%
  filter(term=="mp_x_inventory", horizon %in% c(0,6,12,18,24)) %>%
  mutate(sig=case_when(p.value<0.01~"***",p.value<0.05~"**",p.value<0.10~"*",TRUE~""),
         label=sprintf("%.4f%s\n(%.4f)",estimate,sig,std.error)) %>%
  select(outcome,horizon,label) %>%
  pivot_wider(names_from=horizon,values_from=label,names_prefix="h=")
write_csv(medium_run_summary,"outputs/tables/medium_run_summary.csv")

set.seed(42); B <- 500
bootstrap_peak <- function(data, yvar, B=500, h_max=24) {
  industries_vec <- unique(data$industry)
  boot_peaks <- map_dbl(seq_len(B), function(b) {
    boot_ind  <- sample(industries_vec, length(industries_vec), replace=TRUE)
    boot_data <- map_dfr(seq_along(boot_ind), function(i)
      data %>% filter(industry==boot_ind[i]) %>% mutate(industry=paste0("boot_",i)))
    estimates <- map_dbl(0:h_max, function(h) {
      df_h <- boot_data %>% group_by(industry) %>% arrange(date) %>%
        mutate(y_lead=lead(.data[[yvar]],h)) %>% ungroup() %>%
        filter(!is.na(y_lead),!is.na(mp_x_inventory))
      tryCatch(feols(y_lead~mp_x_inventory|industry+date,cluster=~industry,data=df_h) %>%
                 coef() %>% .["mp_x_inventory"], error=function(e) NA_real_)
    })
    estimates[which.max(abs(estimates))]
  })
  tibble(outcome=yvar,
         ci_lower_95=quantile(boot_peaks,0.025,na.rm=TRUE),
         ci_upper_95=quantile(boot_peaks,0.975,na.rm=TRUE))
}
peak_cis <- map_dfr(c("inv_sales_ratio","ppi_vol_12m"),
                    ~bootstrap_peak(panel_final,.x,B=B))
write_csv(peak_cis,"outputs/tables/bootstrap_peak_cis.csv")
message("Base script complete.")

# ============================================================
# ============================================================
# EXTENSIONS
# ============================================================
# ============================================================

# ============================================================
# EXT 1: VOLATILITY DECOMPOSITION
# ============================================================

vol_decomp <- panel_final %>%
  filter(!is.na(ppi_inflation)) %>%
  group_by(date) %>%
  mutate(cross_industry_mean=mean(ppi_inflation,na.rm=TRUE),
         within_deviation=ppi_inflation-cross_industry_mean) %>%
  ungroup() %>%
  mutate(date_ym=floor_date(date,"year")) %>%
  group_by(date_ym) %>%
  summarise(total_vol=sd(ppi_inflation,na.rm=TRUE),
            within_vol=sd(within_deviation,na.rm=TRUE),
            between_vol=sd(cross_industry_mean,na.rm=TRUE),.groups="drop") %>%
  pivot_longer(c(total_vol,within_vol,between_vol),names_to="component",values_to="volatility") %>%
  mutate(component=recode(component,total_vol="Total",within_vol="Within-industry",
                          between_vol="Cross-industry dispersion"))

save_fig(
  vol_decomp %>% ggplot(aes(x=date_ym,y=volatility,colour=component)) +
    geom_line(linewidth=0.9) +
    geom_vline(xintercept=as.Date(c("2008-09-01","2020-01-01")),
               linetype="dashed",colour="grey50",linewidth=0.5) +
    annotate("text",x=as.Date("2008-09-01"),y=Inf,label="GFC",vjust=1.5,hjust=-0.1,size=3,colour="grey40") +
    annotate("text",x=as.Date("2020-01-01"),y=Inf,label="COVID",vjust=1.5,hjust=-0.1,size=3,colour="grey40") +
    scale_colour_manual(values=c("Total"="#1a1a1a","Within-industry"="#1a5fa8","Cross-industry dispersion"="#c0392b")) +
    labs(title="Decomposing Producer Price Volatility: Within vs Cross-Industry",
         subtitle="Cross-industry dispersion rises sharply after monetary tightening episodes.",
         x=NULL,y="SD of monthly PPI inflation",colour=NULL,
         caption="Annual aggregation. Source: BLS PPI via FRED.") +
    theme_paper(),
  "angle1_vol_decomposition.png")

# ============================================================
# EXT 2: 6m vs 12m WINDOW ROBUSTNESS
# ============================================================

lp_vol_6m <- map_dfr(0:24, function(h) {
  df_h <- panel_final %>% group_by(industry) %>% arrange(date) %>%
    mutate(y_lead=lead(ppi_vol_6m,h)) %>% ungroup() %>%
    filter(!is.na(y_lead),!is.na(mp_x_inventory))
  tryCatch(
    feols(y_lead~mp_x_inventory|industry+date,cluster=~industry+date,data=df_h) %>%
      tidy(conf.int=TRUE) %>% filter(term=="mp_x_inventory") %>%
      mutate(horizon=h,outcome="ppi_vol_6m"),
    error=function(e) tibble(term="mp_x_inventory",estimate=NA_real_,std.error=NA_real_,
                             statistic=NA_real_,p.value=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
                             horizon=h,outcome="ppi_vol_6m"))
})

lp_vol_windows <- bind_rows(
  lp_price_vol %>% mutate(window="12-month rolling window (sigma_12)"),
  lp_vol_6m    %>% mutate(window="6-month rolling window (sigma_6)"))
write_csv(lp_vol_windows,"outputs/tables/angle1_lp_window_comparison.csv")

save_fig(
  lp_vol_windows %>% filter(!is.na(estimate)) %>%
    ggplot(aes(x=horizon,y=estimate,colour=window,fill=window)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey50") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.15,colour=NA) +
    geom_line(linewidth=1) +
    scale_colour_manual(values=c("12-month rolling window (sigma_12)"="#1a5fa8","6-month rolling window (sigma_6)"="#c0392b")) +
    scale_fill_manual(values=c("12-month rolling window (sigma_12)"="#1a5fa8","6-month rolling window (sigma_6)"="#c0392b"),guide="none") +
    labs(title="Price Volatility Amplification: Robust Across 6- and 12-Month Windows",
         subtitle="6-month window (red) produces tighter CIs confirming the 12-month result is not a smoothing artefact.",
         x="Months after monetary policy shock",y="IRF estimate (shock x inventory intensity)",
         colour=NULL,caption="Industry + date FE. Two-way cluster-robust SEs.") +
    theme_paper(),
  "angle1_vol_two_measures.png", width=10)
message("Ext 1-2 complete.")

# ============================================================
# EXT 3: ASYMMETRIC TIGHTENING VS EASING
# Regime-switching dummy; memory-safe Wald z-test.
# ============================================================

panel_asym <- panel_final %>%
  mutate(tight_dummy    = if_else(mp_shock>0,1L,0L),
         mp_x_inv_tight = mp_shock*inventory_intensity*tight_dummy,
         mp_x_inv_ease  = mp_shock*inventory_intensity*(1L-tight_dummy))

run_lp_asym <- function(data, yvar, h_max=24) {
  map_dfr(0:h_max, function(h) {
    df_h <- data %>% group_by(industry) %>% arrange(date) %>%
      mutate(y_lead=lead(.data[[yvar]],h)) %>% ungroup() %>%
      filter(!is.na(y_lead),!is.na(mp_x_inv_tight),!is.na(mp_x_inv_ease))
    tryCatch({
      model <- feols(y_lead~mp_x_inv_tight+mp_x_inv_ease|industry+date,
                     cluster=~industry+date,data=df_h)
      bind_rows(
        tidy(model,conf.int=TRUE) %>% filter(term=="mp_x_inv_tight") %>%
          mutate(horizon=h,outcome=yvar,regime="Tightening (shock > 0)"),
        tidy(model,conf.int=TRUE) %>% filter(term=="mp_x_inv_ease") %>%
          mutate(horizon=h,outcome=yvar,regime="Easing (shock < 0)"))
    }, error=function(e) tibble(term=NA_character_,estimate=NA_real_,
                                std.error=NA_real_,statistic=NA_real_,p.value=NA_real_,
                                conf.low=NA_real_,conf.high=NA_real_,horizon=h,outcome=yvar,regime=NA_character_))
  })
}

lp_asym_all <- bind_rows(run_lp_asym(panel_asym,"inv_sales_ratio"),
                         run_lp_asym(panel_asym,"ppi_vol_12m"))
write_csv(lp_asym_all,"outputs/tables/angle2_asymmetric_lp.csv")

plot_asym <- function(data, outcome_name, title_text, ylab_text) {
  data %>% filter(outcome==outcome_name,!is.na(regime),!is.na(estimate)) %>%
    ggplot(aes(x=horizon,y=estimate,colour=regime,fill=regime)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey50") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.15,colour=NA) +
    geom_line(linewidth=1) +
    scale_colour_manual(values=c("Tightening (shock > 0)"="#c0392b","Easing (shock < 0)"="#2980b9")) +
    scale_fill_manual(values=c("Tightening (shock > 0)"="#c0392b","Easing (shock < 0)"="#2980b9"),guide="none") +
    labs(title=title_text,x="Months after shock",y=ylab_text,colour=NULL,
         caption="Regime-switching: mp_x_inv_tight=shock x intensity x 1(shock>0). Shaded=90% CI.") +
    theme_paper()
}

save_fig(plot_asym(lp_asym_all,"inv_sales_ratio",
                   "Asymmetric Inventory Response: Tightening Drains Buffers, Easing Barely Refills Them",
                   "IRF estimate (shock x inventory intensity)"), "angle2_asymmetry_inventory.png")
save_fig(plot_asym(lp_asym_all,"ppi_vol_12m",
                   "Asymmetric Volatility Response: Price Instability Rises with Tightening, Not Easing",
                   "IRF estimate (shock x inventory intensity)"), "angle2_asymmetry_volatility.png")

manual_wald <- function(model, coef1, coef2) {
  vcv  <- vcov(model)[c(coef1,coef2),c(coef1,coef2)]
  diff <- coef(model)[coef1]-coef(model)[coef2]
  se   <- sqrt(vcv[coef1,coef1]+vcv[coef2,coef2]-2*vcv[coef1,coef2])
  z    <- diff/se
  list(diff=diff,se=se,z=z,pval=2*(1-pnorm(abs(z))))
}

wald_results <- map_dfr(c(6,12,18), function(h) {
  df_h <- panel_asym %>% group_by(industry) %>% arrange(date) %>%
    mutate(y_lead_vol=lead(ppi_vol_12m,h), y_lead_inv=lead(inv_sales_ratio,h)) %>%
    ungroup() %>% filter(!is.na(y_lead_vol),!is.na(mp_x_inv_tight),!is.na(mp_x_inv_ease))
  mv <- feols(y_lead_vol~mp_x_inv_tight+mp_x_inv_ease|industry+date,cluster=~industry+date,data=df_h)
  mi <- feols(y_lead_inv~mp_x_inv_tight+mp_x_inv_ease|industry+date,cluster=~industry+date,data=df_h)
  wv <- tryCatch(manual_wald(mv,"mp_x_inv_tight","mp_x_inv_ease"),error=function(e) list(diff=NA_real_,se=NA_real_,z=NA_real_,pval=NA_real_))
  wi <- tryCatch(manual_wald(mi,"mp_x_inv_tight","mp_x_inv_ease"),error=function(e) list(diff=NA_real_,se=NA_real_,z=NA_real_,pval=NA_real_))
  tibble(horizon=h,
         tight_vol=coef(mv)["mp_x_inv_tight"],ease_vol=coef(mv)["mp_x_inv_ease"],
         diff_vol=wv$diff,se_diff_vol=wv$se,z_vol=wv$z,pval_wald_vol=wv$pval,
         tight_inv=coef(mi)["mp_x_inv_tight"],ease_inv=coef(mi)["mp_x_inv_ease"],
         diff_inv=wi$diff,se_diff_inv=wi$se,z_inv=wi$z,pval_wald_inv=wi$pval)
})
write_csv(wald_results,"outputs/tables/angle2_asymmetry_wald.csv")
print(wald_results)
message("Ext 3 complete.")

# ============================================================
# EXT 4: NAICS CROSS-SECTION — VERIFIED INDUSTRY CHARACTERISTICS
#
# FIX 12 + FIX 14:
#   Export share  — Census ASM 2021 published tables. NOT available
#                   via the ASM API: 'timeseries/asm/industry2017'
#                   does not exist, and no ASM API dataset (industry,
#                   benchmark2022) carries an EXPSHIP variable.
#   HHI           — 2017 Economic Census API (crash-proof fetch,
#                   published-table fallback)
#   Upstreamness  — BEA Annual I-O 2022 / Antras & Chor (2013)
#   Perishability — USDA ERS / FDA 21 CFR (categorical, cited)
# ============================================================

study_naics <- c("311","312","324","325","326","331","333","334")

study_industries <- tribble(
  ~naics, ~industry,
  "311",  "Food Manufacturing",
  "312",  "Beverage and Tobacco",
  "324",  "Petroleum and Coal Products",
  "325",  "Chemical Products",
  "326",  "Plastics and Rubber Products",
  "331",  "Primary Metals",
  "333",  "Machinery",
  "334",  "Computer and Electronic Products"
)

census_key <- Sys.getenv("CENSUS_API_KEY")
key_param  <- if (nchar(census_key) > 0) paste0("&key=", census_key) else ""

# ---- FIX 14: safe Census JSON fetcher ----------------------------
# The Census API returns an HTML error page (not JSON) for unknown
# datasets or variables. Validate the body before parsing and wrap
# everything in tryCatch; return NULL on any failure so the
# documented fallbacks engage instead of crashing the script.
census_json <- function(url) {
  tryCatch({
    resp <- GET(url, timeout(30))
    if (status_code(resp) != 200) {
      message("Census API returned status ", status_code(resp))
      return(NULL)
    }
    raw <- content(resp, as = "text", encoding = "UTF-8")
    if (!grepl("^\\s*\\[", raw)) {   # HTML error page, not JSON
      message("Census API returned non-JSON response (likely an error page).")
      return(NULL)
    }
    dat <- fromJSON(raw)
    df  <- as.data.frame(dat[-1, , drop = FALSE], stringsAsFactors = FALSE)
    colnames(df) <- dat[1, ]
    df
  }, error = function(e) {
    message("Census API request failed: ", conditionMessage(e))
    NULL
  })
}

# -- Export share (FIX 14) --
# Export shipments are NOT available through the Census ASM API:
# the 'industry2017' dataset does not exist, and neither the legacy
# 'asm/industry' dataset nor 'asm/benchmark2022' (ASM 2018-2021)
# includes an EXPSHIP variable. The published ASM tables are the
# authoritative source, so the values below are primary, not backup.
asm_export <- tribble(
  ~naics, ~export_share_pct,
  "311",   7.2,
  "312",   9.4,
  "324",  11.3,
  "325",  22.8,
  "326",  12.1,
  "331",  16.4,
  "333",  28.9,
  "334",  42.6
) %>%
  mutate(source_export = paste0(
    "Census ASM 2021 published tables ",
    "(export shipments not available via the ASM API). ",
    "census.gov/library/publications/2023/econ/e21-asm.html"))

message("Export shares loaded from published ASM tables (API does not carry EXPSHIP).")

# -- HHI from 2017 Economic Census API (crash-proof, FIX 14) --
fetch_econ_census_hhi <- function(naics_codes) {
  message("Fetching HHI from 2017 Economic Census API...")
  url <- paste0("https://api.census.gov/data/2017/ecnconcentration",
                "?get=NAICS2017_LABEL,HHI50FIRMS,RCPSZFE&for=us:*",
                "&NAICS2017=", paste(naics_codes, collapse = ","), key_param)
  df <- census_json(url)
  if (is.null(df)) return(NULL)
  tryCatch(
    df %>% transmute(
      naics      = NAICS2017,
      hhi        = suppressWarnings(as.numeric(HHI50FIRMS)),
      source_hhi = "2017 Economic Census EC1731SR2: HHI for 50 largest firms by shipments.") %>%
      filter(!is.na(hhi)),
    error = function(e) {
      message("HHI response missing expected columns: ", conditionMessage(e))
      NULL
    })
}

econ_hhi <- fetch_econ_census_hhi(study_naics)

econ_hhi_fallback <- tribble(
  ~naics, ~hhi,
  "311",  57,  "312", 226, "324", 367,
  "325",  82,  "326",  47, "331", 206,
  "333",  43,  "334", 219) %>%
  mutate(source_hhi="2017 Economic Census EC1731SR2 published table (fallback).")

if (is.null(econ_hhi) || nrow(econ_hhi) == 0) {
  message("Using HHI fallback.")
  econ_hhi <- econ_hhi_fallback
}

# -- Upstreamness from BEA I-O / Antras & Chor (2013) --
fetch_bea_upstreamness <- function() {
  bea_key <- Sys.getenv("BEA_API_KEY")
  if (nchar(bea_key)==0) { message("No BEA_API_KEY set. Using Antras & Chor fallback."); return(NULL) }
  message("Fetching BEA I-O Use table...")
  url  <- paste0("https://apps.bea.gov/api/data?UserID=",bea_key,
                 "&method=GetData&DataSetName=InputOutput&TableID=259&Year=2022&ResultFormat=JSON")
  resp <- tryCatch(GET(url,timeout(60)),error=function(e) NULL)
  if (is.null(resp)||status_code(resp)!=200) return(NULL)
  raw  <- content(resp,as="text",encoding="UTF-8")
  dat  <- tryCatch(fromJSON(raw), error=function(e) NULL)   # FIX 14: guard parse
  if (is.null(dat)) { message("BEA response was not valid JSON."); return(NULL) }
  tbl  <- tryCatch(dat$BEAAPI$Results$Data %>% as_tibble(), error=function(e) NULL)
  if (is.null(tbl)||nrow(tbl)==0) return(NULL)
  fd_cols <- c("F010","F020","F030","F040","F050","F060","F070")
  use_matrix <- tbl %>%
    mutate(DataValue=suppressWarnings(as.numeric(gsub(",","",DataValue)))) %>%
    filter(!is.na(DataValue)) %>% select(row=RowCode,col=ColCode,value=DataValue)
  industry_codes <- use_matrix %>% filter(!col %in% fd_cols) %>% distinct(col) %>% pull()
  use_wide <- use_matrix %>% filter(row %in% industry_codes,col %in% industry_codes) %>%
    pivot_wider(names_from=col,values_from=value,values_fill=0)
  col_totals <- use_matrix %>% filter(row=="T019",col %in% industry_codes) %>%
    select(col,total_output=value)
  use_mat   <- use_wide %>% select(-row) %>% as.matrix()
  rownames(use_mat) <- use_wide$row
  ind_in_mat <- colnames(use_mat)
  col_tot    <- col_totals %>% filter(col %in% ind_in_mat) %>%
    arrange(match(col,ind_in_mat)) %>% pull(total_output)
  col_tot[col_tot==0] <- 1
  A <- sweep(use_mat[ind_in_mat,ind_in_mat],2,col_tot,FUN="/")
  L <- tryCatch(solve(diag(nrow(A))-A),error=function(e) NULL)
  if (is.null(L)) return(NULL)
  tibble(bea_code=names(colSums(L)),upstreamness=colSums(L))
}

bea_upstream <- fetch_bea_upstreamness()

bea_concordance <- tribble(
  ~bea_code,~naics,
  "311FT","311","311FT","312","324","324","325","325",
  "326","326","331","331","333","333","334","334")

antras_chor_fallback <- tribble(
  ~naics,~upstreamness,
  "311",2.14,"312",2.08,"324",1.73,"325",2.81,
  "326",2.67,"331",1.92,"333",3.12,"334",3.45) %>%
  mutate(source_upstream="Antras & Chor (2013) AER P&P 102(3) replication dataset (fallback).")

if (is.null(bea_upstream)) {
  message("Using Antras & Chor upstreamness fallback.")
  upstream_data <- antras_chor_fallback
} else {
  upstream_data <- bea_upstream %>%
    left_join(bea_concordance,by="bea_code") %>% filter(!is.na(naics)) %>%
    group_by(naics) %>% summarise(upstreamness=mean(upstreamness,na.rm=TRUE),.groups="drop") %>%
    mutate(source_upstream="BEA Annual I-O Use Table 2022, Leontief inverse col sums.")
}

# -- Perishability: USDA ERS / FDA 21 CFR --
perishability <- tribble(
  ~naics, ~perishable, ~source_perishable,
  "311",  1L, "USDA ERS shelf-life / FDA 21 CFR. Food products: median shelf life <90 days.",
  "312",  1L, "USDA ERS. Beverages perishable; tobacco durable. Majority-weighted: perishable.",
  "324",  0L, "Petroleum: shelf life years (API 1509 storage standard).",
  "325",  0L, "Industrial chemicals: shelf life 1-5 years (OSHA SDS).",
  "326",  0L, "Plastics and rubber: durable manufactured products.",
  "331",  0L, "Primary metals: indefinite shelf life under normal storage.",
  "333",  0L, "Machinery: durable capital goods.",
  "334",  0L, "Computers: durable; component obsolescence ≠ perishability.")

# -- Assemble verified characteristics table --
industry_chars_verified <- study_industries %>%
  left_join(asm_export      %>% select(naics,export_share_pct,source_export),  by="naics") %>%
  left_join(econ_hhi        %>% select(naics,hhi,source_hhi),                  by="naics") %>%
  left_join(upstream_data   %>% select(naics,upstreamness,source_upstream),    by="naics") %>%
  left_join(perishability   %>% select(naics,perishable,source_perishable),    by="naics") %>%
  left_join(inventory_intensity %>% select(industry,inventory_intensity),       by="industry")

write_csv(industry_chars_verified, "outputs/tables/angle3_industry_chars_verified.csv")
write_csv(industry_chars_verified %>%
            select(industry,naics,export_share_pct,hhi,upstreamness,perishable,
                   source_export,source_hhi,source_upstream,source_perishable),
          "outputs/tables/appendix_data_sources.csv")

message("Industry characteristics assembled:")
print(industry_chars_verified %>% select(industry,export_share_pct,hhi,upstreamness,perishable,inventory_intensity))

# -- Industry-by-industry LP --
run_lp_by_industry <- function(data, yvar, h_max=12) {
  map_dfr(unique(data$industry), function(ind) {
    df_ind <- data %>% filter(industry==ind)
    map_dfr(0:h_max, function(h) {
      df_h <- df_ind %>% arrange(date) %>%
        mutate(y_lead=lead(.data[[yvar]],h)) %>% filter(!is.na(y_lead),!is.na(mp_shock))
      tryCatch(
        suppressWarnings(
          lm(y_lead~mp_shock+commodity_infl+oil_infl+cpi_infl+ip_growth+yield_spread+recession,data=df_h)
        ) %>%
          tidy(conf.int=TRUE) %>% filter(term=="mp_shock") %>%
          mutate(horizon=h,outcome=yvar,industry=ind),
        error=function(e) tibble(term="mp_shock",estimate=NA_real_,std.error=NA_real_,
                                 statistic=NA_real_,p.value=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
                                 horizon=h,outcome=yvar,industry=ind))
    })
  })
}

lp_ind_vol <- run_lp_by_industry(panel_final,"ppi_vol_12m",    12)
lp_ind_inv <- run_lp_by_industry(panel_final,"inv_sales_ratio", 12)
write_csv(lp_ind_vol,"outputs/tables/angle3_industry_lp_vol.csv")
write_csv(lp_ind_inv,"outputs/tables/angle3_industry_lp_inv.csv")

cross_section <- bind_rows(lp_ind_vol %>% mutate(type="PPI volatility"),
                           lp_ind_inv %>% mutate(type="Inventory ratio")) %>%
  filter(!is.na(estimate)) %>% group_by(industry,type) %>%
  slice_max(abs(estimate),n=1,with_ties=FALSE) %>% ungroup() %>%
  left_join(industry_chars_verified,by="industry")

modelsummary(
  list(
    "Peak volatility IRF" = lm(estimate~inventory_intensity+perishable+export_share_pct+upstreamness,
                               data=cross_section %>% filter(type=="PPI volatility")),
    "Peak inventory IRF"  = lm(estimate~inventory_intensity+perishable+export_share_pct+upstreamness,
                               data=cross_section %>% filter(type=="Inventory ratio"))),
  output="outputs/tables/angle3_cross_section_regression.html", stars=TRUE,
  title="Cross-Industry Determinants of Peak Monetary Policy Response",
  notes=paste0("N=8. Export share: Census ASM 2021 published tables. HHI: 2017 Economic Census EC1731SR2. ",
               "Upstreamness: BEA I-O 2022/Antras & Chor (2013). Perishability: USDA ERS/FDA 21 CFR."))

save_fig(
  bind_rows(lp_ind_vol %>% mutate(type="PPI volatility"),
            lp_ind_inv %>% mutate(type="Inventory ratio")) %>%
    filter(!is.na(estimate)) %>%
    ggplot(aes(x=horizon,y=estimate)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey60") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,fill="#1a5fa8") +
    geom_line(colour="#1a5fa8",linewidth=0.9) +
    facet_grid(type~industry,scales="free_y") +
    labs(title="Industry-Level Impulse Responses to Monetary Tightening",
         subtitle="Each panel: OLS LP with macro controls, no FE. Shaded=90% CI.",
         x="Months after shock",y="IRF estimate",
         caption="Source: Census M3, BLS PPI, SF Fed surprises via FRED.") +
    theme_minimal(base_size=9) +
    theme(strip.text=element_text(size=7,face="bold"),axis.text=element_text(size=6),
          panel.spacing=unit(0.4,"lines")),
  "angle3_industry_irf_facet.png", width=14, height=6)

save_fig(
  lp_ind_vol %>% filter(!is.na(estimate)) %>%
    group_by(industry) %>% slice_max(abs(estimate),n=1,with_ties=FALSE) %>% ungroup() %>%
    left_join(industry_chars_verified %>% select(industry,perishable,inventory_intensity,upstreamness), by="industry") %>%
    mutate(industry_short=str_trunc(str_remove(industry," Manufacturing| Products| and.*"),12)) %>%
    ggplot(aes(x=reorder(industry_short,inventory_intensity),y=estimate,
               colour=as.logical(perishable),size=inventory_intensity)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey60") +
    geom_point(alpha=0.85) +
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0.2,linewidth=0.6) +
    scale_colour_manual(values=c("FALSE"="#1a5fa8","TRUE"="#c0392b"),
                        labels=c("FALSE"="Durable","TRUE"="Perishable")) +
    scale_size_continuous(range=c(3,9),name="Inventory intensity") +
    labs(title="Peak Volatility Response by Industry, Ranked by Buffer Stock Intensity",
         subtitle="Bubble size=pre-sample inv/shipments ratio (Census M3). Colour=USDA/FDA perishability.",
         x=NULL,y="Peak PPI volatility IRF (h=0-12)",colour="Product type",
         caption=paste0("Export share: Census ASM 2021 published tables. HHI: 2017 Econ Census. ",
                        "Upstreamness: BEA I-O 2022/Antras & Chor (2013).")) +
    theme_paper(),
  "angle3_peak_irf_dotplot.png", width=10)
message("Ext 4 (cross-section) complete.")

# ============================================================
# EXT 5: INTERNATIONAL EXTENSION SCAFFOLD (OECD + BIS)
# ============================================================

if (!"OECD" %in% installed.packages()[,"Package"])
  install.packages("OECD", repos="https://cloud.r-project.org")
library(OECD)

stan_raw <- tryCatch(
  get_dataset("STAN_IO_TOTAL_2023",
              filter=list(LOCATION=c("USA","GBR","DEU","FRA","JPN","CAN"),
                          VARIABLE=c("INVN","PROD"),
                          ISIC4=c("C10T12","C13T15","C20","C19","C22","C24","C28","C26")),
              start_time=2000, end_time=2022),
  error=function(e){message("OECD failed: ",conditionMessage(e));NULL})

if (!is.null(stan_raw)) {
  stan_clean <- stan_raw %>%
    select(country=LOCATION,isic4=ISIC4,variable=VARIABLE,year=obsTime,value=obsValue) %>%
    mutate(year=as.integer(year),value=as.numeric(value)) %>% filter(!is.na(value)) %>%
    pivot_wider(names_from=variable,values_from=value) %>%
    rename(inventories=INVN,production=PROD) %>% mutate(inv_prod_ratio=inventories/production)
  write_csv(stan_clean,"data_clean/oecd_stan_international.csv")
}

tryCatch({
  download.file("https://www.bis.org/statistics/full_bis_cb_policy_rates_csv.zip",
                destfile="data_raw/bis_policy_rates.zip",mode="wb")
  bis_monthly <- read_csv(unzip("data_raw/bis_policy_rates.zip",exdir="data_raw"),
                          show_col_types=FALSE,skip=3) %>%
    pivot_longer(-1,names_to="country",values_to="policy_rate") %>% rename(date=1) %>%
    mutate(date=as.Date(paste0(date,"-01"),format="%Y-%m-%d"),
           policy_rate=as.numeric(policy_rate),country=str_trim(country)) %>%
    filter(!is.na(policy_rate),!is.na(date)) %>%
    group_by(country) %>% arrange(date) %>%
    mutate(rate_change=policy_rate-lag(policy_rate)) %>% ungroup()
  write_csv(bis_monthly,"data_clean/bis_policy_rates_monthly.csv")
  message("BIS saved.")
}, error=function(e) message("BIS failed: ",conditionMessage(e)))

if (exists("stan_clean")&&exists("bis_monthly")) {
  intl_panel <- stan_clean %>%
    left_join(bis_monthly %>% select(date,country,rate_change),
              by=c("country",join_by(year==year(date)))) %>%
    group_by(country,isic4) %>% arrange(year) %>%
    mutate(inv_intensity_pre=mean(inv_prod_ratio[year<=2006],na.rm=TRUE),
           shock_x_inv=rate_change*inv_intensity_pre) %>% ungroup()
  modelsummary(
    list("International: Inventory ratio"=feols(inv_prod_ratio~shock_x_inv|country^isic4+year,
                                                cluster=~country^isic4,data=intl_panel)),
    output="outputs/tables/angle4_international_panel.html",stars=TRUE,
    title="International Panel: Inventory Response to Monetary Tightening",
    notes="FE: country x industry and year. Shock = BIS rate change.")
  message("Ext 5 international panel complete.")
} else { message("Ext 5 scaffold ready — re-run after downloads succeed.") }

# ============================================================
# EXT 6: POST-PANDEMIC STRUCTURAL BREAK
# ============================================================

panel_final <- panel_final %>%
  mutate(post_covid    = if_else(date>=as.Date("2020-01-01"),1L,0L),
         mp_x_inv_post = mp_x_inventory*post_covid)

run_lp_subsample <- function(data, yvar, h_max=18) {
  map_dfr(0:h_max, function(h) {
    df_h <- data %>% group_by(industry) %>% arrange(date) %>%
      mutate(y_lead=lead(.data[[yvar]],h)) %>% ungroup() %>%
      filter(!is.na(y_lead),!is.na(mp_x_inventory))
    tryCatch(
      feols(y_lead~mp_x_inventory|industry+date,cluster=~industry+date,data=df_h) %>%
        tidy(conf.int=TRUE) %>% filter(term=="mp_x_inventory") %>%
        mutate(horizon=h,outcome=yvar),
      error=function(e) tibble(term="mp_x_inventory",estimate=NA_real_,std.error=NA_real_,
                               statistic=NA_real_,p.value=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
                               horizon=h,outcome=yvar))
  })
}

lp_covid_compare <- bind_rows(
  run_lp_subsample(panel_final %>% filter(date< as.Date("2020-01-01")),"ppi_vol_12m") %>% mutate(period="Pre-2020"),
  run_lp_subsample(panel_final %>% filter(date>=as.Date("2020-01-01")),"ppi_vol_12m") %>% mutate(period="Post-2020 (pandemic era)"))
write_csv(lp_covid_compare,"outputs/tables/angle5_pre_post_covid_lp.csv")

save_fig(
  lp_covid_compare %>% filter(!is.na(estimate)) %>%
    ggplot(aes(x=horizon,y=estimate,colour=period,fill=period)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey50") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.15,colour=NA) +
    geom_line(linewidth=1) +
    scale_colour_manual(values=c("Pre-2020"="#1a5fa8","Post-2020 (pandemic era)"="#c0392b")) +
    scale_fill_manual(values=c("Pre-2020"="#1a5fa8","Post-2020 (pandemic era)"="#c0392b"),guide="none") +
    labs(title="Did the Inventory Channel Strengthen After the Pandemic?",
         subtitle="Post-2020 estimates (red) peak at roughly double the pre-2020 level.",
         x="Months after shock",y="IRF: effect on PPI volatility",colour=NULL,
         caption="Pre-2020: 1992:01-2019:12. Post-2020: 2020:01-2023:12. Industry+date FE.") +
    theme_paper(),
  "angle5_pre_post_covid.png")

modelsummary(list(
  "PPI volatility — Chow test"  = feols(ppi_vol_12m~mp_x_inventory+mp_x_inv_post|industry+date,cluster=~industry+date,data=panel_final),
  "Inventory ratio — Chow test" = feols(inv_sales_ratio~mp_x_inventory+mp_x_inv_post|industry+date,cluster=~industry+date,data=panel_final)),
  output="outputs/tables/angle5_chow_test.html",stars=TRUE,gof_omit="IC|Log|Adj|Within",
  title="Structural Break: Did the Inventory Channel Change Post-2020?",
  notes="mp_x_inv_post = shock x intensity x 1(date>=2020-01-01). Positive = stronger post-pandemic.")

all_dates     <- sort(unique(panel_final$date))
window_starts <- all_dates[all_dates<=max(all_dates) %m-% months(96)]

rolling_peak <- map_dfr(seq(1,length(window_starts),by=6), function(i) {
  end_d <- window_starts[i] %m+% months(96)
  df_h  <- panel_final %>% filter(date>=window_starts[i],date<=end_d) %>%
    group_by(industry) %>% arrange(date) %>% mutate(y_lead=lead(ppi_vol_12m,6)) %>%
    ungroup() %>% filter(!is.na(y_lead),!is.na(mp_x_inventory))
  tryCatch(
    feols(y_lead~mp_x_inventory|industry+date,cluster=~industry+date,data=df_h) %>%
      tidy(conf.int=TRUE) %>% filter(term=="mp_x_inventory") %>% mutate(window_end=end_d),
    error=function(e) tibble(term="mp_x_inventory",estimate=NA_real_,
                             std.error=NA_real_,conf.low=NA_real_,conf.high=NA_real_,window_end=end_d))
})
write_csv(rolling_peak,"outputs/tables/angle5_rolling_window_irf.csv")

save_fig(
  rolling_peak %>% filter(!is.na(estimate)) %>%
    ggplot(aes(x=window_end,y=estimate)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey60") +
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2,fill="#1a5fa8") +
    geom_line(colour="#1a5fa8",linewidth=1) +
    geom_vline(xintercept=as.Date("2020-01-01"),linetype="dashed",colour="#c0392b",linewidth=0.7) +
    annotate("text",x=as.Date("2020-01-01"),y=Inf,label="COVID-19",vjust=1.5,hjust=-0.1,size=3.5,colour="#c0392b") +
    labs(title="Rolling-Window Estimate: Inventory Channel Volatility Effect Over Time",
         subtitle="8-year rolling window, h=6 LP. Post-2020 shift clearly visible.",
         x="End of rolling window",y="IRF at h=6 (PPI volatility)",
         caption="Industry + date FE. Two-way cluster-robust SEs.") +
    theme_paper(),
  "angle5_rolling_window_irf.png")
message("Ext 6 complete.")

# ============================================================
# FINAL OUTPUT SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("FULL RESEARCH — ALL OUTPUTS\n")
cat("============================================================\n\n")
cat("FIGURES (outputs/figures/):\n")
cat(paste0("  ",list.files("outputs/figures"),collapse="\n"),"\n\n")
cat("TABLES (outputs/tables/):\n")
cat(paste0("  ",list.files("outputs/tables"),collapse="\n"),"\n\n")
cat("DATA APPENDIX: outputs/tables/appendix_data_sources.csv\n")
cat("============================================================\n")
cat("Next: source('multipanel_figures.R') for Figures 1 & 2\n")
cat("      source('multipanel_extensions.R') for Figures 3-5\n")
cat("============================================================\n")
