# Buffer Depletion and Price Instability

**How Monetary Tightening Amplifies Producer Price Volatility Through the Inventory Channel**

An empirical study of the inventory channel of monetary policy transmission in U.S. manufacturing. Using a monthly panel of eight NAICS-3 manufacturing industries (1992–present), the analysis shows that contractionary monetary policy shocks drain inventory buffers and that industries with structurally higher inventory intensity experience larger increases in producer price volatility after tightening.

## Key Features

- **Fully reproducible pipeline** — a single self-contained R script downloads all data from official sources, cleans it, estimates all models, and exports publication-ready figures and tables. No local data files required.
- **Identification** — Bartik-style interaction of high-frequency monetary policy surprises (Bauer & Swanson 2023) with pre-determined (1992–2006) industry inventory intensity, estimated via Jordà local projections with industry and date fixed effects and two-way clustered standard errors.
- **Robustness** — placebo shocks, split samples, crisis-period exclusion, alternative shock measures (first-differenced FFR), alternative volatility windows, wild bootstrap peak confidence intervals.
- **Extensions** — volatility decomposition, tightening/easing asymmetry, cross-industry determinants (export exposure, concentration, upstreamness, perishability), post-pandemic structural break, and an international OECD/BIS scaffold.

## Data Sources

| Data | Source | Access |
|---|---|---|
| Inventories, shipments, inventory-to-shipments ratios | Census M3 Survey | FRED API |
| Producer price indices (NAICS-3) | BLS PPI | FRED API |
| Monetary policy surprises | Bauer & Swanson (2023) | SF Fed website |
| Macro controls (FFR, Treasuries, WTI, CPI, INDPRO) | FRED | FRED API |
| Recession dates | NBER Business Cycle Dating Committee | Hard-coded official dates |
| Export shares | Census ASM 2021 published tables | Hard-coded with citation* |
| Concentration (HHI) | 2017 Economic Census | Census API, published-table fallback |
| Upstreamness | BEA I-O 2022 / Antràs & Chor (2013) | BEA API, replication-data fallback |
| Perishability | USDA ERS / FDA 21 CFR | Categorical, cited |

\* Export shipments are not available through the Census ASM API (verified against the API's dataset and variable listings), so values come from the published ASM tables.

## How to Run

```r
# Optional (raises API rate limits; script runs fine without them):
# Sys.setenv(CENSUS_API_KEY = "your_key")   # api.census.gov/data/key_signup.html
# Sys.setenv(BEA_API_KEY    = "your_key")   # apps.bea.gov/API/signup/

source("buffer_depletion_analysis.R")
```

The script installs any missing packages, creates `data_raw/`, `data_clean/`, and `outputs/` directories relative to the working directory, and writes all figures to `outputs/figures/` and all tables to `outputs/tables/`. A full data-source appendix is exported to `outputs/tables/appendix_data_sources.csv`.

All external API calls (FRED, Census, BEA, OECD, BIS) include retry logic, response validation, and documented fallbacks, so a temporary outage of any single source will not stop the run.

## Requirements

R ≥ 4.2 with: tidyverse, lubridate, zoo, slider, readxl, janitor, fixest, broom, modelsummary, scales, httr, jsonlite (installed automatically if missing). The optional international extension additionally uses the OECD package.

## Selected Outputs

- Local projection IRFs: inventory-to-shipments ratio, inventory growth, PPI inflation, and PPI volatility responses to monetary shocks
- Placebo and split-sample robustness figures
- Tightening vs. easing asymmetry IRFs with Wald tests
- Cross-industry peak-response regression on verified industry characteristics
- Pre/post-2020 structural break and rolling-window estimates

## Author

Paul Baffoe — Economist / Econometrician

## References

- Antràs, P., & Chor, D. (2013). Organizing the Global Value Chain. *Econometrica*.
- Bauer, M. D., & Swanson, E. T. (2023). A Reassessment of Monetary Policy Surprises and High-Frequency Identification. *NBER Macroeconomics Annual*.
- Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local Projections. *American Economic Review*.
- Mian, A., & Sufi, A. (2014). What Explains the 2007–2009 Drop in Employment? *Econometrica*.
