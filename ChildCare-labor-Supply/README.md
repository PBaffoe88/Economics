# Childcare Access and Labor Supply

Empirical analysis of how childcare availability affects labor force participation.

# Data Availability

The input datasets are **not distributed** with this repository. They are built from public sources but include processed microdata extracts that are not redistributed here. To reproduce the analysis, construct the seven files below and place them in `data/` (relative to the repository root). All outputs are written to `output/`.

| File | Contents | Dimensions | Source |
|---|---|---|---|
| `pums_stcloud_analysis.csv` | Individual-level microdata, Stearns & Benton Counties, 2009–2024, with labor supply outcomes and demographics | 103,379 × 57 | ACS PUMS (IPUMS USA / Census Bureau) |
| `zip_year_supply_panel.csv` | ZIP × year panel of licensed childcare slots, providers, openings, closures, 2003–2026 | 666 × 10 | MN DHS / DCYF licensed provider records |
| `fred_monthly_wide.csv` | Monthly St. Cloud MSA labor market series (UR, LF size, employment, LFP, nonfarm) | 257 × 10 | FRED (STCL027UR, STCL027LFN, STCL027EMN, LBSSA27, ENUC410640010SA) |
| `qcew_stcloud_6244_panel.csv` | County-year childcare sector employment, establishments, wages (NAICS 6244) | 29 × 26 | MN DEED QCEW |
| `deed_oews_master.csv` | Childcare worker wage distributions by region (SOC 39-9011) | 23 × 15 | MN DEED OEWS, Q1 2026 |
| `jvs_centralMN_naics62_panel.csv` | Job vacancies, NAICS 62, Central Minnesota | 24 × 10 | MN DEED Job Vacancy Survey |
| `mpls_fed_surveys_complete.csv` | Provider financial stress indicators (insurance costs, tuition changes) | 6 × 20 | Federal Reserve Bank of Minneapolis Annual Childcare Provider Surveys, 2021–2026 |

Key derived variables expected in the PUMS file: `in_lf`, `employed`, `fulltime`, `parent_u5`, `female`, `single_parent`, `low_income`, `immigrant`, `college`, `AGE`, `POVERTY` (income as % of FPL), `log_wage`, `county_name`, `YEAR`.

Questions about data construction are welcome — open an issue on this repository.
