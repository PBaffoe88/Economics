# Crude Oil Momentum & Geopolitical Risk

**A full econometric pipeline linking time-series momentum in crude oil to geopolitical risk shocks.**

This project asks whether geopolitical risk drives crude oil price dynamics — in returns, in volatility, and in trend regimes. It combines the Caldara–Iacoviello (2022, *AER*) daily Geopolitical Risk (GPR) index with WTI and Brent spot prices, and works through momentum construction, volatility modeling, impulse responses, and regime switching in a single reproducible R script.

## Key Features

- **Fully live and reproducible** — no local data files and no API keys. Oil prices come from FRED via `quantmod` and the daily GPR index is downloaded from the authors' website at runtime. Clone, `source()`, done.
- **Time-series momentum** à la Moskowitz, Ooi & Pedersen (2012): 1/3/6/12-month lookback windows, a 12-month TSMOM signal with its annualized Sharpe ratio on WTI, plus classical 50/200-day golden-cross trend regimes.
- **Volatility channel** — EGARCH(1,1) with Student-t errors and the GPR shock as an external regressor in the variance equation, testing whether geopolitical risk raises oil variance beyond what past returns explain.
- **Dynamic causal responses** — Jordà (2005) local projections of cumulative WTI returns to a GPR shock over a 20-trading-day horizon with Newey–West standard errors.
- **Predictive regressions** — does today's GPR shock forecast next-month returns, and does the answer depend on the prevailing trend regime? Estimated with `fixest` and Newey–West VCOV.
- **Regime switching** — a 2-state Markov-switching model (monthly) separating calm and crisis regimes, with smoothed crisis-regime probabilities plotted over time.

## Methodological Details Worth Noting

- The GPR shock is defined as the log index's deviation from its own 252-day rolling mean, so it captures abnormal geopolitical tension rather than the level of the index.
- Non-positive prices are excluded with documentation (WTI settled at −$37.63 on 2020-04-20; log returns are undefined there).
- The GPR spreadsheet's column names occasionally change upstream; the loader detects the date and headline `GPRD` columns programmatically rather than hard-coding positions.
- Volatility is plotted on a log scale because the April 2020 episode pushes annualized volatility above 1000%, which flattens every other episode on a linear axis.

## How to Run

```r
source("oil_momentum_gpr_pipeline.R")
```

The script installs any missing packages, pulls all data live, prints diagnostics (ADF, Lo–MacKinlay variance ratios, EGARCH coefficients, predictive regressions, regime-switching estimates) to the console, and saves six publication-quality figures to `figures/`:

`irf_gpr_oil.png` · `wti_trend_regimes.png` · `wti_momentum_windows.png` · `gpr_vs_wti.png` · `egarch_volatility_gpr.png` · `regime_probabilities.png`

## Requirements

R ≥ 4.2 with: quantmod, readxl, httr, dplyr, tidyr, lubridate, zoo, sandwich, lmtest, rugarch, vrtest, tseries, lpirfs, ggplot2, fixest, MSwM, purrr (installed automatically if missing).

## Author

Paul Baffoe — Economist / Econometrician

## References

- Caldara, D., & Iacoviello, M. (2022). Measuring Geopolitical Risk. *American Economic Review*, 112(4).
- Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local Projections. *American Economic Review*, 95(1).
- Moskowitz, T., Ooi, Y. H., & Pedersen, L. H. (2012). Time Series Momentum. *Journal of Financial Economics*, 104(2).
- Nelson, D. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. *Econometrica*, 59(2).
