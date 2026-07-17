# ============================================================================
# Crude Oil Momentum & Geopolitical Risk: A Full Econometric Pipeline
# ----------------------------------------------------------------------------
# Data (all public, pulled live at runtime):
#   - WTI & Brent spot prices: FRED (DCOILWTICO, DCOILBRENTEU) - no API key needed
#   - Geopolitical Risk (GPR) index: Caldara & Iacoviello (2022, AER),
#     daily series from matteoiacoviello.com
#
# Methodology:
#   1. Time-series momentum construction (Moskowitz, Ooi & Pedersen 2012)
#   2. Trend diagnostics: ADF, variance ratio tests
#   3. Volatility: GARCH(1,1) and EGARCH with GPR in the variance equation
#   4. Local projections (Jorda 2005): dynamic response of oil returns
#      to geopolitical risk shocks, Newey-West SEs
#   5. Predictive regressions: does GPR forecast oil momentum?
#   6. Markov regime-switching trend model (optional, monthly)
# ============================================================================

# ---- 0. Packages -----------------------------------------------------------
pkgs <- c("quantmod", "readxl", "httr", "dplyr", "tidyr", "lubridate",
          "zoo", "sandwich", "lmtest", "rugarch", "vrtest", "tseries",
          "lpirfs", "ggplot2", "fixest", "MSwM", "purrr")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. Live data ingestion ------------------------------------------------

## 1a. Oil prices from FRED (no key required via quantmod)
getSymbols(c("DCOILWTICO", "DCOILBRENTEU"), src = "FRED")

oil <- merge(DCOILWTICO, DCOILBRENTEU) |>
  as.data.frame() |>
  tibble::rownames_to_column("date") |>
  rename(wti = DCOILWTICO, brent = DCOILBRENTEU) |>
  mutate(date = as.Date(date)) |>
  filter(!is.na(wti))

## 1b. Daily Geopolitical Risk index (Caldara & Iacoviello)
gpr_url <- "https://www.matteoiacoviello.com/gpr_files/data_gpr_daily_recent.xls"
tmp <- tempfile(fileext = ".xls")
GET(gpr_url, write_disk(tmp, overwrite = TRUE), timeout(120))
gpr_raw <- read_excel(tmp)

# Column names occasionally change; grab date + headline GPRD robustly
date_col <- names(gpr_raw)[grepl("date", names(gpr_raw), ignore.case = TRUE)][1]
gprd_col <- names(gpr_raw)[grepl("^GPRD$", names(gpr_raw))][1]
gpr <- gpr_raw |>
  transmute(date = as.Date(.data[[date_col]]),
            gprd = as.numeric(.data[[gprd_col]])) |>
  filter(!is.na(gprd))

## 1c. Merge to a daily analysis panel
df <- oil |>
  inner_join(gpr, by = "date") |>
  arrange(date) |>
  # Drop non-positive prices (WTI settled at -$37.63 on 2020-04-20);
  # log returns are undefined there. Documented sample exclusion.
  filter(wti > 0, is.na(brent) | brent > 0) |>
  mutate(
    ret_wti   = 100 * (log(wti) - log(lag(wti))),
    ret_brent = 100 * (log(brent) - log(lag(brent))),
    # Standardized GPR shock: log deviation from a 1-year rolling mean
    log_gpr   = log(gprd + 1),
    gpr_shock = log_gpr - rollmeanr(log_gpr, k = 252, fill = NA)
  ) |>
  filter(!is.na(ret_wti), !is.na(gpr_shock))

# ---- 2. Momentum construction ----------------------------------------------
# Time-series momentum (TSMOM): sign of past k-day cumulative return
# Also classical technical measures for robustness: MA crossover, RSI-style.

mom_windows <- c(21, 63, 126, 252)  # ~1, 3, 6, 12 months of trading days

df <- df |>
  mutate(
    across(all_of("wti"), \(x) x, .names = "px"),
    !!!setNames(
      map(mom_windows, \(k) rlang::expr(100 * (log(px) - log(lag(px, !!k))))),
      paste0("mom_", mom_windows)
    ),
    ma50  = rollmeanr(px, 50,  fill = NA),
    ma200 = rollmeanr(px, 200, fill = NA),
    golden = as.integer(ma50 > ma200),                # trend regime dummy
    tsmom_sig_252 = sign(mom_252)                     # MOP (2012) signal
  )

# Realized momentum "strategy" return: hold sign of 12m momentum, next-day ret
df <- df |> mutate(tsmom_ret = lag(tsmom_sig_252) * ret_wti)

cat("\n--- TSMOM (12m) annualized Sharpe on WTI ---\n")
sr_vec <- na.omit(df$tsmom_ret)
sr <- mean(sr_vec) / sd(sr_vec) * sqrt(252)
print(round(sr, 3))

# ---- 3. Trend & efficiency diagnostics --------------------------------------
cat("\n--- ADF test on log WTI (levels: expect unit root) ---\n")
print(adf.test(na.omit(log(df$px))))

cat("\n--- Lo-MacKinlay variance ratio (returns: VR>1 => momentum) ---\n")
print(Lo.Mac(na.omit(df$ret_wti), kvec = c(2, 5, 10, 20)))

# ---- 4. Volatility: does geopolitical risk drive oil variance? -------------
# EGARCH(1,1) with GPR shock as external regressor in the variance equation.
sub <- na.omit(df[, c("date", "ret_wti", "gpr_shock")])

spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1),
                        external.regressors = as.matrix(sub$gpr_shock)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)
fit_garch <- ugarchfit(spec, data = sub$ret_wti)
cat("\n--- EGARCH(1,1)-t with GPR in variance equation ---\n")
print(fit_garch@fit$matcoef)   # vxreg1 = GPR effect on log-variance

# ---- 5. Local projections: oil response to a GPR shock ----------------------
# Jorda (2005): ret_{t+h} = a_h + b_h * gpr_shock_t + controls + e_{t+h}
# Cumulative IRF over h = 0..20 trading days, Newey-West SEs.

lp_data <- sub |>
  mutate(cum0 = ret_wti) |>
  as.data.frame()

horizons <- 0:20
irf <- map_dfr(horizons, function(h) {
  y <- with(lp_data, rollapply(ret_wti, width = h + 1, sum,
                               align = "left", fill = NA))
  d <- data.frame(y = y, x = lp_data$gpr_shock,
                  l1 = dplyr::lag(lp_data$ret_wti, 1),
                  l2 = dplyr::lag(lp_data$ret_wti, 2),
                  lg = dplyr::lag(lp_data$gpr_shock, 1)) |> na.omit()
  m <- lm(y ~ x + l1 + l2 + lg, data = d)
  v <- NeweyWest(m, lag = h + 1, prewhite = FALSE)
  data.frame(h = h,
             b  = coef(m)["x"],
             se = sqrt(diag(v))["x"])
})

irf <- irf |> mutate(lo = b - 1.96 * se, hi = b + 1.96 * se)

# ---- Shared colorful theme & output folder ---------------------------------
if (!dir.exists("figures")) dir.create("figures")

pal <- c(oil = "#D1495B", gpr = "#00798C", calm = "#3CB371",
         crisis = "#EDAE49", accent = "#6A4C93")

theme_oil <- theme_minimal(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", color = "#1B2A41", size = 15),
    plot.subtitle   = element_text(color = "#5C677D"),
    plot.background = element_rect(fill = "#FDFBF7", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E4DFDA"),
    legend.position = "top",
    strip.text      = element_text(face = "bold", color = "#1B2A41")
  )

p_irf <- ggplot(irf, aes(h, b)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = pal["gpr"], alpha = .25) +
  geom_line(linewidth = 1.3, color = pal["oil"]) +
  geom_point(color = pal["oil"], size = 2) +
  geom_hline(yintercept = 0, linetype = 2, color = "#5C677D") +
  labs(title = "Cumulative WTI return response to a GPR shock",
       subtitle = "Local projections (Jorda), Newey-West 95% bands",
       x = "Trading days after shock", y = "Cumulative log return x100") +
  theme_oil
ggsave("figures/irf_gpr_oil.png", p_irf, width = 8, height = 5, dpi = 300)

# ---- 6. Does GPR predict momentum? ------------------------------------------
# Predictive regression of forward k-day returns on current GPR shock,
# conditioning on the prevailing trend regime (golden cross).

pred <- df |>
  mutate(fwd_21 = lead(mom_21, 21)) |>          # next-month return
  select(date, fwd_21, gpr_shock, golden, mom_252) |>
  na.omit()

m_pred <- feols(fwd_21 ~ gpr_shock * golden + mom_252, data = pred,
                vcov = NW(lag = 21) ~ date)
cat("\n--- Predictive regression: 1m-ahead returns on GPR x trend regime ---\n")
print(summary(m_pred))

# ---- 7. (Optional) Markov-switching trend model, monthly -------------------
mth <- df |>
  mutate(ym = floor_date(date, "month")) |>
  group_by(ym) |>
  summarise(ret = sum(ret_wti, na.rm = TRUE),
            gpr = mean(gprd, na.rm = TRUE)) |>
  na.omit()

base <- lm(ret ~ gpr, data = mth)
ms <- tryCatch(
  msmFit(base, k = 2, sw = c(TRUE, TRUE, TRUE)),
  error = function(e) NULL
)
if (!is.null(ms)) {
  cat("\n--- 2-state Markov-switching: GPR effect by regime ---\n")
  print(summary(ms))
}

# ---- 8. Publication figures (all saved as PNG in ./figures) -----------------

## 8a. WTI price with 50/200-day MAs, shaded by trend regime (golden cross)
# Collapse consecutive same-regime days into contiguous blocks (avoids
# the vertical-striping artifact from per-day geom_tile shading)
trend_df <- df |> filter(!is.na(ma200))
reg_blocks <- trend_df |>
  mutate(block = cumsum(golden != lag(golden, default = dplyr::first(golden)))) |>
  group_by(block, golden) |>
  summarise(xmin = min(date), xmax = max(date), .groups = "drop")

p_trend <- ggplot(trend_df, aes(date)) +
  geom_rect(data = reg_blocks, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                fill = factor(golden)), alpha = .18) +
  geom_line(aes(y = px,   color = "WTI spot"),   linewidth = .4) +
  geom_line(aes(y = ma50, color = "50-day MA"),  linewidth = .8) +
  geom_line(aes(y = ma200, color = "200-day MA"), linewidth = .8) +
  scale_color_manual(NULL, values = c("WTI spot" = "#5C677D",
                                      "50-day MA" = pal[["oil"]],
                                      "200-day MA" = pal[["gpr"]])) +
  scale_fill_manual(NULL, values = c("0" = pal[["crisis"]], "1" = pal[["calm"]]),
                    labels = c("0" = "Downtrend", "1" = "Uptrend"),
                    guide = guide_legend(override.aes = list(alpha = .4))) +
  labs(title = "WTI crude: price, moving averages, and trend regimes",
       subtitle = "Shading = 50/200-day golden-cross regime",
       x = NULL, y = "USD per barrel") +
  theme_oil
ggsave("figures/wti_trend_regimes.png", p_trend, width = 10, height = 5.5, dpi = 300)

## 8b. Momentum across horizons (1m / 3m / 6m / 12m)
p_mom <- df |>
  select(date, starts_with("mom_")) |>
  pivot_longer(-date, names_to = "window", values_to = "mom") |>
  mutate(window = factor(window, levels = paste0("mom_", mom_windows),
                         labels = c("1 month", "3 months", "6 months", "12 months"))) |>
  filter(!is.na(mom)) |>
  ggplot(aes(date, mom, color = window)) +
  geom_hline(yintercept = 0, linetype = 2, color = "#5C677D") +
  geom_line(linewidth = .5, show.legend = FALSE) +
  facet_wrap(~window, ncol = 2, scales = "free_y") +
  scale_color_manual(values = unname(pal[1:4])) +
  labs(title = "WTI time-series momentum by lookback window",
       subtitle = "Trailing cumulative log returns (x100)",
       x = NULL, y = "Momentum (%)") +
  theme_oil
ggsave("figures/wti_momentum_windows.png", p_mom, width = 10, height = 6.5, dpi = 300)

## 8c. GPR index vs WTI (dual panel, avoids dual-axis distortion)
p_gpr <- df |>
  select(date, `WTI (USD/bbl)` = px, `Geopolitical Risk index` = gprd) |>
  pivot_longer(-date) |>
  ggplot(aes(date, value, color = name)) +
  geom_line(linewidth = .5, show.legend = FALSE) +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(`Geopolitical Risk index` = pal[["gpr"]],
                                `WTI (USD/bbl)` = pal[["oil"]])) +
  labs(title = "Crude oil and geopolitical risk",
       subtitle = "Daily GPR index (Caldara-Iacoviello) and WTI spot",
       x = NULL, y = NULL) +
  theme_oil
ggsave("figures/gpr_vs_wti.png", p_gpr, width = 10, height = 6.5, dpi = 300)

## 8d. EGARCH conditional volatility, colored by GPR shock intensity
# Log scale: the April 2020 episode pushes annualized vol above 1000%,
# which flattens every other episode on a linear axis.
vol_df <- sub |>
  mutate(sigma = as.numeric(sigma(fit_garch)) * sqrt(252))  # annualized %
p_vol <- ggplot(vol_df, aes(date, sigma, color = gpr_shock)) +
  geom_line(linewidth = .8) +
  scale_y_log10(breaks = c(10, 25, 50, 100, 250, 500, 1000, 2000),
                labels = scales::comma) +
  scale_color_gradient2(name = "GPR shock",
                        low = pal[["calm"]], mid = "#BFC9CA",
                        high = pal[["oil"]], midpoint = 0) +
  labs(title = "EGARCH conditional volatility of WTI returns",
       subtitle = "Log scale; line colored by contemporaneous geopolitical-risk shock",
       x = NULL, y = "Annualized volatility (%, log scale)") +
  theme_oil
ggsave("figures/egarch_volatility_gpr.png", p_vol, width = 10, height = 5, dpi = 300)

## 8e. Markov-switching smoothed regime probabilities (crisis regime)
if (!is.null(ms)) {
  reg_df <- mth |>
    mutate(p_crisis = ms@Fit@smoProb[-1, which.max(ms@std)]) # high-vol regime
  p_reg <- ggplot(reg_df, aes(ym)) +
    geom_area(aes(y = p_crisis), fill = pal[["crisis"]], alpha = .55) +
    geom_line(aes(y = p_crisis), color = pal[["oil"]], linewidth = .6) +
    labs(title = "Probability of the high-volatility (crisis) regime",
         subtitle = "Smoothed probabilities, 2-state Markov-switching model",
         x = NULL, y = "P(crisis regime)") +
    theme_oil
  ggsave("figures/regime_probabilities.png", p_reg, width = 10, height = 4.5, dpi = 300)
}

cat("\nDone. Figures saved to ./figures:\n",
    " - irf_gpr_oil.png\n",
    " - wti_trend_regimes.png\n",
    " - wti_momentum_windows.png\n",
    " - gpr_vs_wti.png\n",
    " - egarch_volatility_gpr.png\n",
    " - regime_probabilities.png\n")
