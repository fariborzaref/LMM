############################################################
# LINEAR MIXED-EFFECTS MODELING (LMM) — DR. FARIBORZ AREF
# CONTEXT: Inequality dynamics across 27 OECD countries, 2002–2021
# PURPOSE: Separate within-country (time) from between-country variation
# MODELS: Random intercepts + random slopes; interaction with health_ineq
# OUTPUTS: Diagnostics, influence analysis, bootstrap CIs, LOCO-CV, plots
############################################################

# = 0) PACKAGES =========================
required <- c(
  "lme4","lmerTest","performance","see","checkmate","data.table",
  "tidyverse","janitor","broom.mixed","patchwork","ggplot2","influence.ME"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(required, library, character.only = TRUE))

set.seed(2025)

# - GLOBAL PLOT STYLE ---------------------
theme_set(theme_minimal(base_size = 13))
options(dplyr.summarise.inform = FALSE)

# = 1) LOAD & VALIDATE DATA ===================
# EXPECTED CSV COLUMNS (case-insensitive):
#   country, year, gini_income, health_ineq, labor_ineq, gdp_pc, unemployment, openness
path <- "oecd_inequality_2002_2021.csv"
assert_file_exists <- function(p) if (!file.exists(p)) {
  stop("File not found: ", p, "\nPlace the OECD CSV next to this script.", call. = FALSE)
}
assert_file_exists(path)

dat <- fread(path) |> janitor::clean_names()

need <- c("country","year","gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "))

# Filter to window and tidy
dat <- dat |>
  filter(year >= 2002, year <= 2021) |>
  mutate(
    country    = as.factor(country),
    year_c     = year - 2002,         # time centered at 2002
    gdp_pc_k   = gdp_pc / 1000,       # scale GDP per capita
    openness01 = openness / 100       # convert 0–1
  )

# Report missingness and simple outlier flags
cat("\n=== DATA HEALTH CHECK ===\n")
na_rate <- sapply(dat[, ..need], function(x) mean(is.na(x)))
print(round(na_rate, 3))

flag_out <- dat |>
  mutate(
    gini_flag  = gini_income  < 15 | gini_income  > 65,
    unemp_flag = unemployment <  1 | unemployment > 30,
    gdp_flag   = gdp_pc       <  5e3 | gdp_pc     > 120e3
  ) |>
  summarise(across(ends_with("flag"), sum, na.rm = TRUE))
cat("\nOutlier flags (count):\n"); print(flag_out)

# Balancedness
panel_count <- dat |> count(country)
if (length(unique(panel_count$n)) != 1) {
  warning("Panel is not balanced; proceeding with REML/ML which tolerate unbalanced data.")
}

# = 2) MODEL SPECIFICATION ===================
# BASELINE: Random intercepts by country
form_base <- gini_income ~ year_c + gdp_pc_k + unemployment + openness01 + (1 | country)

# EXTENDED: Random slopes for gdp_pc_k and unemployment
form_slopes <- gini_income ~ year_c + gdp_pc_k + unemployment + openness01 +
  (gdp_pc_k + unemployment | country)

# INTERACTION: Does health_ineq moderate the GDP→Inequality link?
form_interact <- gini_income ~ year_c + gdp_pc_k * health_ineq + labor_ineq +
  unemployment + openness01 + (1 | country)

# Fit with ML for comparability; switch to REML for final estimates if desired
m1 <- lmer(form_base,   data = dat, REML = FALSE)
m2 <- lmer(form_slopes, data = dat, REML = FALSE)

cat("\n=== MODEL COMPARISON (ML) ===\n")
print(anova(m1, m2))

# Random-effects identifiability (principal components of RE covariance)
cat("\n=== RANDOM-EFFECTS PCA (m2) ===\n")
print(performance::check_random_effects(m2))

# Collinearity check
cat("\n=== COLLINEARITY (FIXED EFFECTS) ===\n")
print(performance::check_collinearity(m2))

# Residual diagnostics
cat("\n=== MODEL DIAGNOSTICS (m2) ===\n")
print(performance::check_model(m2))

cat("\nR2 (Nakagawa):\n"); print(performance::r2_nakagawa(m2))
cat("\nICC:\n");        print(performance::icc(m2))

# = 3) CORE INTERPRETATIONS (m2) ==================
cat("\n=== FIXED EFFECTS (m2) ===\n")
fx <- broom.mixed::tidy(m2, effects = "fixed", conf.int = TRUE)
print(fx |> select(term, estimate, conf.low, conf.high, p.value))

cat("\nINTERPRETATION NOTES:\n")
cat("• year_c > 0  → upward trend in inequality (per year since 2002).\n")
cat("• gdp_pc_k < 0 → higher income levels associated with lower inequality (holding others).\n")
cat("• unemployment > 0 → slack labor markets raise inequality.\n")
cat("• Random slopes: GDP/Unemployment effects vary across countries (heterogeneous responses).\n")

# = 4) INTERACTION MODEL (m3) ==================
m3 <- lmer(form_interact, data = dat, REML = FALSE)
cat("\n=== INTERACTION MODEL (m3: gdp_pc_k × health_ineq) ===\n")
fx3 <- broom.mixed::tidy(m3, effects = "fixed", conf.int = TRUE)
print(fx3 |> select(term, estimate, conf.low, conf.high, p.value))

cat("\nINTERPRETATION (m3):\n")
cat("• gdp_pc_k:health_ineq > 0 → economic gains reduce inequality less where health inequality is high.\n")

# = 5) PARAMETRIC BOOTSTRAP CIs ==================
cat("\n=== PARAMETRIC BOOTSTRAP (m2) — 1000 DRAWS ===\n")
# (Adjust nsim for speed vs. precision)
boot_m2 <- bootMer(
  m2, FUN = function(fit) fixef(fit),
  nsim = 1000, seed = 2025, use.u = TRUE, type = "parametric", verbose = FALSE
)
boot_ci <- t(apply(boot_m2$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)))
print(round(boot_ci, 3))

# = 6) INFLUENCE BY COUNTRY (Δ COEF) ==============
# Leave-one-country-out change in key coefficients (year_c, gdp_pc_k)
cat("\n=== INFLUENCE ANALYSIS: LEAVE-ONE-COUNTRY-OUT (m2) ===\n")
infl <- influence.ME::influence(m2, group = "country")
s_infl <- summary(infl)
print(s_infl)

# Compute Δ for selected terms
base_fx <- fixef(m2)
delta_terms <- c("year_c","gdp_pc_k","unemployment")
deltas <- lapply(delta_terms, function(term) {
  dd <- influence.ME::fixef(infl)[, term] - base_fx[term]
  tibble(country = names(dd), term = term, delta = as.numeric(dd))
}) |> bind_rows()

# = 7) LOCO PREDICTIVE CHECK ==================
# Leave-One-Country-Out CV: refit without a country, predict its trajectory
cat("\n=== LOCO (LEAVE-ONE-COUNTRY-OUT) PREDICTIVE CHECK (m2) ===\n")
loco_countries <- unique(dat$country)
loco_results <- vector("list", length(loco_countries))
names(loco_results) <- as.character(loco_countries)

for (cc in loco_countries) {
  train <- filter(dat, country != cc)
  test  <- filter(dat, country == cc)
  # Refit quickly with ML
  fit_cc <- lmer(form_slopes, data = train, REML = FALSE)
  pred   <- predict(fit_cc, newdata = test, re.form = NULL)
  loco_results[[as.character(cc)]] <- tibble(country = cc, year = test$year, gini_true = test$gini_income, gini_pred = pred)
}
loco_df <- bind_rows(loco_results) |>
  mutate(abs_err = abs(gini_true - gini_pred))

cat("Mean absolute error by country (LOCO):\n")
print(loco_df |>
        group_by(country) |>
        summarise(MAE = mean(abs_err, na.rm = TRUE)) |>
        arrange(desc(MAE)) |>
        head(10))

# = 8) VISUALIZATIONS ======================
# 8a) Predicted inequality trajectories for selected countries
countries_sel <- c("United States","Germany","Japan","Norway","Mexico")
pred_dat <- dat |>
  filter(country %in% countries_sel) |>
  mutate(pred = predict(m2, newdata = ., re.form = NULL))

p_traj <- ggplot(pred_dat, aes(year, pred, color = country)) +
  geom_line(size = 1.1) +
  labs(title = "Predicted Income Inequality Trajectories (LMM, 2002–2021)",
       x = "Year", y = "Predicted Gini")
print(p_traj)

# 8b) Influence plot: Δ coefficient when dropping each country
p_infl <- deltas |>
  ggplot(aes(x = reorder(country, delta), y = delta)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 2) +
  facet_wrap(~ term, scales = "free_y") +
  coord_flip() +
  labs(title = "Leave-One-Country-Out Influence on Key Coefficients",
       x = "Country", y = expression(Delta~Coefficient))
print(p_infl)

# 8c) LOCO error distribution
p_err <- loco_df |>
  ggplot(aes(x = country, y = abs_err)) +
  geom_boxplot(outlier.alpha = 0.2) +
  coord_flip() +
  labs(title = "LOCO Predictive Error by Country",
       x = "Country", y = "Absolute Error (|Gini_true − Gini_pred|)")
print(p_err)

# Save figures
if (!dir.exists("LMM/figs")) dir.create("LMM/figs", recursive = TRUE)
ggsave("LMM/figs/lmm_pred_trajectories.png", p_traj, width = 8, height = 4.6, dpi = 300)
ggsave("LMM/figs/lmm_influence_delta.png",   p_infl, width = 8, height = 5.6, dpi = 300)
ggsave("LMM/figs/lmm_loco_error.png",        p_err, width = 8, height = 6.0, dpi = 300)

# = 9) REPORT TABLE ========================
# (Note: sjPlot::tab_model is convenient; here we produce a tidy table as CSV)
tbl_all <- bind_rows(
  broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE) |> mutate(model = "Random-Intercept"),
  broom.mixed::tidy(m2, effects = "fixed", conf.int = TRUE) |> mutate(model = "Random-Slopes"),
  broom.mixed::tidy(m3, effects = "fixed", conf.int = TRUE) |> mutate(model = "Interaction")
) |> select(model, term, estimate, conf.low, conf.high, p.value)

if (!dir.exists("LMM/out")) dir.create("LMM/out", recursive = TRUE)
fwrite(tbl_all, "LMM/out/lmm_fixed_effects_tables.csv")

# = 10) SAVE MODEL ARTIFACTS ==================
saveRDS(list(
  LMM_base = m1, LMM_slopes = m2, LMM_interact = m3,
  bootstrap_ci = boot_ci, loco_errors = loco_df, influence = deltas
), file = "LMM/lmm_oecd_models.rds")

cat("\n=== DONE ===\nArtifacts saved in LMM/figs and LMM/out; models saved to LMM/lmm_oecd_models.rds\n")
#############################
# END OF SCRIPT
#############################

