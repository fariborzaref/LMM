############################################################
# LINEAR MIXED EFFECTS MODELING (LMM)
# Author: Dr. Fariborz Aref
# Context: Inequality dynamics across 27 OECD countries, 2002 to 2021
# Purpose: Separate within country time from between country variation
# Models: Random intercepts and random slopes, test interaction with health_ineq
# Outputs: Diagnostics, influence analysis, bootstrap CIs, LOCO CV, plots
############################################################

# 0) Packages
required <- c(
  "lme4","lmerTest","performance","see","checkmate","data.table",
  "tidyverse","janitor","broom.mixed","patchwork","ggplot2","influence.ME"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))

set.seed(2025)
options(dplyr.summarise.inform = FALSE, stringsAsFactors = FALSE)

# Global plot style, compact academic typography
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 11, base_family = "serif") +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title   = ggplot2::element_text(size = 10),
      axis.text    = ggplot2::element_text(size = 9),
      strip.text   = ggplot2::element_text(size = 10, face = "bold"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
)

# 1) Load and validate data
# Expected columns, case insensitive:
# country, year, gini_income, health_ineq, labor_ineq, gdp_pc, unemployment, openness
path <- "oecd_inequality_2002_2021.csv"
if (!file.exists(path)) stop("File not found: ", path, "\nPlace the OECD CSV next to this script.", call. = FALSE)

dat <- data.table::fread(path) |> janitor::clean_names()

need <- c("country","year","gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "))

# Filter window and tidy
dat <- dat |>
  dplyr::filter(year >= 2002, year <= 2021) |>
  dplyr::mutate(
    country    = as.factor(country),
    year_c     = year - 2002,      # time centered at 2002
    gdp_pc_k   = gdp_pc / 1000,    # scale GDP per capita
    openness01 = openness / 100    # convert 0 to 1
  )

# Data health check
cat("\n=== DATA HEALTH CHECK ===\n")
na_rate <- sapply(dat[, need, with = FALSE], function(x) mean(is.na(x)))
print(round(na_rate, 3))

flag_out <- dat |>
  dplyr::mutate(
    gini_flag  = gini_income  < 15 | gini_income  > 65,
    unemp_flag = unemployment <  1 | unemployment > 30,
    gdp_flag   = gdp_pc       <  5000 | gdp_pc    > 120000
  ) |>
  dplyr::summarise(dplyr::across(dplyr::ends_with("flag"), \(x) sum(x, na.rm = TRUE)))
cat("\nOutlier flags count:\n"); print(flag_out)

panel_count <- dat |> dplyr::count(country)
if (length(unique(panel_count$n)) != 1) {
  warning("Panel is not balanced. Proceed with ML or REML which tolerate unbalanced data.")
}

# 2) Model specification
form_base <- gini_income ~ year_c + gdp_pc_k + unemployment + openness01 + (1 | country)
form_slopes <- gini_income ~ year_c + gdp_pc_k + unemployment + openness01 + (gdp_pc_k + unemployment | country)
form_interact <- gini_income ~ year_c + gdp_pc_k * health_ineq + labor_ineq + unemployment + openness01 + (1 | country)

# Fit with ML for comparability, switch to REML for final estimates if desired
m1 <- lme4::lmer(form_base,   data = dat, REML = FALSE)
m2 <- lme4::lmer(form_slopes, data = dat, REML = FALSE)

cat("\n=== MODEL COMPARISON ML ===\n")
print(anova(m1, m2))

# Random effects identifiability
cat("\n=== RANDOM EFFECTS PCA m2 ===\n")
print(performance::check_random_effects(m2))

# Collinearity
cat("\n=== COLLINEARITY FIXED EFFECTS m2 ===\n")
print(performance::check_collinearity(m2))

# Residual diagnostics
cat("\n=== MODEL DIAGNOSTICS m2 ===\n")
print(performance::check_model(m2))

cat("\nR2 Nakagawa:\n"); print(performance::r2_nakagawa(m2))
cat("\nICC:\n");        print(performance::icc(m2))

# 3) Core interpretations
cat("\n=== FIXED EFFECTS m2 ===\n")
fx <- broom.mixed::tidy(m2, effects = "fixed", conf.int = TRUE)
print(fx |> dplyr::select(term, estimate, conf.low, conf.high, p.value))

cat("\nINTERPRETATION NOTES\n")
cat("year_c positive suggests upward trend in inequality per year since 2002.\n")
cat("gdp_pc_k negative suggests higher income levels relate to lower inequality, all else equal.\n")
cat("unemployment positive suggests slack labor markets raise inequality.\n")
cat("random slopes indicate GDP and unemployment effects vary by country.\n")

# 4) Interaction model
m3 <- lme4::lmer(form_interact, data = dat, REML = FALSE)
cat("\n=== INTERACTION MODEL m3 gdp_pc_k x health_ineq ===\n")
fx3 <- broom.mixed::tidy(m3, effects = "fixed", conf.int = TRUE)
print(fx3 |> dplyr::select(term, estimate, conf.low, conf.high, p.value))

cat("\nINTERPRETATION m3\n")
cat("gdp_pc_k:health_ineq positive suggests economic gains reduce inequality less where health inequality is high.\n")

# 5) Parametric bootstrap CIs
cat("\n=== PARAMETRIC BOOTSTRAP m2 1000 DRAWS ===\n")
boot_m2 <- lme4::bootMer(
  m2, FUN = function(fit) lme4::fixef(fit),
  nsim = 1000, seed = 2025, use.u = TRUE, type = "parametric", verbose = FALSE
)
boot_ci <- t(apply(boot_m2$t, 2, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)))
print(round(boot_ci, 3))

# 6) Influence by country
cat("\n=== INFLUENCE ANALYSIS LEAVE ONE COUNTRY OUT m2 ===\n")
infl <- influence.ME::influence(m2, group = "country")
print(summary(infl))

base_fx <- lme4::fixef(m2)
delta_terms <- c("year_c","gdp_pc_k","unemployment")
deltas <- lapply(delta_terms, function(term) {
  dd <- influence.ME::fixef(infl)[, term] - base_fx[term]
  tibble::tibble(country = names(dd), term = term, delta = as.numeric(dd))
}) |> dplyr::bind_rows()

# 7) LOCO predictive check
cat("\n=== LOCO PREDICTIVE CHECK m2 ===\n")
loco_countries <- unique(dat$country)
loco_results <- vector("list", length(loco_countries))
names(loco_results) <- as.character(loco_countries)

for (cc in loco_countries) {
  train <- dplyr::filter(dat, country != cc)
  test  <- dplyr::filter(dat, country == cc)
  fit_cc <- lme4::lmer(form_slopes, data = train, REML = FALSE)
  pred   <- stats::predict(fit_cc, newdata = test, re.form = NULL)
  loco_results[[as.character(cc)]] <- tibble::tibble(
    country = cc, year = test$year, gini_true = test$gini_income, gini_pred = pred
  )
}
loco_df <- dplyr::bind_rows(loco_results) |>
  dplyr::mutate(abs_err = abs(gini_true - gini_pred))

cat("Mean absolute error by country LOCO:\n")
print(loco_df |>
        dplyr::group_by(country) |>
        dplyr::summarise(MAE = mean(abs_err, na.rm = TRUE)) |>
        dplyr::arrange(dplyr::desc(MAE)) |>
        dplyr::slice_head(n = 10))

# 8) Visualizations
countries_sel <- c("United States","Germany","Japan","Norway","Mexico")
pred_dat <- dat |>
  dplyr::filter(country %in% countries_sel) |>
  dplyr::mutate(pred = stats::predict(m2, newdata = ., re.form = NULL))

p_traj <- ggplot2::ggplot(pred_dat, ggplot2::aes(year, pred, color = country)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(title = "Predicted income inequality trajectories, 2002 to 2021",
                x = "Year", y = "Predicted Gini")

p_infl <- deltas |>
  ggplot2::ggplot(ggplot2::aes(x = reorder(country, delta), y = delta)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 3) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~ term, scales = "free_y") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "Leave one country out influence on key coefficients",
                x = "Country", y = expression(Delta~Coefficient))

p_err <- loco_df |>
  ggplot2::ggplot(ggplot2::aes(x = country, y = abs_err)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "LOCO predictive error by country",
                x = "Country", y = "Absolute error")

# Save outputs
if (!dir.exists("LMM/figs")) dir.create("LMM/figs", recursive = TRUE)
if (!dir.exists("LMM/out"))  dir.create("LMM/out",  recursive = TRUE)

ggplot2::ggsave("LMM/figs/lmm_pred_trajectories.png", p_traj, width = 7.5, height = 4.4, dpi = 300)
ggplot2::ggsave("LMM/figs/lmm_influence_delta.png",   p_infl, width = 7.5, height = 5.2, dpi = 300)
ggplot2::ggsave("LMM/figs/lmm_loco_error.png",        p_err,  width = 7.5, height = 5.6, dpi = 300)

tbl_all <- dplyr::bind_rows(
  broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE) |> dplyr::mutate(model = "Random intercept"),
  broom.mixed::tidy(m2, effects = "fixed", conf.int = TRUE) |> dplyr::mutate(model = "Random slopes"),
  broom.mixed::tidy(m3, effects = "fixed", conf.int = TRUE) |> dplyr::mutate(model = "Interaction")
) |> dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

if (!dir.exists("LMM/out")) dir.create("LMM/out", recursive = TRUE)
data.table::fwrite(tbl_all, "LMM/out/lmm_fixed_effects_tables.csv")

saveRDS(list(
  LMM_base = m1, LMM_slopes = m2, LMM_interact = m3,
  bootstrap_ci = boot_ci, loco_errors = loco_df, influence = deltas
), file = "LMM/lmm_oecd_models.rds")

cat("\n=== DONE ===\nArtifacts saved in LMM/figs and LMM/out. Models saved to LMM/lmm_oecd_models.rds\n")
