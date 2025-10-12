############################################################
# Linear Mixed-Effects Modeling (LMM) – Dr. Fariborz Aref
# Context: Inequality dynamics across 27 OECD countries, 2002–2021
# Goal: Model within-country (temporal) and between-country variation
# Method: Random intercepts & random slopes for economic covariates
############################################################

# ---- 0) Packages ----
required <- c("lme4", "lmerTest", "sjPlot", "performance",
              "tidyverse", "data.table")
to_install <- setdiff(required, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(required, library, character.only = TRUE)

# ---- 1) Load & Prepare Data ----
# Expected CSV: oecd_inequality_2002_2021.csv
# Variables: country, year, gini_income, health_ineq, labor_ineq,
#            gdp_pc, unemployment, openness
dat <- fread("oecd_inequality_2002_2021.csv") |> janitor::clean_names()

dat <- dat |>
  filter(year >= 2002, year <= 2021) |>
  mutate(
    country = as.factor(country),
    year_c = year - 2002,                   # centered time
    gdp_pc_k = gdp_pc / 1000,               # scale GDP per capita
    openness = openness / 100               # convert to 0-1
  )

# ---- 2) Model Specification ----
# Baseline model: random intercepts by country
m1 <- lmer(gini_income ~ year_c + gdp_pc_k + unemployment + openness +
             (1 | country), data = dat, REML = FALSE)

# Extended model: random slopes for GDP growth & unemployment
m2 <- lmer(gini_income ~ year_c + gdp_pc_k + unemployment + openness +
             (gdp_pc_k + unemployment | country), data = dat, REML = FALSE)

# Model comparison
anova(m1, m2)

# ---- 3) Model Diagnostics ----
check_model(m2)       # visual diagnostics
r2_nakagawa(m2)       # marginal & conditional R²
icc(m2)               # intra-class correlation

# ---- 4) Interpretation ----
summary(m2)

# Extract fixed effects table
coefs <- summary(m2)$coefficients |> as.data.frame()
coefs$Variable <- rownames(coefs)
rownames(coefs) <- NULL
coefs <- coefs |> relocate(Variable)
print(coefs)

# Random-effects structure
re_var <- as.data.frame(VarCorr(m2))
print(re_var)

cat("\nInterpretation Guide:\n")
cat("• Significant positive coefficient on 'year_c' → rising inequality over time.\n")
cat("• GDP_pc_k negative → richer countries show lower Gini, controlling for others.\n")
cat("• Random slopes indicate that GDP-inequality relationship varies by country.\n")
cat("• ICC quantifies share of total variance due to between-country differences.\n")

# ---- 5) Visualization ----
# Predicted inequality trajectories for selected countries
countries_sel <- c("United States", "Germany", "Japan", "Norway", "Mexico")

pred_dat <- dat |>
  filter(country %in% countries_sel) |>
  mutate(pred = predict(m2, newdata = ., re.form = NULL))

ggplot(pred_dat, aes(x = year, y = pred, color = country)) +
  geom_line(size = 1.1) +
  labs(title = "Predicted Income Inequality Trajectories (LMM, 2002–2021)",
       y = "Predicted Gini", x = "Year") +
  theme_minimal(base_size = 13)

# ---- 6) Cross-Domain Interaction Example ----
# Test whether health inequality moderates effect of GDP on income inequality
m3 <- lmer(gini_income ~ year_c + gdp_pc_k * health_ineq + labor_ineq +
             unemployment + (1 | country), data = dat, REML = FALSE)
summary(m3)

cat("\nInterpretation:\n")
cat("• The interaction term gdp_pc_k:health_ineq tests whether economic growth\n")
cat("  reduces inequality less in countries with higher health inequality.\n")
cat("• Significant positive coefficient → social disparities dilute economic gains.\n")

# ---- 7) Reporting Table ----
tab_model(m1, m2, m3, show.re.var = TRUE,
          title = "Linear Mixed-Effects Models for OECD Income Inequality (2002–2021)",
          dv.labels = c("Random-Intercept Model", "Random-Slope Model", "Interaction Model"))

# ---- 8) Save Artifacts ----
saveRDS(list(LMM_base = m1, LMM_slope = m2, LMM_int = m3),
        file = "LMM/lmm_oecd_models.rds")

############################################################
# End of Script
############################################################
