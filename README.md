# Linear Mixed Effects Modeling, OECD Inequality

**Author:** Dr. Fariborz Aref  
**Discipline:** Quantitative Sociology and Inequality Research  
**License:** MIT  

### Purpose
Model income inequality with country level random effects to separate within country change over time from between country differences. Includes model comparison, diagnostics, influence, bootstrap confidence intervals, and leave one country out prediction.

### Structure

### Key Methods
- Random intercepts and random slopes by country  
- Interaction of GDP per capita with health inequality  
- Diagnostics with `performance` package  
- Parametric bootstrap confidence intervals  
- Leave one country out predictive check

### Quick example
```r
library(lme4)
m <- lmer(gini_income ~ year_c + gdp_pc_k + unemployment + openness01 + (gdp_pc_k + unemployment | country),
          data = dat, REML = FALSE)
summary(m)


