suppressPackageStartupMessages({
  library(here)
  library(splines)
  library(dplyr)
  library(arrow)
})

fit_age_eigengene_models <- function(df, eigengene, covariates = character(), spline_df = 4) {
  covar <- paste(covariates, collapse = " + ")
  rhs_linear <- paste(c("Age", covar), collapse = " + ")
  rhs_spline <- paste(c(sprintf("ns(Age, df = %d)", spline_df), covar), collapse = " + ")
  linear <- lm(as.formula(paste(eigengene, "~", rhs_linear)), data = df)
  spline <- lm(as.formula(paste(eigengene, "~", rhs_spline)), data = df)
  cmp <- anova(linear, spline)
  tibble(
    eigengene = eigengene,
    linear_age_p = coef(summary(linear))["Age", "Pr(>|t|)"],
    spline_vs_linear_p = cmp$`Pr(>F)`[2],
    linear_r2 = summary(linear)$r.squared,
    spline_r2 = summary(spline)$r.squared,
    delta_r2 = summary(spline)$r.squared - summary(linear)$r.squared
  )
}

message("Use this script after module eigengenes are exported to parquet.")
