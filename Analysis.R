# =============================================================================
# 0. Libraries
# =============================================================================
library(tidyverse)
library(mice)
library(gtsummary)
library(broom)
library(furrr)
library(boot)
library(jtools)

# =============================================================================
# 1. Construct analytic dataset (demo.csv)
#    - Impute cognition
#    - Merge brain volumes, biomarkers, cognition
#    - Derive CSVD level, Abeta_status, vascular risk factors
# =============================================================================

# Raw data
data.raw <- read_csv("data/raw/data.csv")

# Cognition block (MMSE:CDT) – treat 0 as NA, impute with random forest
data.raw.cognition <- data.raw %>%
  select(id, MMSE:CDT)

data.raw.cognition[data.raw.cognition == 0] <- NA_real_

data.raw.cognition <- data.raw.cognition %>%
  filter(!(rowSums(is.na(.)) == (ncol(.) - 1) & !is.na(id)))

if (file.exists("data/raw/cognition_raw_imputed.csv")) {
  cognition.imputed <- read_csv("data/raw/cognition_raw_imputed.csv")
} else {
  cognition.rf <- mice::mice(
    data.raw.cognition[, -1],
    method = "rf",
    m = 100,
    maxit = 100
  )
  cognition.imputed <- cbind(
    id = data.raw.cognition$id,
    complete(cognition.rf)
  )
  write_csv(cognition.imputed, "data/raw/cognition_raw_imputed.csv")
}

cognition.imputed <- read_csv("data/raw/cognition_raw_imputed.csv")

# Brain volumes (already vscale-adjusted in your original script)
data.bv <- read_csv("data/raw/whole_brain_volume_adjusted.csv") %>%
  select(id, vscale, GMV, WMV, CSFV, ICV, WMH) %>%
  mutate(
    GMV = GMV * vscale,
    WMV = WMV * vscale,
    WMHV = WMH * vscale,
    CSFV = CSFV * vscale,
    ICV = ICV * vscale
  ) %>%
  select(id, GMV, WMV, CSFV, ICV, WMHV)

# Demo + vascular + biomarkers (“demo_all.csv” assumed to contain
# hypertension, diabetes, Abeta, etc.)
data.cvr <- read_csv("data/raw/demo_all.csv") %>%
  select(
    id, location, date_birth, date_MRI,
    sex, education,
    hypertension, diabetes,
    vascular_burden,
    Aβ42, Aβ40, ptau181, GFAP, NFL,
    MMSE:CDT
  ) %>%
  left_join(data.bv, by = "id") %>%
  mutate(
    Abeta_status = if_else(Aβ42 >= 5.13, 1L, 0L),
    `Aβ42/Aβ40` = Aβ42 / Aβ40,
    sex = if_else(sex == "M", 1L, 0L),
    date_MRI = lubridate::as_date(as.character(date_MRI)),
    date_birth = lubridate::as_date(as.character(date_birth)),
    age = as.numeric((date_MRI - date_birth) / 365)
  )

# Replace original cognition with imputed cognition
data.cvr <- data.cvr %>%
  select(-MMSE:CDT) %>%
  left_join(cognition.imputed, by = "id")

# Clean cognition (winsorisation) and compute composite domains
data.demo <- data.cvr %>%
  mutate(
    `CSVD level` = case_when(
      vascular_burden > 1 ~ "Advanced",
      vascular_burden == 1 ~ "Moderate",
      vascular_burden == 0 ~ "Mild",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(location == "GULOU") %>%
  filter(NFL <= 179) %>%
  mutate(
    VRDR  = pmin(VRDR, 14),
    VRC   = pmin(VRC, 14),
    VST1  = pmin(VST1, 60),
    VST2  = pmin(VST2, 90),
    VST3  = pmin(VST3, 120),
    TMTA  = pmin(TMTA, 360),
    TMTB  = pmin(TMTB, 720),
    BNT   = pmin(BNT, 60),
    AVLTDR = pmin(AVLTDR, 24)
  ) %>%
  mutate(
    TMTA = 1 / TMTA,
    TMTB = 1 / TMTB,
    VST1 = 1 / VST1,
    VST2 = 1 / VST2,
    VST3 = 1 / VST3
  ) %>%
  mutate(across(CVF:CDT, ~ scale(.)[, ])) %>%
  mutate(
    EM  = (AVLTDR + VRDR) / 2,
    ISP = (TMTA + VST1 + VST2) / 3,
    LF  = (BNT + CVF) / 2,
    EF  = (TMTB + VST3) / 2,
    VS  = (CDT + VRC) / 2
  ) %>%
  mutate(`ptau181/Aβ42` = ptau181 / Aβ42) %>%
  select(
    id, age, sex, education, `CSVD level`,
    hypertension, diabetes,
    GMV, WMV, CSFV, ICV, WMHV,
    Aβ42, Aβ40, `Aβ42/Aβ40`, ptau181, `ptau181/Aβ42`,
    NFL, GFAP,
    Abeta_status,
    MMSE, MoCA = MOCA, EM, ISP, LF, EF, VS
  ) %>%
  distinct() %>%
  drop_na()

write_csv(data.demo, "data/demo.csv")

# =============================================================================
# 2. Demographic table (in R only)
# =============================================================================

data.demo %>%
  select(-id) %>%
  mutate(
    `CSVD level` = factor(`CSVD level`, levels = c("Mild", "Moderate", "Advanced")),
    sex = if_else(sex == 1, "Male", "Female"),
    GMV = GMV / 1000,
    WMV = WMV / 1000,
    CSFV = CSFV / 1000,
    ICV = ICV / 1000,
    WMHV = WMHV / 1000
  ) %>%
  tbl_summary(
    by = `CSVD level`,
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    digits = all_continuous() ~ 2,
    missing_text = "Missing"
  ) %>%
  add_p(test = all_continuous() ~ "aov") %>%
  add_overall(last = TRUE)

# =============================================================================
# 3. Partial correlations & raw correlations
#    - residualise markers on age + sex + education + hypertension + diabetes + ICV + Abeta_status
# =============================================================================

data.all <- read_csv("data/demo.csv") %>% distinct()

markers <- c(
  "GMV", "WMV", "WMHV",
  "Aβ42", "Aβ40", "Aβ42/Aβ40",
  "GFAP", "NFL", "ptau181",
  "MMSE", "MoCA", "EM", "ISP", "LF", "EF", "VS"
)

# Residuals (partialled out covariates)
data.resid <- data.all %>%
  pivot_longer(all_of(markers), names_to = "marker_label", values_to = "marker") %>%
  group_by(marker_label) %>%
  mutate(
    residuals = unname(
      resid(
        lm(
          marker ~ age + sex + education +
            hypertension + diabetes + ICV + Abeta_status,
          data = cur_data()
        )
      )
    )
  ) %>%
  ungroup() %>%
  select(-marker) %>%
  pivot_wider(names_from = "marker_label", values_from = "residuals")

# Marker combinations
marker_combinations <- combinat::combn(markers, 2) %>%
  t() %>%
  as_tibble() %>%
  setNames(c("marker1", "marker2"))

# Partial correlations (by CSVD level)
result.partial.cor <-
  marker_combinations %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(burden_level = rep(c("Mild", "Moderate", "Advanced"), nrow(marker_combinations))) %>%
  mutate(
    stats = pmap(
      list(burden_level, marker1, marker2),
      function(level_, m1_, m2_) {
        df <- data.resid %>% filter(`CSVD level` == level_)
        x  <- df %>% pull(m1_)
        y  <- df %>% pull(m2_)
        coef.cor <- cor.test(x, y) %>% broom::tidy()

        tibble(
          N     = nrow(df),
          r     = coef.cor$estimate,
          z     = 0.5 * log((1 + r) / (1 - r)), # Fisher r-to-z
          t     = coef.cor$statistic,
          p     = coef.cor$p.value,
          df    = coef.cor$parameter,
          r_low = coef.cor$conf.low,
          r_high = coef.cor$conf.high
        )
      }
    )
  ) %>%
  unnest(stats)

write_csv(result.partial.cor, "results/200-cor_diff/cor_partial.csv")

# Raw correlations (no residualisation)
result.cor <-
  marker_combinations %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(burden_level = rep(c("Mild", "Moderate", "Advanced"), nrow(marker_combinations))) %>%
  mutate(
    stats = pmap(
      list(burden_level, marker1, marker2),
      function(level_, m1_, m2_) {
        df <- data.all %>% filter(`CSVD level` == level_)
        x  <- df %>% pull(m1_)
        y  <- df %>% pull(m2_)
        coef.cor <- cor.test(x, y) %>% broom::tidy()

        tibble(
          N     = nrow(df),
          r     = coef.cor$estimate,
          z     = 0.5 * log((1 + r) / (1 - r)),
          t     = coef.cor$statistic,
          p     = coef.cor$p.value,
          df    = coef.cor$parameter,
          r_low = coef.cor$conf.low,
          r_high = coef.cor$conf.high
        )
      }
    )
  ) %>%
  unnest(stats)

write_csv(result.cor, "results/200-cor_diff/cor_raw.csv")

# =============================================================================
# 4. Differences between correlations (Fisher z)
# =============================================================================

diff_r <- function(z1, z2, n1, n2, prefix = "") {
  se_diff <- sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
  z_diff  <- (z1 - z2) / se_diff
  p_value <- 2 * (1 - pnorm(abs(z_diff)))

  tibble(
    !!str_c(prefix, "se_diff") := se_diff,
    !!str_c(prefix, "z_diff")  := z_diff,
    !!str_c(prefix, "p_value") := p_value
  )
}

data.cor <- read_csv("results/200-cor_diff/cor_raw.csv")
diff.nest <- data.cor %>%
  group_by(marker1, marker2) %>%
  nest()

result.diff <- pmap_dfr(
  list(diff.nest$marker1, diff.nest$marker2, diff.nest$data),
  function(m1, m2, df) {
    df <- df %>%
      mutate(burden_level = factor(burden_level,
                                   levels = c("Mild", "Moderate", "Advanced")))

    advanced <- df %>% filter(burden_level == "Advanced")
    moderate <- df %>% filter(burden_level == "Moderate")
    mild     <- df %>% filter(burden_level == "Mild")

    cbind(
      marker1 = m1,
      marker2 = m2,
      diff_r(advanced$z, moderate$z, advanced$N, moderate$N, "A_Mo_"),
      diff_r(advanced$z, mild$z,     advanced$N, mild$N,     "A_Mi_"),
      diff_r(moderate$z, mild$z,     moderate$N, mild$N,     "Mo_Mi_")
    )
  }
)

write_csv(result.diff, "results/200-cor_diff/diff_cor_raw.csv")

# =============================================================================
# 5. Interaction linear models (global markers × CSVD/WMHV; bootstrap slopes)
# =============================================================================

data.all <- read_csv("data/demo.csv") %>%
  distinct() %>%
  mutate(
    `CSVD level` = case_when(
      `CSVD level` == "Advanced" ~ 2L,
      `CSVD level` == "Moderate" ~ 1L,
      `CSVD level` == "Mild"     ~ 0L,
      TRUE ~ NA_integer_
    ),
    `WMHV level` = ntile(WMHV, 3L) - 1L
  )

burden_type <- c("CSVD level", "WMHV level")
markers_ad <- c("Aβ42", "Aβ40", "Aβ42/Aβ40", "ptau181")
markers_response <- c(
  "GMV", "WMV", "GFAP", "WMHV", "NFL",
  "MMSE", "MoCA", "EM", "ISP", "LF", "EF", "VS"
)

# Function used by boot
get_slopes <- function(data, indices) {
  d <- data[indices, ]
  model <- lm(
    Y ~ age + sex + education + hypertension + diabetes + ICV +
      Abeta_status + X * burden,
    data = d
  )
  tidy_model <- broom::tidy(model)

  intercept_X <- tidy_model %>% filter(term == "X") %>% pull(estimate)
  slope_burden0 <- intercept_X
  slope_burden1 <- intercept_X + tidy_model %>% filter(term == "X:burden1") %>% pull(estimate)
  slope_burden2 <- intercept_X + tidy_model %>% filter(term == "X:burden2") %>% pull(estimate)

  c(slope_burden0, slope_burden1, slope_burden2)
}

result.comb <- expand.grid(
  burden_type      = burden_type,
  markers_ad       = markers_ad,
  markers_response = markers_response
) %>%
  as_tibble()

future::plan(future::multisession, workers = 6)

result.slope <- furrr::future_pmap_dfr(
  list(result.comb$burden_type, result.comb$markers_ad, result.comb$markers_response),
  function(iBurden, iAD, iResponse) {

    df <- data.all %>%
      select(
        Y = all_of(iResponse),
        X = all_of(iAD),
        age, sex, education,
        hypertension, diabetes, ICV, Abeta_status,
        burden = all_of(iBurden)
      ) %>%
      filter(
        abs((Y - mean(Y, na.rm = TRUE)) / sd(Y, na.rm = TRUE)) <= 3 |
          abs((X - mean(X, na.rm = TRUE)) / sd(X, na.rm = TRUE)) <= 3
      ) %>%
      mutate(
        X         = scale(X)[, ],
        Y         = scale(Y)[, ],
        age       = scale(age)[, ],
        education = scale(education)[, ],
        ICV       = scale(ICV)[, ],
        burden    = factor(burden)
      )

    # Original slopes from the full model
    model <- lm(
      Y ~ age + sex + education + hypertension + diabetes + ICV +
        Abeta_status + X * burden,
      data = df
    )
    tidy_model <- broom::tidy(model)

    intercept_X <- tidy_model %>% filter(term == "X") %>% pull(estimate)
    slope_burden0 <- intercept_X
    slope_burden1 <- intercept_X + tidy_model %>% filter(term == "X:burden1") %>% pull(estimate)
    slope_burden2 <- intercept_X + tidy_model %>% filter(term == "X:burden2") %>% pull(estimate)

    # Bootstrap
    set.seed(1234)
    results <- boot::boot(df, get_slopes, R = 5000)

    bootstrap_slopes <- results$t
    colnames(bootstrap_slopes) <- c("burden0", "burden1", "burden2")

    bootstrap_diff <- as_tibble(bootstrap_slopes) %>%
      mutate(
        diff_burden1_burden0 = burden1 - burden0,
        diff_burden2_burden0 = burden2 - burden0,
        diff_burden2_burden1 = burden2 - burden1
      )

    summary_stats <- bootstrap_diff %>%
      summarise(
        mean_burden0 = mean(burden0),
        lower_ci_burden0 = quantile(burden0, 0.025),
        upper_ci_burden0 = quantile(burden0, 0.975),

        mean_burden1 = mean(burden1),
        lower_ci_burden1 = quantile(burden1, 0.025),
        upper_ci_burden1 = quantile(burden1, 0.975),

        mean_burden2 = mean(burden2),
        lower_ci_burden2 = quantile(burden2, 0.025),
        upper_ci_burden2 = quantile(burden2, 0.975),

        mean_diff_burden1_burden0 = mean(diff_burden1_burden0),
        lower_ci_diff_burden1_burden0 = quantile(diff_burden1_burden0, 0.025),
        upper_ci_diff_burden1_burden0 = quantile(diff_burden1_burden0, 0.975),

        mean_diff_burden2_burden0 = mean(diff_burden2_burden0),
        lower_ci_diff_burden2_burden0 = quantile(diff_burden2_burden0, 0.025),
        upper_ci_diff_burden2_burden0 = quantile(diff_burden2_burden0, 0.975),

        mean_diff_burden2_burden1 = mean(diff_burden2_burden1),
        lower_ci_diff_burden2_burden1 = quantile(diff_burden2_burden1, 0.025),
        upper_ci_diff_burden2_burden1 = quantile(diff_burden2_burden1, 0.975)
      )

    tibble(
      vascular_type = iBurden,
      marker_AD = iAD,
      marker_response = iResponse,
      slope_mild      = slope_burden0,
      slope_moderate  = slope_burden1,
      slope_advanced  = slope_burden2
    ) %>%
      bind_cols(summary_stats)
  }
)

write_csv(result.slope, "results/201-slope-diff/result_slope.csv")

# =============================================================================
# 6. Regional linear models (no plots, FDR-adjusted)
# =============================================================================

# Base AD + confounders
data.base <- read_csv("data/demo.csv") %>%
  distinct() %>%
  select(
    id, age, sex, education, WMHV,
    Aβ42, Aβ40, `Aβ42/Aβ40`, ptau181,
    hypertension, diabetes, ICV, Abeta_status
  ) %>%
  pivot_longer(c("Aβ42", "Aβ40", "Aβ42/Aβ40", "ptau181"),
               values_to = "AD", names_to = "AD_marker")

# Regional IDPs
data.regional <- read_csv("data/idp_regional.csv")

data.all.reg <- data.base %>%
  left_join(data.regional, relationship = "many-to-many") %>%
  drop_na()

data.nest.reg <- data.all.reg %>%
  group_by(AD_marker, measure, region) %>%
  nest()

future::plan(future::multisession, workers = 6)

result.lm.reg <- furrr::future_pmap_dfr(
  list(
    data.nest.reg$AD_marker,
    data.nest.reg$measure,
    data.nest.reg$region,
    data.nest.reg$data
  ),
  function(iMarker, iMeasure, iRegion, iDF) {

    # Standardised variables
    iDF.std <- jtools::standardize(
      iDF,
      vars = c("value", "age", "education", "WMHV", "AD", "ICV")
    )

    lm.fit.std <- lm(
      value ~ age + sex + education + hypertension + diabetes +
        ICV + Abeta_status + WMHV * AD,
      data = iDF.std
    )

    lm.fit <- lm(
      value ~ age + sex + education + hypertension + diabetes +
        ICV + Abeta_status + WMHV * AD,
      data = iDF
    )

    lm.coef     <- broom::tidy(lm.fit)
    lm.coef.std <- broom::tidy(lm.fit.std)

    tibble(
      marker = iMarker,
      measure = iMeasure,
      region = iRegion,

      beta_AD           = lm.coef.std %>% filter(term == "AD")        %>% pull(estimate),
      beta_se_AD        = lm.coef.std %>% filter(term == "AD")        %>% pull(std.error),
      beta_t_AD         = lm.coef.std %>% filter(term == "AD")        %>% pull(statistic),
      beta_p_AD         = lm.coef.std %>% filter(term == "AD")        %>% pull(p.value),

      beta_WMHV         = lm.coef.std %>% filter(term == "WMHV")      %>% pull(estimate),
      beta_se_WMHV      = lm.coef.std %>% filter(term == "WMHV")      %>% pull(std.error),
      beta_t_WMHV       = lm.coef.std %>% filter(term == "WMHV")      %>% pull(statistic),
      beta_p_WMHV       = lm.coef.std %>% filter(term == "WMHV")      %>% pull(p.value),

      beta_interaction  = lm.coef.std %>% filter(term == "WMHV:AD")   %>% pull(estimate),
      beta_se_interaction = lm.coef.std %>% filter(term == "WMHV:AD") %>% pull(std.error),
      beta_t_interaction  = lm.coef.std %>% filter(term == "WMHV:AD") %>% pull(statistic),
      beta_p_interaction  = lm.coef.std %>% filter(term == "WMHV:AD") %>% pull(p.value),

      B_AD              = lm.coef %>% filter(term == "AD")        %>% pull(estimate),
      B_se_AD           = lm.coef %>% filter(term == "AD")        %>% pull(std.error),
      B_t_AD            = lm.coef %>% filter(term == "AD")        %>% pull(statistic),
      B_p_AD            = lm.coef %>% filter(term == "AD")        %>% pull(p.value),

      B_WMHV            = lm.coef %>% filter(term == "WMHV")      %>% pull(estimate),
      B_se_WMHV         = lm.coef %>% filter(term == "WMHV")      %>% pull(std.error),
      B_t_WMHV          = lm.coef %>% filter(term == "WMHV")      %>% pull(statistic),
      B_p_WMHV          = lm.coef %>% filter(term == "WMHV")      %>% pull(p.value),

      B_interaction     = lm.coef %>% filter(term == "WMHV:AD")   %>% pull(estimate),
      B_se_interaction  = lm.coef %>% filter(term == "WMHV:AD")   %>% pull(std.error),
      B_t_interaction   = lm.coef %>% filter(term == "WMHV:AD")   %>% pull(statistic),
      B_p_interaction   = lm.coef %>% filter(term == "WMHV:AD")   %>% pull(p.value)
    )
  }
)

# FDR-adjustment on p-values, excluding cerebellar HO regions
result.lm.reg.adj <- result.lm.reg %>%
  filter(region < 111) %>%
  group_by(measure, marker) %>%
  mutate(
    across(
      c(
        beta_p_AD, beta_p_WMHV, beta_p_interaction,
        B_p_AD, B_p_WMHV, B_p_interaction
      ),
      ~ p.adjust(.x, method = "fdr"),
      .names = "{col}_fdr"
    )
  ) %>%
  ungroup()

write_csv(result.lm.reg.adj, "results/300-regional linear models/results_all_FDR.csv")