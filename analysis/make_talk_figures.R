## make_talk_figures.R -------------------------------------------------------
## Regenerates the exploratory figures used in the Senior Student Presentation
## deck. The modelling figures (beta surface, predicted survival, frailty
## histogram) come from nhanes_analysis_interaction.Rmd.
##
##   source(here::here("analysis", "make_talk_figures.R"))
##
## Outputs (all 8.0 x 4.4 in, 300 dpi, deck palette):
##   analysis/figures_talk/km_by_activity.png
##   analysis/figures_talk/aft_age_interaction.png
##   analysis/figures_talk/mims_profile_by_status.png
##   analysis/figures_talk/mims_profile_by_status_50_80.png
##   analysis/figures_talk/mims_profile_difference.png
##   analysis/figures_talk/mims_profile_two_panel.png   (8.0 x 6.4)
## ---------------------------------------------------------------------------

library(here)
library(dplyr)
library(survival)
library(survminer)

out_dir <- here("analysis", "figures_talk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data <- readRDS(here("data", "nhanes_with_demo.rds"))

data2 <- data %>%
  filter(
    age >= 50,
    if_all(c(time, event, age, gender, race, BMI, CHD, TMIMS), ~ !is.na(.)),
    time > 0,
    CHD %in% c("Yes", "No")
  ) %>%
  mutate(
    time = time / 12,
    MIMS_group = cut(
      TMIMS,
      breaks = quantile(TMIMS, c(0, 1/3, 2/3, 1)),
      labels = c("Low", "Medium", "High"),
      include.lowest = TRUE
    )
  )

## ---- numbers quoted on the data slide -------------------------------------
cat("\n--- confirm these against the deck ---\n")
cat("N analysed        :", nrow(data2), "\n")
cat("Deaths            :", sum(data2$event), sprintf("(%.1f%%)", 100 * mean(data2$event)), "\n")
cat("Design clusters   :", dplyr::n_distinct(data2$STRAPSU), "\n")
cat("Max follow-up (y) :", round(max(data2$time), 1), "\n\n")

## ---- Kaplan-Meier by tertile of total activity ----------------------------
km_fit <- survfit(Surv(time, event) ~ MIMS_group, data = data2)

p <- ggsurvplot(
  km_fit,
  data        = data2,
  conf.int    = TRUE,
  palette     = c("#C0392B", "#B58500", "#202C64"),  # deck palette
  legend.title = "Activity tertile",
  legend.labs = c("Low", "Medium", "High"),
  legend      = c(0.15, 0.25),
  xlab        = "Years since accelerometer wear",
  ylab        = "Survival probability",
  ylim        = c(0.55, 1),
  censor      = FALSE,
  ggtheme     = theme_minimal(base_size = 13)
)

ggplot2::ggsave(
  file.path(out_dir, "km_by_activity.png"),
  plot = p$plot, width = 8.0, height = 4.4, dpi = 300, bg = "white"
)

cat("Wrote", file.path(out_dir, "km_by_activity.png"), "\n")


## ===========================================================================
## Exploratory: does the PA-mortality association change with age?
## ===========================================================================
## Ordinary (non-functional, non-frailty) lognormal AFT with total daily
## activity as a single scalar predictor, fitted separately within age bands.
## If the activity effect moves across bands, that is the empirical motivation
## for letting beta vary over age in beta(s, a).
##
## Effect measure: survival-time ratio per +1 SD of TMIMS. On the AFT scale
## exp(coef) multiplies survival time directly, so 1.5 means "1.5x longer
## median survival per SD of activity". Reference line at 1 = no association.
##
## The figure is a subgroup forest plot: one row per age band plus an overall
## row. Read the interaction off whether the dots shift left or right as you
## move down the rows. No fitted trend line is drawn, deliberately -- an
## earlier version overlaid one and it disagreed with the band estimates,
## because a LINEAR age main effect lets the TMIMS x age product soak up age
## curvature (TMIMS is itself strongly age-related). The interaction test in
## the subtitle therefore uses ns(age, 3) as the main effect, so it is
## identified from how the activity slope moves.
##
## Intervals are model-based and treat participants as independent; this is an
## exploratory motivating figure, and the frailty model in the main analysis is
## what handles the within-MVU correlation properly.
##
## CAVEAT: public-use NHANES top-codes age at 80, so the "80+" band is a
## heterogeneous group recorded as exactly 80.
## ---------------------------------------------------------------------------

library(ggplot2)
library(splines)

AGE_REF <- 65   # centring value, for numerical stability only

dat_aft <- data2 %>%
  mutate(
    TMIMS_z = as.numeric(scale(TMIMS)),
    age_c   = age - AGE_REF,
    age_grp = cut(age, breaks = c(50, 60, 70, 80, Inf), right = FALSE,
                  labels = c("50-59", "60-69", "70-79", "80+"))
  )

BANDS <- levels(dat_aft$age_grp)

## ---- activity effect within each age band, and overall --------------------
## The overall row MUST adjust for age. Within a 10-year band age is nearly
## constant, so leaving it out changes almost nothing; pooled across 50-80+ it
## is a strong confounder -- older participants are both less active and closer
## to death -- and an age-unadjusted pooled coefficient absorbs that, landing
## ABOVE every band estimate. That is an aggregation artefact, not an effect.
##
## Bands get a linear age term, the pooled row gets ns(age, 3). The 80+ band is
## skipped automatically: NHANES top-codes age at 80, so age is constant there
## and the term is not identified.
fit_row <- function(dd, lab, age_adj = "none") {
  if (sum(dd$event) < 20) return(NULL)                 # too few events to fit
  rhs   <- "TMIMS_z + gender + BMI + CHD"
  n_age <- length(unique(dd$age_c))
  if (age_adj == "spline" && n_age > 10) rhs <- paste(rhs, "+ ns(age_c, df = 3)")
  if (age_adj == "linear" && n_age >  3) rhs <- paste(rhs, "+ age_c")
  f <- tryCatch(
    survreg(as.formula(paste("Surv(time, event) ~", rhs)),
            data = dd, dist = "lognormal"),
    error = function(e) NULL
  )
  if (is.null(f)) return(NULL)
  est <- coef(f)[["TMIMS_z"]]
  se  <- sqrt(vcov(f)["TMIMS_z", "TMIMS_z"])
  data.frame(row = lab, n = nrow(dd), events = sum(dd$event),
             ratio = exp(est), lower = exp(est - 1.96 * se),
             upper = exp(est + 1.96 * se))
}

forest_df <- bind_rows(
  lapply(BANDS, function(g)
    fit_row(dat_aft[dat_aft$age_grp == g, , drop = FALSE], g, age_adj = "linear")),
  list(fit_row(dat_aft, "All ages 50+", age_adj = "spline"))
) %>%
  mutate(
    kind  = ifelse(row == "All ages 50+", "overall", "band"),
    label = sprintf("%s\n%s deaths / %s", row,
                    format(events, big.mark = ","), format(n, big.mark = ",")),
    label = factor(label, levels = rev(label))
  )

## same pooled fit with no age term, to size the confounding for the console
naive <- fit_row(dat_aft, "pooled, age-unadjusted", age_adj = "none")

## ---- formal interaction test, flexible age main effect --------------------
fit_int <- survreg(
  Surv(time, event) ~ TMIMS_z + ns(age_c, df = 3) + TMIMS_z:age_c +
                      gender + BMI + CHD,
  data = dat_aft, dist = "lognormal"
)
b_int <- coef(fit_int)[["TMIMS_z:age_c"]]
p_int <- summary(fit_int)$table["TMIMS_z:age_c", "p"]

## ---- console summary ------------------------------------------------------
cat("\n--- exploratory AFT: does the TMIMS effect move with age? ---\n")
cat("Survival-time ratio per +1 SD TMIMS, lognormal AFT,\n")
cat("adjusted for sex, BMI, prevalent CHD and age.\n\n")
print(forest_df %>%
        transmute(row, n, events, ratio = round(ratio, 2),
                  CI = sprintf("[%.2f, %.2f]", lower, upper)),
      row.names = FALSE)
cat(sprintf("\nPooled estimate, age-unadjusted : %.2f [%.2f, %.2f]  <- confounded by age\n",
            naive$ratio, naive$lower, naive$upper))
k <- nrow(forest_df)
cat(sprintf("Pooled estimate, ns(age,3)      : %.2f [%.2f, %.2f]  <- in the figure\n",
            forest_df$ratio[k], forest_df$lower[k], forest_df$upper[k]))
cat("The first sits above every band because older participants are both less\n")
cat("active and closer to death; the bands are age-controlled by construction.\n")

cat(sprintf("\nInteraction (ns(age,3) main effect, linear TMIMS x age):\n"))
cat(sprintf("  %+.5f per year of age, p = %s\n\n",
            b_int, format.pval(p_int, digits = 3, eps = 1e-4)))
if (p_int > 0.05) {
  cat("NOTE: interaction is not significant once the age main effect is\n")
  cat("flexible. Present this figure as 'the effect is substantial at every\n")
  cat("age', not as an age trend.\n\n")
}

## ---- plot -----------------------------------------------------------------
p_age <- ggplot(forest_df, aes(x = ratio, y = label, colour = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 1.5, colour = "grey85") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.18, linewidth = 0.9) +
  geom_point(size = 3.4) +
  scale_colour_manual(values = c(band = "#202C64", overall = "#B58500"),
                      guide = "none") +
  scale_x_log10(breaks = c(1, 1.25, 1.5, 1.75, 2, 2.5)) +
  labs(
    x = "Survival-time ratio per +1 SD of total daily activity",
    y = NULL,
    subtitle = sprintf(
      "TMIMS x age interaction p = %s",
      format.pval(p_int, digits = 2, eps = 1e-4)
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.subtitle      = element_text(colour = "grey35", size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(colour = "grey20", lineheight = 1.05)
  )

ggsave(file.path(out_dir, "aft_age_interaction.png"),
       plot = p_age, width = 8.0, height = 4.4, dpi = 300, bg = "white")

cat("Wrote", file.path(out_dir, "aft_age_interaction.png"), "\n")


## ===========================================================================
## Exploratory: average activity profile by mortality status, within age band
## ===========================================================================
## Shows the functional covariate itself, which the KM and forest figures only
## summarise as a scalar. If the gap between decedents and survivors sits at a
## different time of day, or has a different shape, in the 70s than in the 60s,
## that is the visual case for letting beta vary over BOTH s and age.
##
## Descriptive and unadjusted on purpose: no covariates, and within a band the
## decedents are still the older half, so read this as "what the data look
## like", not as an effect estimate.
## ---------------------------------------------------------------------------

MINUTES    <- ncol(data2$MIMS)          # 1440 minute-level MIMS values
HOD        <- seq_len(MINUTES) / (MINUTES / 24)   # hour of day, 1/60 .. 24

## Two comparisons are produced: the adjacent bands 60-69 vs 70-79, and the
## extremes 50-59 vs 80+, where any age shift in the profile should be clearest.
## ---------------------------------------------------------------------------

## rows for a set of age bands; specs is a named list of c(lo, hi), inclusive.
## Index-based rather than a bind_rows of data2 subsets, because data2$MIMS is
## a matrix column and bind_rows does not carry those reliably.
band_data <- function(specs) {
  idx  <- lapply(specs, function(r) which(data2$age >= r[1] & data2$age <= r[2]))
  dd   <- data2[unlist(idx, use.names = FALSE), , drop = FALSE]
  dd$age_band <- factor(rep(names(specs), lengths(idx)), levels = names(specs))
  dd$status   <- factor(dd$event, levels = c(0, 1),
                        labels = c("Alive", "Deceased"))
  dd
}

STATUS_COLS <- c("Alive" = "#202C64", "Deceased" = "#C0392B")

profile_fig <- function(dd, file) {

  prof_df <- lapply(
    split(seq_len(nrow(dd)), list(dd$age_band, dd$status), drop = TRUE),
    function(ix) {
      M <- dd$MIMS[ix, , drop = FALSE]
      data.frame(
        hour     = HOD,
        mean     = colMeans(M, na.rm = TRUE),
        se       = apply(M, 2, sd, na.rm = TRUE) / sqrt(length(ix)),
        n        = length(ix),
        age_band = dd$age_band[ix[1]],
        status   = dd$status[ix[1]]
      )
    }
  ) %>% bind_rows()

  cat("\n--- average MIMS profile by mortality status:", file, "---\n")
  print(table(dd$age_band, dd$status))

  p <- ggplot(prof_df, aes(x = hour, y = mean, colour = status, fill = status)) +
    geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                colour = NA, alpha = 0.18) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ age_band) +
    scale_x_continuous(breaks = seq(0, 24, 6),
                       labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM"),
                       expand = expansion(mult = 0.01)) +
    scale_colour_manual(values = STATUS_COLS) +
    scale_fill_manual(values = STATUS_COLS) +
    labs(x = "Time of day", y = "Average MIMS", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", colour = "#202C64")
    )

  ggsave(file.path(out_dir, file), plot = p,
         width = 8.0, height = 4.4, dpi = 300, bg = "white")
  cat("Wrote", file.path(out_dir, file), "\n")
  invisible(prof_df)
}

## adjacent bands
profile_fig(band_data(list("Age 60-69" = c(60, 69), "Age 70-79" = c(70, 79))),
            "mims_profile_by_status.png")

## the extremes; NHANES top-codes age at 80, so "80+" is everyone at the code
profile_fig(band_data(list("Age 50-59" = c(50, 59), "Age 80+" = c(80, Inf))),
            "mims_profile_by_status_50_80.png")

cat("\nCurves are minute-level means: one point per minute of the day.\n")
cat("Bands are pointwise 95% intervals for the MEAN curve, not for an\n")
cat("individual's profile, so they are narrow by construction.\n\n")


## ===========================================================================
## Exploratory: difference activity profile (Alive - Deceased), all ages
## ===========================================================================
## One curve per age band on a single panel: the minute-by-minute gap between
## the average survivor and the average decedent. Differencing removes the
## overall diurnal shape that dominates the raw profiles and is common to both
## groups, so what is left is the part that tracks mortality. Positive means
## survivors were more active at that minute.
##
## Raw minute-level differences, no smoothing and no uncertainty bands. A
## difference of two means carries the noise of both, so expect a busy trace.
##
## Still descriptive and unadjusted: within a band the decedents are the older,
## sicker half, so this is not an effect estimate.
## ---------------------------------------------------------------------------

ALL_BANDS <- list("Age 50-59" = c(50, 59), "Age 60-69" = c(60, 69),
                  "Age 70-79" = c(70, 79), "Age 80+"   = c(80, Inf))

diff_df <- lapply(names(ALL_BANDS), function(lab) {
  dd <- band_data(ALL_BANDS[lab])
  mu_a <- colMeans(dd$MIMS[dd$status == "Alive",    , drop = FALSE], na.rm = TRUE)
  mu_d <- colMeans(dd$MIMS[dd$status == "Deceased", , drop = FALSE], na.rm = TRUE)
  data.frame(
    hour     = HOD,
    diff     = mu_a - mu_d,
    age_band = lab,
    n_alive  = sum(dd$status == "Alive"),
    n_dead   = sum(dd$status == "Deceased")
  )
}) %>%
  bind_rows() %>%
  mutate(age_band = factor(age_band, levels = names(ALL_BANDS)))

## cool -> warm encodes young -> old
AGE_COLS <- c("Age 50-59" = "#8FA8D8", "Age 60-69" = "#3C5BA9",
              "Age 70-79" = "#B58500", "Age 80+"   = "#C0392B")

## ---- console summary ------------------------------------------------------
cat("\n--- difference activity profile (Alive - Deceased) ---\n")
peak_tbl <- diff_df %>%
  group_by(age_band) %>%
  summarise(
    n_alive   = first(n_alive),
    n_dead    = first(n_dead),
    mean_gap  = round(mean(diff), 1),
    peak_gap  = round(max(diff), 1),
    peak_hour = sprintf("%02d:%02d",
                        floor(hour[which.max(diff)]),
                        round((hour[which.max(diff)] %% 1) * 60)),
    .groups   = "drop"
  )
print(peak_tbl, row.names = FALSE)
cat("\nRaw minute-level differences; peak_hour is the single noisiest kind of\n")
cat("summary there is, so read mean_gap and the shape of the curve first.\n\n")

## ---- plot -----------------------------------------------------------------
p_diff <- ggplot(diff_df, aes(x = hour, y = diff, colour = age_band)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0, 24, 6),
                     labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM"),
                     expand = expansion(mult = 0.01)) +
  scale_colour_manual(values = AGE_COLS) +
  labs(x = "Time of day",
       y = "Average MIMS: Alive - Deceased",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "mims_profile_difference.png"),
       plot = p_diff, width = 8.0, height = 4.4, dpi = 300, bg = "white")

cat("Wrote", file.path(out_dir, "mims_profile_difference.png"), "\n")

## ===========================================================================
## Smoothed two-panel figure (patchwork)
## ===========================================================================
## Panel A: average activity profile of survivors vs decedents, all ages 50+
##          pooled -- establishes the diurnal shape both groups share and the
##          level gap between them.
## Panel B: the survivor-minus-decedent difference, split into two age groups,
##          50-69 and 70-80+. Differencing removes the shape that dominates
##          panel A, so B is where a change in TIMING with age shows up.
##
## Built as two independent ggplots and stacked with patchwork, so each panel
## keeps its own legend and its own y scale.
##
## Both use geom_smooth with a CYCLIC cubic spline in time of day (bs = "cc"),
## so each fit joins across midnight instead of running free at both ends --
## the same basis the model uses for s in beta(s, a). se = FALSE.
## ---------------------------------------------------------------------------

library(mgcv)        # geom_smooth(method = "gam") needs it
library(patchwork)   # install.packages("patchwork") if missing

K_SMOOTH <- 15

TOD_SCALE <- scale_x_continuous(
  breaks = seq(0, 24, 6),
  labels = c("12 AM", "6 AM", "12 PM", "6 PM", "12 AM"),
  expand = expansion(mult = 0.01)
)

TALK_THEME <- theme_minimal(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

## ---- panel A: pooled profiles by mortality status -------------------------
dd_all <- band_data(list("All" = c(50, Inf)))

profA <- lapply(c("Alive", "Deceased"), function(st) {
  M <- dd_all$MIMS[dd_all$status == st, , drop = FALSE]
  data.frame(hour = HOD, value = colMeans(M, na.rm = TRUE), status = st)
}) %>%
  bind_rows() %>%
  mutate(status = factor(status, levels = c("Alive", "Deceased")))

p_top <- ggplot(profA, aes(x = hour, y = value, colour = status)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = K_SMOOTH),
              se = FALSE, linewidth = 1.1) +
  TOD_SCALE +
  scale_colour_manual(values = c("Alive" = "#202C64", "Deceased" = "#C0392B")) +
  labs(x = NULL, y = "Mean MIMS", colour = NULL) +
  TALK_THEME

## ---- panel B: difference, two lumped age groups ---------------------------
LUMP <- list("Age 50-69" = c(50, 69), "Age 70-80+" = c(70, Inf))

profB <- lapply(names(LUMP), function(lab) {
  dd   <- band_data(LUMP[lab])
  mu_a <- colMeans(dd$MIMS[dd$status == "Alive",    , drop = FALSE], na.rm = TRUE)
  mu_d <- colMeans(dd$MIMS[dd$status == "Deceased", , drop = FALSE], na.rm = TRUE)
  data.frame(hour = HOD, value = mu_a - mu_d, age_group = lab)
}) %>%
  bind_rows() %>%
  mutate(age_group = factor(age_group, levels = names(LUMP)))

p_bot <- ggplot(profB, aes(x = hour, y = value, colour = age_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cc", k = K_SMOOTH),
              se = FALSE, linewidth = 1.1) +
  TOD_SCALE +
  scale_colour_manual(values = c("Age 50-69"  = "#2E8B74",
                                 "Age 70-80+" = "#B58500")) +
  labs(x = "Time of day", y = "Alive - Deceased (MIMS)", colour = NULL) +
  TALK_THEME

## ---- console summary ------------------------------------------------------
cat("\n--- smoothed two-panel figure ---\n")
cat("Panel A, all ages 50+:\n")
print(table(dd_all$status))
cat("\nPanel B, group sizes:\n")
print(do.call(rbind, lapply(names(LUMP), function(lab) {
  dd <- band_data(LUMP[lab])
  data.frame(group = lab, alive = sum(dd$status == "Alive"),
             deceased = sum(dd$status == "Deceased"))
})), row.names = FALSE)
cat(sprintf("\nSmoother: cyclic cubic spline, k = %d.\n\n", K_SMOOTH))

## ---- combine and save -----------------------------------------------------
p_two <- p_top / p_bot + plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "mims_profile_two_panel.png"),
       plot = p_two, width = 8.0, height = 6.4, dpi = 300, bg = "white")

cat("Wrote", file.path(out_dir, "mims_profile_two_panel.png"), "\n")
