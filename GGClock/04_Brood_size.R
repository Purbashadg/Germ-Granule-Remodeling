## ============================================================
## BATCH 1 — Reproductive output and the GG-clock (Days 1-6)
## ============================================================
pk <- c("glmmTMB", "emmeans", "ggplot2", "patchwork", "readxl")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(glmmTMB); library(emmeans); library(ggplot2); library(patchwork) })
set.seed(19)

## ---- PATHS ----
BASE       <- "/Users/purbasha/Desktop/Germ granules/Age clock/Analysis_raw image"
RUN_TAG    <- "run2"   
PROG_FILE  <- file.path(BASE, "Progeny count.xlsx")
B1_OUT_DIR <- file.path(BASE, paste0("results_v6_B1_D1to6_", RUN_TAG),
                        "Batch1_results_D1to6")
CLK_FILE   <- file.path(B1_OUT_DIR, "Batch1_LOOCV_predictions.csv")
PROG_OUT   <- file.path(BASE, "progeny_results_B1_D1to6")
if (!file.exists(CLK_FILE)) stop("Clock file not found:\n  ", CLK_FILE, call. = FALSE)
if (!dir.exists(PROG_OUT)) dir.create(PROG_OUT, recursive = TRUE)

## ---- 1. LOAD PROGENY (Batch 1 only) ----
raw <- as.data.frame(readxl::read_excel(PROG_FILE))
pick <- function(pat) {
  j <- grep(pat, names(raw), ignore.case = TRUE)[1]
  if (is.na(j)) stop("No column matching ", pat, ". Columns: ", paste(names(raw), collapse = ", "))
  j
}
prog1 <- data.frame(
  worm  = as.character(raw[[pick("^worm")]]),
  day   = as.integer(gsub("[^0-9]", "", raw[[pick("^age")]])),
  prog  = as.numeric(raw[[pick("^progeny")]]),
  batch = gsub("[^0-9]", "", as.character(raw[[pick("^batch")]])),
  stringsAsFactors = FALSE)
prog1 <- prog1[prog1$batch == "1" & is.finite(prog1$prog) & is.finite(prog1$day), ]
prog1$worm <- factor(prog1$worm)
ALL_DAYS <- sort(unique(prog1$day))
cat(sprintf("Batch 1 progeny: %d rows, %d worms, days %s\n",
            nrow(prog1), nlevels(prog1$worm), paste(ALL_DAYS, collapse = ",")))

summ <- function(df, v) do.call(rbind, lapply(split(df, df$day), function(x)
  data.frame(day = x$day[1], m = mean(x[[v]]), sd = sd(x[[v]]), n = nrow(x))))

## ---- 2. DAILY DESCRIPTIVES ----
desc <- do.call(rbind, lapply(split(prog1, prog1$day), function(x)
  data.frame(day = x$day[1], n = nrow(x), mean = mean(x$prog), sd = sd(x$prog),
             median = median(x$prog), se = sd(x$prog) / sqrt(nrow(x)))))
cat("\n--- Daily progeny per worm ---\n"); print(round(desc, 2), row.names = FALSE)
write.csv(desc, file.path(PROG_OUT, "B1_progeny_descriptives.csv"), row.names = FALSE)

## ---- 3. CUMULATIVE OUTPUT ----
prog1 <- prog1[order(prog1$worm, prog1$day), ]
prog1$cum <- ave(prog1$prog, prog1$worm, FUN = cumsum)
n_days   <- tapply(prog1$day, prog1$worm, function(x) length(unique(x)))
complete <- names(n_days)[n_days == length(ALL_DAYS)]
cumdat   <- prog1[prog1$worm %in% complete, ]
cat(sprintf("\nWorms with complete Days %d-%d records: %d of %d\n",
            min(ALL_DAYS), max(ALL_DAYS), length(complete), nlevels(prog1$worm)))

brood <- aggregate(prog ~ worm, cumdat, sum); names(brood)[2] <- "brood"
cat(sprintf("Lifetime brood: mean = %.1f | SD = %.1f | median = %.1f | range %d-%d (n = %d)\n",
            mean(brood$brood), sd(brood$brood), median(brood$brood),
            min(brood$brood), max(brood$brood), nrow(brood)))

cum_d <- summ(cumdat, "cum")
cum_d$pct_of_brood <- 100 * cum_d$m / mean(brood$brood)
cat("\n--- Cumulative progeny per worm ---\n"); print(round(cum_d, 1), row.names = FALSE)
write.csv(cum_d, file.path(PROG_OUT, "B1_cumulative_progeny.csv"), row.names = FALSE)
write.csv(brood, file.path(PROG_OUT, "B1_lifetime_brood_per_worm.csv"), row.names = FALSE)

bl <- setNames(brood$brood, as.character(brood$worm))
cumdat$brood <- bl[as.character(cumdat$worm)]
if (any(cumdat$brood == 0)) cat("Note: excluded", length(unique(cumdat$worm[cumdat$brood == 0])),
                                "worm(s) with zero lifetime brood from the % panel\n")
cumpct <- cumdat[cumdat$brood > 0, ]
cumpct$pct <- 100 * cumpct$cum / cumpct$brood
pct_d <- summ(cumpct, "pct")
cat("\n--- Cumulative progeny, % of own lifetime brood ---\n")
print(round(pct_d, 1), row.names = FALSE)

PLATEAU_PCT   <- 95
plateau_start <- min(pct_d$day[pct_d$m >= PLATEAU_PCT])
cat(sprintf("Plateau (>= %d%% of lifetime brood) reached on Day %d\n", PLATEAU_PCT, plateau_start))

## ---- 4. GLMM on daily progeny: Days 1-5 (Day 6 all zeros) ----
MODEL_DAYS <- 1:5
pm <- prog1[prog1$day %in% MODEL_DAYS, ]
pm$day_f <- factor(pm$day, levels = MODEL_DAYS)
m1 <- glmmTMB(prog ~ day_f + (1 | worm), family = nbinom2, data = pm)
m0 <- glmmTMB(prog ~ 1      + (1 | worm), family = nbinom2, data = pm)
cat("\n--- Negative binomial GLMM, Days 1-5 ---\n")
cat("Converged (pdHess):", m1$sdr$pdHess, "\n")
lrt <- anova(m0, m1); print(lrt)
p_lrt <- lrt$`Pr(>Chisq)`[2]
cat("\nPairwise vs Day 1 (Dunnett, response ratios):\n")
em <- emmeans(m1, ~ day_f)
print(summary(contrast(em, "trt.vs.ctrl", ref = 1), type = "response", infer = TRUE))
write.csv(as.data.frame(summary(em, type = "response")),
          file.path(PROG_OUT, "B1_progeny_GLMM_emmeans.csv"), row.names = FALSE)
p_txt <- ifelse(p_lrt < 0.001, "< 0.001", sprintf("= %.3f", p_lrt))
# LRT
lrt_df <- data.frame(
  test   = "Day effect (LRT vs intercept-only), Days 1-5",
  Chisq  = lrt$Chisq[2],
  df     = lrt$`Df`[2] - lrt$`Df`[1],
  P      = p_lrt)
write.csv(lrt_df, file.path(PROG_OUT, "B1_progeny_GLMM_LRT.csv"), row.names = FALSE)

# Dunnett contrasts vs Day 1 (response ratios, 95% CI)
dun <- as.data.frame(summary(contrast(em, "trt.vs.ctrl", ref = 1),
                             type = "response", infer = TRUE))
write.csv(dun, file.path(PROG_OUT, "B1_progeny_GLMM_Dunnett_vs_Day1.csv"), row.names = FALSE)
print(lrt_df); print(dun)

## ---- 5. CLOCK (Batch 1 LOO-CV, Days 1-6) ----
clk <- read.csv(CLK_FILE)
clk <- data.frame(day = clk$actual_day, pred = clk$predicted_day)
clk <- clk[is.finite(clk$day) & is.finite(clk$pred), ]

CLK_DAYS <- sort(unique(clk$day))
miss_clk <- setdiff(ALL_DAYS, CLK_DAYS)
if (length(miss_clk))
  stop("Clock file has no data for day(s): ", paste(miss_clk, collapse = ", "),
       ". Re-run the Batch 1 pipeline with those days included.", call. = FALSE)
cat("Clock days present:", paste(CLK_DAYS, collapse = ","), "\n")

r2   <- cor(clk$day, clk$pred)^2
mae  <- mean(abs(clk$pred - clk$day))
rmse <- sqrt(mean((clk$pred - clk$day)^2))
rho  <- cor(clk$day, clk$pred, method = "spearman")
cat(sprintf("\n--- Clock Batch 1 LOO-CV: Pearson r\u00b2 = %.3f | MAE = %.2f d | RMSE = %.2f d | rho = %.3f | n = %d ---\n",
            r2, mae, rmse, rho, nrow(clk)))
ck_d <- summ(clk, "pred")
cat("\n--- Predicted age per day (mean, SD, n) ---\n"); print(round(ck_d, 2), row.names = FALSE)

LAST_DAY <- max(CLK_DAYS)
Y_CLK <- max(6, ceiling(max(clk$pred, ck_d$m + ck_d$sd, na.rm = TRUE)))   # clock axis max (days)

## ============================================================
## 6. FIGURE 
## ============================================================
FIG_FONT <- "Times New Roman"
FS       <- 12
ANN      <- FS / .pt          # in-plot text size (mm) = 12 pt
LW       <- 0.5
PT_SZ    <- 1.8

pub_theme <- theme_classic(base_size = FS, base_family = FIG_FONT) +
  theme(
    text        = element_text(family = FIG_FONT, size = FS, face = "bold", colour = "black"),
    axis.title  = element_text(family = FIG_FONT, size = FS, face = "bold", colour = "black"),
    axis.text   = element_text(family = FIG_FONT, size = FS, face = "bold", colour = "black"),
    axis.line   = element_line(linewidth = LW, colour = "black"),
    axis.ticks  = element_line(linewidth = LW, colour = "black"),
    legend.text = element_text(family = FIG_FONT, size = FS, face = "bold", colour = "black"),
    plot.title = element_blank(), plot.subtitle = element_blank(),
    plot.caption = element_blank())
theme_set(pub_theme)
update_geom_defaults("text", list(family = FIG_FONT, size = ANN, colour = "black", fontface = "bold"))

W_MAIN <- 180; H_MAIN <- 130

save_fig <- function(p, name, w = W_MAIN, h = H_MAIN) {
  ggsave(file.path(PROG_OUT, paste0(name, ".tiff")), p, device = "tiff",
         width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
  ggsave(file.path(PROG_OUT, paste0(name, ".pdf")), p, device = cairo_pdf,
         width = w, height = h, units = "mm")
}

LEG_FILE <- file.path(PROG_OUT, "Figure_legends.txt")
writeLines(character(0), LEG_FILE)
register_legend <- function(n, name, text) {
  cat(sprintf("Figure %d (%s). %s\n\n", n, name,
              paste(strwrap(text, 100), collapse = "\n")),
      file = LEG_FILE, append = TRUE)
}

cP <- "firebrick"; cK <- "steelblue"
XLIM <- c(min(ALL_DAYS) - 0.4, max(ALL_DAYS) + 0.4)
SEP_TXT <- paste0("Progeny and clock data come from separate animals; the comparison is ",
                  "cohort-level and no correlation was computed.")
DAY_RNG <- sprintf("D%d\u2013D%d", min(CLK_DAYS), LAST_DAY)

## ============================================================
## ============================================================
## FIGURE 1: Daily progeny vs chronological age
## ============================================================
fig_daily <- ggplot(prog1, aes(day, prog)) +
  geom_jitter(width = 0.15, height = 0, size = PT_SZ, alpha = 0.5, colour = "grey50") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12, colour = cP, linewidth = LW) +
  stat_summary(fun = mean, geom = "line",  colour = cP, linewidth = 0.6) +
  stat_summary(fun = mean, geom = "point", colour = cP, size = 3.5) +
  scale_x_continuous(breaks = ALL_DAYS) +
  coord_cartesian(xlim = XLIM) +
  labs(x = "Chronological age (days)", y = "Progeny per worm per day")
save_fig(fig_daily, "Figure_1_daily_progeny_vs_chronological_age")
register_legend(1, "daily progeny", paste0(
  "Daily progeny by chronological age (Batch 1). Gray points: individual worms; red: mean \u00b1 SEM. ",
  "Negative binomial GLMM (Days 1\u20135, worm as random intercept): day effect likelihood-ratio *P* ", p_txt,
  ". Day 6 (all zeros) is shown descriptively and excluded from the model."))
## ============================================================
## FIGURE 2: Cumulative progeny vs chronological age
## ============================================================
fig_cum <- ggplot() +
  geom_line(data = cumdat, aes(day, cum, group = worm), colour = "grey60", alpha = 0.45, linewidth = 0.4) +
  geom_ribbon(data = cum_d, aes(day, ymin = pmax(m - sd, 0), ymax = m + sd), fill = cP, alpha = 0.15) +
  geom_line(data = cum_d, aes(day, m), colour = cP, linewidth = 1.1) +
  geom_point(data = cum_d, aes(day, m), colour = cP, size = 3) +
  scale_x_continuous(breaks = ALL_DAYS) +
  coord_cartesian(xlim = XLIM) +
  labs(x = "Chronological age (days)", y = "Cumulative progeny per worm")
save_fig(fig_cum, "Figure_2_cumulative_progeny_vs_chronological_age")
register_legend(2, "cumulative progeny", sprintf(paste0(
  "Cumulative progeny by chronological age (Batch 1). Grey lines: individual worms; red line and point: mean; ",
  "shaded band: mean \u00b1 SD. Lifetime brood = %.0f \u00b1 %.0f (mean \u00b1 SD), n = %d worms with complete records."),
  mean(brood$brood), sd(brood$brood), nrow(brood)))

## ============================================================
## FIGURE 3: Clock performance (LOO-CV, Days 1-6)
## ============================================================
clock_txt <- sprintf("n = %d\nMAE = %.2f d\nRMSE = %.2f d\nPearson r\u00b2 = %.2f\n\u03c1 = %.2f",
                     nrow(clk), mae, rmse, r2, rho)
ax_max <- max(ALL_DAYS, Y_CLK)

fig_clock <- ggplot(clk, aes(day, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = LW) +
  geom_jitter(width = 0.15, height = 0, size = PT_SZ, alpha = 0.5, colour = "grey50") +
  stat_summary(fun = median, geom = "point", size = 3.5, colour = cP) +
  annotate("text", x = -Inf, y = Inf, label = clock_txt, hjust = -0.08, vjust = 1.1,
           family = FIG_FONT, fontface = "bold", size = ANN, lineheight = 0.95) +
  scale_x_continuous(breaks = ALL_DAYS) + scale_y_continuous(breaks = 1:ax_max) +
  coord_cartesian(xlim = c(0.7, max(ALL_DAYS) + 0.3), ylim = c(0.7, ax_max + 0.5)) +
  labs(x = "Chronological age (days)", y = "Predicted age (days)")
save_fig(fig_clock, "Figure_3_clock_predicted_vs_chronological_age")
register_legend(3, "clock performance", sprintf(paste0(
  "Germ-granule age clock, Batch 1, leave-one-out cross-validation (Days %s). Grey points: individual worms; red: median; ",
  "dashed line: predicted equals chronological age. MAE = %.2f d, RMSE = %.2f d, r\u00b2 = %.2f, ",
  "Spearman \u03c1 = %.2f, n = %d."),
  paste(range(CLK_DAYS), collapse = "-"), mae, rmse, r2, rho, nrow(clk)))

## ============================================================
## Stacked figures: output (top), predicted age (bottom)
## ============================================================
stacked_plot <- function(ind, mean_df, ylab_top, top_lines = FALSE,
                         band_start = NULL, ylim_top = NULL) {
  band <- if (is.null(band_start)) NULL else
    annotate("rect", xmin = band_start, xmax = XLIM[2], ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.35)
  p_top <- ggplot() + band +
    (if (top_lines) geom_line(data = ind, aes(day, y, group = worm), colour = "grey60", alpha = 0.4, linewidth = 0.4)
     else geom_jitter(data = ind, aes(day, y), width = 0.12, height = 0, size = PT_SZ, alpha = 0.4, colour = "grey50")) +
    geom_line(data = mean_df,  aes(day, m), colour = cP, linewidth = 1) +
    geom_point(data = mean_df, aes(day, m), colour = cP, size = 3) +
    scale_x_continuous(breaks = ALL_DAYS) +
    coord_cartesian(xlim = XLIM, ylim = ylim_top) +
    labs(y = ylab_top) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p_bot <- ggplot() + band +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = LW) +
    geom_jitter(data = clk, aes(day, pred), width = 0.12, height = 0, size = PT_SZ, alpha = 0.4, colour = "lightsteelblue3") +
    geom_errorbar(data = ck_d, aes(day, ymin = m - sd, ymax = m + sd), width = 0.12, colour = cK, linewidth = LW) +
    geom_line(data = ck_d,  aes(day, m), colour = cK, linewidth = 1) +
    geom_point(data = ck_d, aes(day, m), colour = cK, size = 3, shape = 17) +
    scale_x_continuous(breaks = ALL_DAYS) +
    coord_cartesian(xlim = XLIM, ylim = c(0.5, Y_CLK)) +
    labs(x = "Chronological age (days)", y = "Clock-predicted age (days)")
  
  (p_top / p_bot) + plot_layout(heights = c(1, 1))
}

## FIGURE 4: Daily progeny + predicted age
ind_daily <- data.frame(day = prog1$day, y = prog1$prog, worm = prog1$worm)
fig_daily_stack <- stacked_plot(ind_daily, summ(prog1, "prog"), "Progeny per worm per day")
save_fig(fig_daily_stack, "Figure_4_daily_progeny_and_predicted_age_stacked")
register_legend(4, "daily progeny and predicted age", paste0(
  "Daily progeny and clock-predicted age (Batch 1). Top: individual worms (grey) and mean (red). ",
  "Bottom: individual imaged worms (light blue) and mean \u00b1 SD (blue); dashed line: predicted equals chronological age. ",
  SEP_TXT))

## FIGURE 5: Cumulative % of brood + predicted age
ind_cum <- data.frame(day = cumpct$day, y = cumpct$pct, worm = cumpct$worm)
fig_cum_stack <- stacked_plot(ind_cum, pct_d, "Cumulative progeny (% of brood)",
                              top_lines = TRUE, band_start = plateau_start, ylim_top = c(0, 105))
save_fig(fig_cum_stack, "Figure_5_cumulative_progeny_and_predicted_age_stacked")
register_legend(5, "cumulative progeny and predicted age", paste0(
  "Cumulative progeny and clock-predicted age (Batch 1). Top: each worm's cumulative progeny as a percentage of its own ",
  "lifetime brood (grey) and mean (red). Bottom: individual imaged worms (light blue) and mean \u00b1 SD (blue); ",
  "dashed line: predicted equals chronological age. Shaded region: Day ", plateau_start,
  " onward, when mean cumulative output is \u2265 ", PLATEAU_PCT, "% of lifetime brood. ", SEP_TXT))

## ---- trajectory plot ----
traj_plot <- function(prog_d, ylab) {
  j <- merge(prog_d[, c("day", "m")], ck_d[, c("day", "m")], by = "day", suffixes = c("_prog", "_pred"))
  ggplot(j, aes(m_pred, m_prog)) +
    geom_path(colour = "grey50", linewidth = 0.7,
              arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    geom_point(size = 3.5, colour = cP) +
    geom_text(aes(label = paste0("D", day)), nudge_y = max(j$m_prog) * 0.05,
              family = FIG_FONT, fontface = "bold", size = ANN, colour = "black") +
    coord_cartesian(clip = "off") +
    labs(x = "Mean clock-predicted age (days)", y = ylab)
}

## FIGURE 6
fig_daily_traj <- traj_plot(summ(prog1, "prog"), "Mean progeny per worm per day")
save_fig(fig_daily_traj, "Figure_6_daily_progeny_vs_predicted_age")
register_legend(6, "daily progeny vs predicted age", paste0(
  "Mean daily progeny versus mean clock-predicted age (Batch 1). Each point is one chronological day (",
  DAY_RNG, "); arrows show the order of days. Descriptive cohort-level comparison; no correlation computed."))

## FIGURE 7
fig_cum_traj <- traj_plot(cum_d, "Mean cumulative progeny per worm")
save_fig(fig_cum_traj, "Figure_7_cumulative_progeny_vs_predicted_age")
register_legend(7, "cumulative progeny vs predicted age", paste0(
  "Mean cumulative progeny versus mean clock-predicted age (Batch 1). Each point is one chronological day (",
  DAY_RNG, "); arrows show the order of days. Descriptive cohort-level comparison; no correlation computed."))

## ============================================================
## FIGURE 8: Single panel, dual axis 
## ============================================================
sem_df <- function(d) { d$se <- d$sd / sqrt(d$n); d }
pct_s <- sem_df(pct_d); ck_s <- sem_df(ck_d)
SC  <- 100 / Y_CLK
off <- 0.12

## ============================================================
## FIGURE 9: Change in mean clock-predicted age (forest plot)
## Difference in means (later - earlier day); 95% percentile bootstrap CI
## ============================================================
N_BOOT  <- 5000
BOX_LOW <- FALSE    
set.seed(19)

delta_row <- function(a, b) {
  x <- clk$pred[clk$day == a]            # earlier day 
  y <- clk$pred[clk$day == b]            # later day
  est <- mean(y) - mean(x)
  bs  <- replicate(N_BOOT, mean(sample(y, replace = TRUE)) -
                     mean(sample(x, replace = TRUE)))
  ci  <- quantile(bs, c(0.025, 0.975), names = FALSE)
  mw  <- suppressWarnings(wilcox.test(y, x, alternative = "greater",
                                      exact = FALSE))
  data.frame(from = a, to = b, label = sprintf("D%d to D%d", a, b),
             est = est, lo = ci[1], hi = ci[2], p = mw$p.value,
             U = unname(mw$statistic), n_from = length(x), n_to = length(y))
}
fmt_lab <- function(est, p)
  sprintf("%+.2f d (p %s)", est, ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))

d_consec  <- do.call(rbind, lapply(3:5, function(a) delta_row(a, a + 1)))
d_overall <- delta_row(4, 6)
all_d <- rbind(d_consec, d_overall)
write.csv(all_d, file.path(PROG_OUT, "B1_clock_predicted_age_changes.csv"), row.names = FALSE)
cat("\n--- Change in mean predicted age (bootstrap 95% CI; one-sided MWU) ---\n")
print(all_d, row.names = FALSE, digits = 3)
cat(sprintf("Overall D4->D6: %+.2f d, 95%% CI %.2f-%.2f\n",
            d_overall$est, d_overall$lo, d_overall$hi))

prep <- function(d) {
  d$txt   <- fmt_lab(d$est, d$p)
  d$label <- factor(d$label, levels = rev(d$label)); d
}
d_consec <- prep(d_consec); d_overall <- prep(d_overall)
X_LIM_D  <- c(-0.5, max(all_d$hi) + 1.4)

delta_panel <- function(df, header, show_x = FALSE) {
  p <- ggplot(df, aes(est, label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = LW) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0.18, colour = cK, linewidth = 0.8) +
    geom_point(shape = 17, size = 3.5, colour = cK) +
    geom_text(aes(x = hi + 0.06, label = txt), hjust = 0,
              family = FIG_FONT, fontface = "bold", size = ANN) +
    scale_x_continuous(breaks = 0:3) +
    coord_cartesian(xlim = X_LIM_D, clip = "off") +
    labs(title = header, x = NULL, y = NULL) +
    theme(plot.title = element_text(family = FIG_FONT, size = FS, face = "bold",
                                    hjust = 0, colour = "grey30"))
  if (!show_x)
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

p_top <- delta_panel(d_consec, "Consecutive days")
p_bot <- delta_panel(d_overall, "Overall change, day 4 to day 6", show_x = TRUE) +
  labs(x = "Change in mean predicted age (days)")
if (BOX_LOW)
  p_bot <- p_bot + theme(
    plot.background = element_rect(colour = "black", fill = NA, linewidth = LW),
    plot.margin = margin(4, 8, 4, 4))

fig_delta <- (p_top / p_bot) + plot_layout(heights = c(3, 1.3))
save_fig(fig_delta, "Figure_9_change_in_mean_predicted_age")
register_legend(9, "change in clock-predicted age", sprintf(paste0(
  "Change in mean predicted age between days for the days 1-6 clock, calculated as later day minus earlier day. ",
  "Points show the difference in means, and bars show 95%% bootstrap confidence intervals (%s resamples of animals ",
  "within each day). The upper panel shows consecutive days, and the lower panel shows the overall change from day 4 ",
  "to day 6 (%+.2f d, 95%% CI %.2f-%.2f). P values are from one-sided Mann-Whitney U tests."),
  format(N_BOOT, big.mark = ","), d_overall$est, d_overall$lo, d_overall$hi))
print(fig_delta)

fig_cum_single <- ggplot() +
  annotate("rect", xmin = plateau_start, xmax = XLIM[2],
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35) +
  geom_jitter(data = cumpct, aes(day - off, pct), width = 0.04, height = 0,
              size = PT_SZ, alpha = 0.35, colour = "grey50") +
  geom_jitter(data = clk, aes(day + off, pred * SC), width = 0.04, height = 0,
              size = PT_SZ, alpha = 0.4, colour = "lightsteelblue3") +
  geom_errorbar(data = pct_s, aes(day - off, ymin = m - se, ymax = m + se), width = 0.08, colour = cP, linewidth = LW) +
  geom_line(data = pct_s,  aes(day - off, m), colour = cP, linewidth = 1) +
  geom_point(data = pct_s, aes(day - off, m), colour = cP, size = 3) +
  geom_errorbar(data = ck_s, aes(day + off, ymin = (m - se) * SC, ymax = (m + se) * SC), width = 0.08, colour = cK, linewidth = LW) +
  geom_line(data = ck_s,  aes(day + off, m * SC), colour = cK, linewidth = 1, linetype = "dashed") +
  geom_point(data = ck_s, aes(day + off, m * SC), colour = cK, size = 3, shape = 17) +
  scale_x_continuous(breaks = ALL_DAYS) +
  scale_y_continuous(name = "Cumulative progeny (% of brood)", breaks = seq(0, 100, 25),
                     sec.axis = sec_axis(~ . / SC, breaks = 0:Y_CLK, name = "Clock-predicted age (days)")) +
  coord_cartesian(xlim = XLIM, ylim = c(0, 105)) +
  labs(x = "Chronological age (days)") +
  theme(axis.title.y.left  = element_text(colour = cP), axis.text.y.left  = element_text(colour = cP),
        axis.title.y.right = element_text(colour = cK), axis.text.y.right = element_text(colour = cK))
save_fig(fig_cum_single, "Figure_8_cumulative_progeny_and_predicted_age_single_panel")
register_legend(8, "cumulative progeny and predicted age, single panel", paste0(
  "Cumulative progeny and clock-predicted age (Batch 1). Red (left axis): cumulative progeny as a percentage of each ",
  "worm's lifetime brood; blue (right axis): clock-predicted age (leave-one-out cross-validation). Points: individual ",
  "animals; lines: mean \u00b1 SEM. The right axis is scaled for display only. Shaded region: Day ", plateau_start,
  " onward (mean cumulative output \u2265 ", PLATEAU_PCT, "% of brood). ", SEP_TXT))

print(fig_daily); print(fig_cum); print(fig_clock)
print(fig_daily_stack); print(fig_cum_stack)
print(fig_daily_traj); print(fig_cum_traj); print(fig_cum_single)
cat("\nSaved 8 figures (.tiff/.pdf) and Figure_legends.txt to:\n", PROG_OUT, "\n")