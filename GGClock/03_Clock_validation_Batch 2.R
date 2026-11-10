## ====
## BATCH 2 — EXTERNAL VALIDATION OF A BATCH-1 GG CLOCK
## ====

options(warn = 1)
library(ggplot2)

MAIN_DAY  <- 4                 
FIG_FONT  <- "Times New Roman"

## ====
## A: DEPENDENCY CHECK
## ====

required_functions <- c(
  "run_pair", "extract_features_real", "write_feature_tables",
  "clock_matrix", "clock_apply", "fit_clock_missing",
  "predict_clock_missing", "finite_mean", "finite_median",
  "finite_sd", "safe_cor", "validate_noise_reference",
  "clean_numeric_nonfinite", "save_fig", "register_legend"
)
missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1),
          mode = "function", inherits = TRUE)
]

required_objects <- c(
  "feature_cols", "BLEED_COEFF", "LOCAL_CONTRAST_Z_MIN",
  "pub_theme", "trait_bold", "feature_channel", "feat_display",
  "chan_pal", "col_day", "CLOCK_ALPHA", "CLOCK_SEED",
  "PT_SZ", "LW", "LW_AXIS", "SYM_RHO", "SYM_R2"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(c(missing_functions, missing_objects))) {
  stop("Load the Batch 1 definitions first. Missing: ",
       paste(c(missing_functions, missing_objects), collapse = ", "))
}
stopifnot(CLOCK_ALPHA == 0.5)
feat_display[] <- sub("area \\(px\\)", "area (px\u00B2)", feat_display)

## ====
## A2: Times New Roman 12 pt, all text bold
## ====

FS_B2  <- 12
ANN_B2 <- FS_B2 / .pt          # in-plot geom_text/annotate size (mm) = 12 pt

pub_theme <- pub_theme +
  theme(
    text        = element_text(family = FIG_FONT, size = FS_B2,
                               face = "bold", colour = "black"),
    axis.title  = element_text(family = FIG_FONT, size = FS_B2,
                               face = "bold", colour = "black"),
    axis.text   = element_text(family = FIG_FONT, size = FS_B2,
                               face = "bold", colour = "black"),
    legend.text = element_text(family = FIG_FONT, size = FS_B2,
                               face = "bold", colour = "black"),
    strip.text  = element_text(family = FIG_FONT, size = FS_B2,
                               face = "bold", colour = "black"),
    axis.line   = element_line(linewidth = LW_AXIS, colour = "black"),
    axis.ticks  = element_line(linewidth = LW_AXIS, colour = "black")
  )


trait_bold <- theme(
  axis.text.y = element_text(family = FIG_FONT, size = FS_B2,
                             face = "bold", colour = "black")
)

theme_set(pub_theme)
update_geom_defaults("text", list(family = FIG_FONT, size = ANN_B2,
                                  colour = "black", fontface = "bold"))

## ====
## B: PATHS
## ====

B2_ROOT <- "/Users/purbasha/Desktop/Germ granules/Age clock/Analysis_raw image"

B2_OUT_DIR <- file.path(
  B2_ROOT, "results_v5g_B2_real",
  "Batch2_results_v5g_ALLDAYS", "external_validation_B1"
)
B2_FIG_DIR <- file.path(B2_OUT_DIR, "figures")
B2_OVL_DIR <- file.path(B2_OUT_DIR, "qc_overlays")
B1_OUT_DIR <- file.path(
  B2_ROOT, "results_v5g_B1_real", "Batch1_results_v5g_ALLDAYS"
)

for (d in c(B2_OUT_DIR, B2_FIG_DIR, B2_OVL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

B2_DAY_REGEX <- "^B2_D[0-9]+$"
AL_REGEX <- "^AL_B2_D[0-9]+_W[0-9]+$"
IF_REGEX <- "^IF_B2_D[0-9]+_W[0-9]+$"

## ---- figure numbering: Figures 1-3 are reserved for the main figures;
##      supplementary figures start at Figure_4 ----
.FIG <- new.env()
.FIG$n <- 3L
next_fig <- function() {
  .FIG$n <- .FIG$n + 1L
  .FIG$n
}

plab <- function(p) {
  if (!is.finite(p)) "P unavailable"
  else if (p < 0.0001) "P < 0.0001"
  else sprintf("P = %.3f", p)
}

## typographic minus for legend numbers
mn <- function(x, d = 2) sub("-", "\u2212", sprintf(paste0("%.", d, "f"), x))

## ====
## C: INPUT 
## ====

read_external_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

check_worm_table <- function(d, label) {
  needed <- c("worm_id", "condition", "day", "n_nuclei", feature_cols)
  if (!all(needed %in% names(d))) {
    stop(label, " is missing: ",
         paste(setdiff(needed, names(d)), collapse = ", "))
  }
  stopifnot(
    nrow(d) > 0L, !anyNA(d$worm_id), !anyDuplicated(d$worm_id),
    is.numeric(d$day), all(is.finite(d$day))
  )
  for (f in feature_cols) {
    if (is.logical(d[[f]]) && all(is.na(d[[f]]))) d[[f]] <- as.numeric(d[[f]])
    if (!is.numeric(d[[f]])) stop(label, ": non-numeric feature ", f)
    d[[f]][!is.finite(d[[f]])] <- NA_real_
  }
  d
}

## ====
## D: FILE DISCOVERY — B2 NAMING CONVENTIONS
## ====

find_nuclei <- function(folder, ch, pat) {
  files <- list.files(folder, pattern = "\\.(tif|tiff)$",
                      full.names = TRUE, ignore.case = TRUE,
                      recursive = FALSE)
  files <- files[!dir.exists(files)]
  if (!length(files)) {
    return(data.frame(nucleus = numeric(0), path = character(0)))
  }
  ids <- vapply(basename(files), function(nm) {
    m <- regmatches(nm, regexec(pat, nm, ignore.case = TRUE))[[1]]
    if (length(m) < 2L) NA_real_ else as.numeric(m[2])
  }, numeric(1))
  ok <- is.finite(ids) & ids >= 1
  if (anyDuplicated(ids[ok])) stop("Duplicate ", ch, " nucleus IDs in ", folder)
  data.frame(nucleus = ids[ok], path = files[ok], stringsAsFactors = FALSE)
}

find_nuclei_AL <- function(folder, ch) {
  find_nuclei(folder, ch, paste0("_", ch, "n?([0-9]+)\\.(tif|tiff)$"))
}
find_nuclei_IF <- function(folder, ch) {
  find_nuclei(folder, ch, paste0(ch, "n([0-9]+)\\.(tif|tiff)$"))
}

get_day <- function(wf) {
  h <- regmatches(basename(wf), regexpr("D([0-9]+)", basename(wf)))
  if (!length(h)) NA_integer_ else as.integer(sub("D", "", h))
}
get_worm <- function(wf) {
  h <- regmatches(basename(wf), regexpr("W([0-9]+)$", basename(wf)))
  if (!length(h)) NA_character_ else sub("W", "", h)
}

collect_pairs_B2 <- function(wf, cond) {
  stopifnot(cond %in% c("AL", "IF"))
  day <- get_day(wf)
  wn <- get_worm(wf)
  if (is.na(day) || is.na(wn)) return(NULL)
  
  finder <- if (cond == "AL") find_nuclei_AL else find_nuclei_IF
  red <- finder(wf, "R")
  grn <- finder(wf, "G")
  ids <- intersect(red$nucleus, grn$nucleus)
  if (!length(ids)) {
    warning("No matched red/green nucleus IDs in ", wf)
    return(NULL)
  }
  data.frame(
    worm_id = paste0(cond, "_D", day, "_W", wn),
    condition = cond, day = as.integer(day), worm = wn, nucleus = ids,
    base_id = paste0(cond, "_B2_D", day, "_W", wn, "_n", ids),
    red_path = red$path[match(ids, red$nucleus)],
    green_path = grn$path[match(ids, grn$nucleus)],
    stringsAsFactors = FALSE
  )
}

collect_group_pairs <- function(root, cond, folder_regex) {
  day_dirs <- list.dirs(root, recursive = FALSE)
  day_dirs <- day_dirs[grepl(B2_DAY_REGEX, basename(day_dirs))]
  if (!length(day_dirs)) {
    stop("No folders matching ", B2_DAY_REGEX, " under ", root)
  }
  out <- list()
  for (dd in day_dirs) {
    wfs <- list.dirs(dd, recursive = FALSE)
    wfs <- wfs[grepl(folder_regex, basename(wfs))]
    for (wf in wfs) {
      pr <- collect_pairs_B2(wf, cond)
      if (!is.null(pr)) out[[length(out) + 1L]] <- pr
    }
  }
  tb <- do.call(rbind, out)
  if (is.null(tb) || !nrow(tb)) {
    cat("[", cond, "] no matched nucleus pairs found.\n")
    return(NULL)
  }
  if (anyDuplicated(tb$base_id)) stop("Duplicate base_id in ", cond)
  cat(sprintf("[%s] %d nucleus pairs across %d worms; days: %s\n",
              cond, nrow(tb), length(unique(tb$worm_id)),
              paste(sort(unique(tb$day)), collapse = ",")))
  tb
}

## ====
## E: B1 TRAINING INPUTS AND SHARED NOISE REFERENCE
## ====

B1_per_worm <- check_worm_table(
  read_external_csv(file.path(B1_OUT_DIR, "Batch1_per_worm.csv")),
  "B1 per-worm table"
)

stopifnot(
  nrow(B1_per_worm) >= 4L,
  length(unique(B1_per_worm$day)) >= 2L,
  !anyNA(B1_per_worm$condition),
  all(tolower(trimws(B1_per_worm$condition)) %in% c("al", "adlib"))
)

no_foci_B1 <- !is.na(B1_per_worm$green_n) & B1_per_worm$green_n == 0
if (any(!is.na(as.matrix(B1_per_worm[
  no_foci_B1, c("green_size", "green_elong", "green_app_net"), drop = FALSE
])))) {
  stop("B1 contains unmasked no-foci traits. ",
       "Regenerate the conditional B1 table from image-level data.")
}

b1_ref_check <- read_external_csv(
  file.path(B1_OUT_DIR, "Batch1_noise_reference.csv")
)
if (!"batch_ref" %in% names(b1_ref_check)) {
  stop("Batch1_noise_reference.csv is missing column: batch_ref")
}
if (nrow(b1_ref_check) != 1L) {
  stop("Batch1_noise_reference.csv must contain exactly one reference row.")
}
B1_REF <- suppressWarnings(as.numeric(as.character(b1_ref_check$batch_ref)))
if (length(B1_REF) != 1L || !is.finite(B1_REF) || B1_REF <= 0) {
  stop("The saved Batch 1 noise reference must be finite and positive.")
}
validate_noise_reference(B1_REF, "B1 AL noise reference")

write.csv(data.frame(B1_REF = B1_REF),
          file.path(B2_OUT_DIR, "Batch1_reference_used.csv"),
          row.names = FALSE)
cat(sprintf("[REFERENCE] B1 AL local_bg_mad = %.6g\n", B1_REF))

## ====
## F: B2 EXTRACTION USING THE B1 NOISE REFERENCE
## ====

process_group_B2 <- function(pair_tbl, tag) {
  stopifnot(!is.null(pair_tbl), nrow(pair_tbl) > 0L,
            !anyDuplicated(pair_tbl$base_id))
  
  mad_i <- vapply(seq_len(nrow(pair_tbl)), function(i) {
    tryCatch({
      p <- run_pair(pair_tbl$red_path[i], pair_tbl$green_path[i],
                    bleed_coeff = BLEED_COEFF,
                    contrast_z_min = LOCAL_CONTRAST_Z_MIN)
      finite_median(p$green_diagnostic$table$local_bg_mad)
    }, error = function(e) {
      warning("Noise measurement failed: ", pair_tbl$base_id[i],
              ": ", e$message)
      NA_real_
    })
  }, numeric(1))
  
  pair_tbl$img_bg_mad <- mad_i
  pair_tbl$noise_ratio <- mad_i / B1_REF
  pair_tbl$noise_fallback <- !is.finite(pair_tbl$noise_ratio)
  pair_tbl$noise_ratio[pair_tbl$noise_fallback] <- 1.0
  pair_tbl$adaptive_z_min <- LOCAL_CONTRAST_Z_MIN * pair_tbl$noise_ratio
  
  if (any(pair_tbl$noise_fallback)) {
    warning(sum(pair_tbl$noise_fallback), " ", tag,
            " images lack a noise estimate; using noise_ratio = 1.")
  }
  
  write.csv(pair_tbl,
            file.path(B2_OUT_DIR, paste0("Batch2_", tag, "_noise_and_pairs.csv")),
            row.names = FALSE)
  cat(sprintf("[%s] Own median MAD=%.6g (diagnostic only); B1 reference=%.6g\n",
              tag, finite_median(mad_i), B1_REF))
  
  group_overlay <- file.path(B2_OVL_DIR, tag)
  dir.create(group_overlay, recursive = TRUE, showWarnings = FALSE)
  
  feats <- lapply(seq_len(nrow(pair_tbl)), function(i) {
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf("  [%s %d/%d] %s\n", tag, i, nrow(pair_tbl),
                  pair_tbl$base_id[i]))
    }
    meta <- as.list(pair_tbl[i, c("condition", "day", "worm_id",
                                  "nucleus", "base_id")])
    tryCatch(
      extract_features_real(
        pair_tbl$red_path[i], pair_tbl$green_path[i], meta,
        overlay_dir = group_overlay,
        contrast_z_min = pair_tbl$adaptive_z_min[i]
      ),
      error = function(e) {
        warning("SKIP ", pair_tbl$base_id[i], ": ", e$message)
        NULL
      }
    )
  })
  
  success <- vapply(feats, function(x) !is.null(x) && nrow(x) > 0L,
                    logical(1))
  audit <- pair_tbl[, c("base_id", "worm_id", "day", "nucleus"), drop = FALSE]
  audit$extraction_success <- success
  write.csv(audit,
            file.path(B2_OUT_DIR, paste0("Batch2_", tag, "_extraction_audit.csv")),
            row.names = FALSE)
  
  if (!any(success)) stop("No successfully extracted images for ", tag)
  if (any(!success)) {
    warning(sum(!success), " failed ", tag, " images; see extraction audit.")
  }
  
  per_image <- do.call(rbind, feats[success])
  
  ## Saves raw zerosIncluded tables, masks green_size/green_elong/
  ## green_app_net when green_n == 0, then aggregates per worm.
  tables <- write_feature_tables(per_image, B2_OUT_DIR,
                                 prefix = paste0("Batch2_", tag))
  tables$per_worm
}

## ====
## G: FIT B1 ONLY; APPLY TO BOTH B2 GROUPS
## ====

train_B1_clock <- function(B1_pw, seed = 19) {
  set.seed(seed)
  model <- fit_clock_missing(clock_matrix(B1_pw), as.numeric(B1_pw$day))
  
  if (model$boundary) {
    warning("B1 tuning selected a lambda-grid boundary. ",
            "Inspect the tuning results; do not tune using B2 performance.")
  }
  
  model$training_worm_ids <- B1_pw$worm_id
  model$training_days <- B1_pw$day
  model$noise_reference <- B1_REF
  model$training_source <- file.path(B1_OUT_DIR, "Batch1_per_worm.csv")
  model$seed <- seed
  
  saveRDS(model, file.path(B2_OUT_DIR,
                           "B1_only_clock_used_for_external_validation.rds"))
  write.csv(
    data.frame(lambda = model$lambda, CV_MSE = model$cvm,
               selected = model$lambda == model$lambda.min),
    file.path(B2_OUT_DIR, "B1_clock_tuning.csv"), row.names = FALSE
  )
  write.csv(
    data.frame(
      trait = feature_cols,
      imputation_median = unname(model$prep$med[feature_cols]),
      disabled = unname(model$prep$disabled[feature_cols]),
      coefficient = unname(model$coefficient[feature_cols])
    ),
    file.path(B2_OUT_DIR, "B1_clock_preprocessing_and_coefficients.csv"),
    row.names = FALSE
  )
  model
}

project_B2 <- function(clock_obj, pw, label) {
  pw <- check_worm_table(pw, paste("B2", label))
  predicted <- predict_clock_missing(clock_obj, clock_matrix(pw))
  stopifnot(length(predicted) == nrow(pw), all(is.finite(predicted)))
  data.frame(worm_id = pw$worm_id, group = label, day = pw$day,
             n_nuclei = pw$n_nuclei, predicted = predicted,
             stringsAsFactors = FALSE)
}

## ====
## H1: B2 AL EXTERNAL VALIDATION  
## ====

r2_sse <- function(y, p) {
  sst <- sum((y - mean(y))^2)
  if (sst == 0) NA_real_ else 1 - sum((y - p)^2) / sst
}

evaluate_B2_AL <- function(al_tbl) {
  stopifnot(nrow(al_tbl) > 0L, all(is.finite(al_tbl$predicted)),
            all(is.finite(al_tbl$day)))
  
  al_tbl$error <- al_tbl$predicted - al_tbl$day
  
  metrics <- data.frame(
    n = nrow(al_tbl),
    rho = safe_cor(al_tbl$day, al_tbl$predicted, "spearman"),
    r2_pearson = safe_cor(al_tbl$day, al_tbl$predicted)^2,
    R2_SSE = r2_sse(al_tbl$day, al_tbl$predicted),
    MAE = mean(abs(al_tbl$error)),
    RMSE = sqrt(mean(al_tbl$error^2)),
    bias = mean(al_tbl$error)
  )
  
  cat("\n==== B2 AL: EXTERNAL VALIDATION ====\n")
  print(metrics, digits = 4)
  write.csv(metrics, file.path(B2_OUT_DIR, "Batch2_AL_external_metrics.csv"),
            row.names = FALSE)
  
  by_day <- aggregate(
    cbind(MAE = abs(error), RMSE = error^2, bias = error) ~ day,
    data = al_tbl, FUN = mean
  )
  by_day$RMSE <- sqrt(by_day$RMSE)
  write.csv(by_day, file.path(B2_OUT_DIR, "Batch2_AL_external_MAE_by_day.csv"),
            row.names = FALSE)
  print(by_day, digits = 3)
  
  dl <- sort(unique(al_tbl$day))
  al_tbl$day_f <- factor(al_tbl$day, levels = dl)
  
  use_day_cols <- all(as.character(dl) %in% names(col_day))
  pt_layer <- if (use_day_cols) {
    list(geom_jitter(aes(colour = day_f), width = 0.08, height = 0,
                     size = PT_SZ, alpha = 0.85, stroke = 0.3),
         scale_colour_manual(values = col_day))
  } else {
    list(geom_jitter(width = 0.08, height = 0, size = PT_SZ, alpha = 0.85,
                     stroke = 0.3, colour = "grey30"))
  }
  
  ## in-plot clock metrics (Figure 1 only)
  stat_txt <- paste0(
    "n = ", metrics$n, "\n",
    SYM_RHO, " = ", sprintf("%.2f", metrics$rho), "\n",
    "r\u00b2 = ", sprintf("%.2f", metrics$r2_pearson), "\n",
    "MAE = ", sprintf("%.2f", metrics$MAE), " d\n",
    "RMSE = ", sprintf("%.2f", metrics$RMSE), " d"
  )
  
  ## coord_cartesian (not scale limits) so no points are dropped.
  ## If the text overlaps points, use x = Inf, y = -Inf, hjust = 1.08,
  ## vjust = -0.15 (bottom-right corner).
  p_ext <- ggplot(al_tbl, aes(day, predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey50", linewidth = LW) +
    pt_layer +
    annotate("text", x = -Inf, y = Inf, label = stat_txt,
             hjust = -0.08, vjust = 1.15, family = FIG_FONT,
             fontface = "bold", size = ANN_B2, lineheight = 0.95) +
    scale_x_continuous(breaks = dl, labels = paste0("D", dl)) +
    scale_y_continuous(breaks = dl, labels = paste0("D", dl)) +
    coord_cartesian(
      xlim = range(dl) + c(-0.5, 0.5),
      ylim = range(c(dl, al_tbl$predicted)) + c(-0.5, 0.5)
    ) +
    labs(x = "Chronological age (days)", y = "Predicted age (days)") +
    pub_theme +
    theme(legend.position = "none")
  
  invisible(list(al = al_tbl, metrics = metrics, by_day = by_day,
                 plot = p_ext))
}

## ====
## H2: IF vs AL CLOCK PROJECTION
## ====

project_IF_vs_AL_B2 <- function(al_tbl, if_tbl,
                                out_dir = B2_OUT_DIR,
                                fig_dir = B2_FIG_DIR,
                                tag = "_B1ref") {
  if (is.null(if_tbl) || nrow(if_tbl) == 0L) {
    cat("\n[IF] Skipping comparison: no B2 IF data.\n")
    return(NULL)
  }
  
  G_AL <- "AL (external)"
  G_IF <- "IF (projected)"
  
  al2 <- data.frame(worm_id = al_tbl$worm_id, group = G_AL,
                    day = al_tbl$day, predicted = al_tbl$predicted,
                    stringsAsFactors = FALSE)
  if2 <- data.frame(worm_id = if_tbl$worm_id, group = G_IF,
                    day = if_tbl$day, predicted = if_tbl$predicted,
                    stringsAsFactors = FALSE)
  
  stopifnot(all(is.finite(al2$predicted)), all(is.finite(if2$predicted)),
            all(is.finite(if2$day)))
  
  if_days_all <- sort(unique(if2$day))
  if_days <- if_days_all[if_days_all %in% al2$day]
  if (length(if_days) < length(if_days_all)) {
    warning("IF day(s) without same-day B2 AL controls: ",
            paste(setdiff(if_days_all, if_days), collapse = ", "),
            ". Kept in the projection table; excluded from group inference.")
  }
  if (!length(if_days)) stop("No shared AL/IF days in B2.")
  
  al_err <- al2$predicted - al2$day
  al_bias <- tapply(al_err, al2$day, mean)
  raw_accel <- if2$predicted - if2$day
  same_day_bias <- as.numeric(al_bias[as.character(if2$day)])
  
  if_out <- data.frame(
    worm_id = if_tbl$worm_id, n_nuclei = if_tbl$n_nuclei,
    true_day = if2$day, predicted_age = if2$predicted,
    age_accel_raw = raw_accel, al_bias_same_day = same_day_bias,
    age_accel_corrected = raw_accel - same_day_bias,
    stringsAsFactors = FALSE
  )
  write.csv(if_out,
            file.path(out_dir, paste0("Batch2_IF_projected_onto_B1_clock",
                                      tag, ".csv")),
            row.names = FALSE)
  cat("\n==== B2 IF PROJECTED ONTO B1 CLOCK", tag, "====\n")
  print(if_out, row.names = FALSE, digits = 4)
  
  comb <- rbind(al2[al2$day %in% if_days, , drop = FALSE],
                if2[if2$day %in% if_days, , drop = FALSE])
  comb$error <- comb$predicted - comb$day
  comb$group <- factor(comb$group, levels = c(G_AL, G_IF))
  write.csv(comb,
            file.path(out_dir, paste0("Batch2_IF_vs_AL_combined", tag, ".csv")),
            row.names = FALSE)
  
  per_day <- do.call(rbind, lapply(if_days, function(d) {
    a <- comb$error[comb$group == G_AL & comb$day == d]
    f <- comb$error[comb$group == G_IF & comb$day == d]
    ok <- length(a) >= 2L && length(f) >= 2L
    pw <- function(alt) {
      p <- if (ok) {
        suppressWarnings(
          wilcox.test(f, a, alternative = alt, exact = FALSE)$p.value)
      } else NA_real_
      if (!is.finite(p)) NA_real_ else p
    }
    data.frame(day = d, n_AL = length(a), n_IF = length(f),
               mean_err_AL = mean(a), mean_err_IF = mean(f),
               IF_minus_AL = mean(f) - mean(a),
               p_one_sided = pw("less"), p_two_sided = pw("two.sided"))
  }))
  per_day$p_one_sided_BH <- p.adjust(per_day$p_one_sided, method = "BH")
  
  write.csv(per_day,
            file.path(out_dir, paste0("Batch2_IF_vs_AL_by_day", tag, ".csv")),
            row.names = FALSE)
  cat("\n--- Same-day prediction errors (B2) ---\n")
  print(per_day, row.names = FALSE, digits = 3)
  
  fit <- if (length(unique(comb$day)) > 1L) {
    lm(error ~ factor(day) + group, data = comb)
  } else {
    lm(error ~ group, data = comb)
  }
  cf <- summary(fit)$coefficients
  rn <- grep("^group", rownames(cf), value = TRUE)
  
  est <- p_two <- p_less <- NA_real_
  ci <- c(NA_real_, NA_real_)
  if (length(rn) == 1L) {
    est <- unname(cf[rn, 1])
    tval <- unname(cf[rn, 3])
    p_two <- unname(cf[rn, 4])
    if (fit$df.residual > 0L && !is.na(tval)) {
      p_less <- pt(tval, fit$df.residual)
      ci <- as.numeric(suppressWarnings(confint(fit)[rn, ]))
    }
  }
  if (!is.finite(p_two)) p_two <- NA_real_
  if (!is.finite(p_less)) p_less <- NA_real_
  ci[!is.finite(ci)] <- NA_real_
  
  cat(sprintf(paste0("\n[B2 IF vs AL%s] Adjusted IF effect: %+.3f d ",
                     "(95%% CI %+.3f to %+.3f)\n",
                     "Two-sided p=%.4f; one-sided IF-younger p=%.4f\n"),
              tag, est, ci[1], ci[2], p_two, p_less))
  
  write.csv(
    data.frame(effect_days = est, ci_low = ci[1], ci_high = ci[2],
               p_two_sided = p_two, p_one_sided_IF_younger = p_less,
               n_AL = sum(comb$group == G_AL),
               n_IF = sum(comb$group == G_IF)),
    file.path(out_dir, paste0("Batch2_IF_effect_summary", tag, ".csv")),
    row.names = FALSE
  )
  
  per_day_txt <- paste(
    sprintf("D%d: %+.2f d (BH-adjusted %s; n = %d AL, %d IF)",
            per_day$day, per_day$IF_minus_AL, plab(per_day$p_one_sided_BH),
            per_day$n_AL, per_day$n_IF),
    collapse = "; ")
  
  day_levels <- sort(unique(c(al2$day, if_days)))
  day_labels <- paste0("D", day_levels)
  cols <- setNames(c("#4D4D4D", "#D7191C"), c(G_AL, G_IF))
  
  set.seed(1)
  off <- ifelse(comb$group == G_IF, 0.16, -0.16)
  comb$x_plot <- comb$day + off + runif(nrow(comb), -0.05, 0.05)
  
  sm <- do.call(rbind, lapply(
    split(comb, list(comb$group, comb$day), drop = TRUE),
    function(df) {
      data.frame(group = df$group[1], day = df$day[1],
                 xc = df$day[1] + ifelse(df$group[1] == G_IF, 0.16, -0.16),
                 m = mean(df$predicted), s = finite_sd(df$predicted))
    }
  ))
  
  ## ---- Scatter: predicted vs chronological age ----
  p_sc <- ggplot(comb, aes(x = x_plot, y = predicted, colour = group,
                           shape = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey50", linewidth = LW) +
    geom_point(size = PT_SZ, alpha = 0.7) +
    geom_crossbar(data = sm,
                  aes(x = xc, y = m, ymin = m - s, ymax = m + s,
                      colour = group),
                  inherit.aes = FALSE, width = 0.26, linewidth = LW,
                  fill = NA, na.rm = TRUE) +
    scale_colour_manual(values = cols) +
    scale_shape_manual(values = c(16, 17)) +
    scale_x_continuous(breaks = day_levels, labels = day_labels) +
    scale_y_continuous(breaks = day_levels, labels = day_labels) +
    coord_cartesian(xlim = range(day_levels) + c(-0.6, 0.6),
                    ylim = range(c(day_levels, comb$predicted)) + c(-0.3, 0.3)) +
    labs(x = "Chronological age (days)", y = "Predicted age (days)") +
    pub_theme +
    theme(legend.position = "bottom")
  
  fn1 <- next_fig()
  save_fig(p_sc, fn1, fig_dir, "small", h = 80)
  register_legend(fn1, sprintf(
    paste0("Batch 1 clock applied to independent Batch 2 AL and IF worms by ",
           "chronological age. Points: individual worms; crossbars: mean \u00b1 ",
           "SD. Dashed line: perfect prediction. Same-day IF minus AL mean ",
           "prediction error (one-sided Wilcoxon test, IF younger; ",
           "BH-adjusted across days) - %s. Pooled day-adjusted IF effect = ",
           "%+.2f d (95%% CI %+.2f to %+.2f; %s, one-sided). Inference is ",
           "conditional on the fixed Batch 1 clock."),
    per_day_txt, est, ci[1], ci[2], plab(p_less)), fig_dir)
  
  ## ---- Prediction error by age ----
  comb$day_f <- factor(comb$day, levels = day_levels)
  
  p_err <- ggplot(comb, aes(x = day_f, y = error, fill = group,
                            colour = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = LW) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.3,
                 linewidth = LW, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                               dodge.width = 0.75),
               size = PT_SZ, alpha = 0.75) +
    scale_x_discrete(labels = setNames(day_labels, as.character(day_levels))) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    labs(x = "Chronological age (days)", y = "Prediction error (days)") +
    pub_theme +
    theme(legend.position = "bottom")
  
  fn2 <- next_fig()
  save_fig(p_err, fn2, fig_dir, "small", h = 80)
  register_legend(fn2, sprintf(
    paste0("Prediction error (predicted minus chronological age) for Batch 2 ",
           "AL and IF worms by chronological age. Boxes: median and ",
           "interquartile range; whiskers: 1.5 x IQR; points: individual ",
           "worms. Dashed line: zero error. Same-day IF minus AL difference ",
           "(one-sided Wilcoxon test, IF younger; BH-adjusted) - %s. Pooled ",
           "day-adjusted IF effect = %+.2f d (%s, one-sided)."),
    per_day_txt, est, plab(p_less)), fig_dir)
  
  list(if_table = if_out, comparison = comb, per_day = per_day, lm_fit = fit,
       effect = est, ci = ci, p_two = p_two, p_less = p_less)
}

## ====
## H3: CLOCK ATTRIBUTION
## ====

## 
## 
attr_plot <- function(res, show_pct = FALSE, font = FIG_FONT) {
  pl <- res[res$in_clock, , drop = FALSE]
  pl$sig <- is.finite(pl$p_boot_BH) & pl$p_boot_BH < 0.05
  pl$lab <- paste0(pl$display, ifelse(pl$sig, " *", ""))
  ord <- order(pl$contribution_days, decreasing = TRUE)
  pl$lab <- factor(pl$lab, levels = pl$lab[ord])
  pl$channel <- factor(pl$channel, c("RED", "GREEN", "BOTH"),
                       c("PGL-1 (red)", "MUT-16 (green)", "Both channels"))
  
  p <- ggplot(pl, aes(contribution_days, lab, fill = channel)) +
    geom_vline(xintercept = 0, linewidth = LW) +
    geom_col(width = 0.65, colour = "black", linewidth = LW) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                  orientation = "y", width = 0.25, linewidth = LW)
  
  if (show_pct) {
    pl$pct_lab <- sprintf("%.0f%%", pl$share_of_total_pct)
    pl$x_lab <- ifelse(pl$contribution_days >= 0, pl$ci_high, pl$ci_low)
    pl$hj <- ifelse(pl$contribution_days >= 0, -0.2, 1.2)
    p <- p + geom_text(data = pl,
                       aes(x = x_lab, y = lab, label = pct_lab, hjust = hj),
                       family = font, fontface = "bold", size = ANN_B2,
                       colour = "black", inherit.aes = FALSE)
  }
  
  p +
    scale_fill_manual(values = c("PGL-1 (red)" = "#C0392B",
                                 "MUT-16 (green)" = "#3FA66B",
                                 "Both channels" = "#7D3C98"),
                      drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    coord_cartesian(clip = "off") +
    labs(x = "Contribution to clock age shift (days)", y = NULL) +
    pub_theme + trait_bold +
    theme(axis.text.x = element_text(face = "bold", colour = "black"),
          axis.text.y = element_text(face = "bold", colour = "black"),
          plot.margin = margin(5, 8, 5, 12),
          legend.position = "bottom")
}

clock_attribution_B2 <- function(pw_AL, pw_IF, clock_obj, day,
                                 out_dir = B2_OUT_DIR,
                                 fig_dir = B2_FIG_DIR,
                                 tag = "_B1ref",
                                 n_boot = 5000,
                                 focus = c("green_n", "red_radial_cv"),
                                 seed = CLOCK_SEED) {
  al <- pw_AL[pw_AL$day == day, , drop = FALSE]
  fi <- pw_IF[pw_IF$day == day, , drop = FALSE]
  stopifnot(nrow(al) >= 2L, nrow(fi) >= 2L)
  
  XA <- as.matrix(clock_apply(clock_matrix(al), clock_obj$prep))
  XI <- as.matrix(clock_apply(clock_matrix(fi), clock_obj$prep))
  stopifnot(all(feature_cols %in% colnames(XA)),
            all(feature_cols %in% colnames(XI)))
  XA <- XA[, feature_cols, drop = FALSE]
  XI <- XI[, feature_cols, drop = FALSE]
  
  beta <- unname(clock_obj$coefficient[feature_cols])
  beta[!is.finite(beta)] <- 0
  in_clock <- beta != 0
  
  contrib_fun <- function(ia, ii) {
    beta * (colMeans(XI[ii, , drop = FALSE]) - colMeans(XA[ia, , drop = FALSE]))
  }
  
  nA <- nrow(XA); nI <- nrow(XI)
  obs <- contrib_fun(seq_len(nA), seq_len(nI))
  total_obs <- sum(obs)
  
  actual_shift <- mean(predict_clock_missing(clock_obj, clock_matrix(fi))) -
    mean(predict_clock_missing(clock_obj, clock_matrix(al)))
  stopifnot(isTRUE(all.equal(total_obs, actual_shift, tolerance = 1e-8)))
  
  set.seed(seed)
  boot <- t(vapply(seq_len(n_boot), function(b) {
    contrib_fun(sample.int(nA, replace = TRUE), sample.int(nI, replace = TRUE))
  }, numeric(length(feature_cols))))
  colnames(boot) <- feature_cols
  boot_total <- rowSums(boot)
  
  q <- function(x) quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  p_younger <- function(x) (1 + sum(x >= 0)) / (length(x) + 1)
  p_older   <- function(x) (1 + sum(x <= 0)) / (length(x) + 1)
  
  ci_mat <- t(apply(boot, 2, q))
  
  sds <- if (!is.null(clock_obj$prep$sds)) {
    unname(clock_obj$prep$sds[feature_cols])
  } else rep(NA_real_, length(feature_cols))
  
  mean_AL <- vapply(feature_cols, function(f) finite_mean(al[[f]]), numeric(1))
  mean_IF <- vapply(feature_cols, function(f) finite_mean(fi[[f]]), numeric(1))
  shift_raw <- mean_IF - mean_AL
  
  p_boot_younger <- vapply(feature_cols, function(f) p_younger(boot[, f]),
                           numeric(1))
  p_boot_two <- vapply(feature_cols, function(f) {
    min(1, 2 * min(p_younger(boot[, f]), p_older(boot[, f])))
  }, numeric(1))
  p_boot_two[!in_clock] <- NA_real_
  
  res <- data.frame(
    trait = feature_cols,
    display = unname(feat_display[feature_cols]),
    channel = unname(feature_channel[feature_cols]),
    in_clock = in_clock, clock_beta = beta, training_SD = sds,
    mean_AL = unname(mean_AL), mean_IF = unname(mean_IF),
    shift_raw = unname(shift_raw),
    shift_training_SD = unname(shift_raw) / sds,
    contribution_days = obs,
    ci_low = ci_mat[, 1], ci_high = ci_mat[, 2],
    share_of_total_pct = 100 * obs / total_obs,
    p_boot_younger = p_boot_younger,
    p_boot_two = p_boot_two,
    stringsAsFactors = FALSE
  )
  res$p_boot_BH <- NA_real_
  res$p_boot_BH[res$in_clock] <- p.adjust(res$p_boot_two[res$in_clock],
                                          method = "BH")
  res$reliable <- res$in_clock & (res$ci_high < 0 | res$ci_low > 0)
  res$direction <- ifelse(res$contribution_days < 0, "younger",
                          ifelse(res$contribution_days > 0, "older", "none"))
  res <- res[order(res$contribution_days), , drop = FALSE]
  rownames(res) <- NULL
  
  total_ci <- q(boot_total)
  focus <- intersect(focus, feature_cols)
  focus_obs <- sum(obs[feature_cols %in% focus])
  focus_boot <- rowSums(boot[, focus, drop = FALSE])
  focus_ci <- q(focus_boot)
  focus_share_boot <- 100 * focus_boot / boot_total
  focus_share_boot <- focus_share_boot[is.finite(focus_share_boot)]
  focus_share_ci <- q(focus_share_boot)
  
  summ <- data.frame(
    day = day, n_AL = nA, n_IF = nI, n_boot = n_boot,
    total_shift_days = total_obs,
    total_ci_low = total_ci[1], total_ci_high = total_ci[2],
    p_boot_total_younger = p_younger(boot_total),
    focus_traits = paste(focus, collapse = "+"),
    focus_contribution_days = focus_obs,
    focus_ci_low = focus_ci[1], focus_ci_high = focus_ci[2],
    focus_share_pct = 100 * focus_obs / total_obs,
    focus_share_ci_low = focus_share_ci[1],
    focus_share_ci_high = focus_share_ci[2],
    pct_boot_same_sign_total = 100 * mean(sign(focus_boot) == sign(boot_total))
  )
  
  write.csv(res, file.path(out_dir,
                           sprintf("Batch2_IF_vs_AL_clock_attribution_D%d%s.csv", day, tag)),
            row.names = FALSE)
  write.csv(data.frame(boot_id = seq_len(n_boot), boot, total = boot_total),
            file.path(out_dir,
                      sprintf("Batch2_IF_vs_AL_clock_attribution_boot_D%d%s.csv", day, tag)),
            row.names = FALSE)
  write.csv(summ, file.path(out_dir,
                            sprintf("Batch2_IF_vs_AL_clock_attribution_summary_D%d%s.csv", day, tag)),
            row.names = FALSE)
  
  cat(sprintf(paste0(
    "\n[ATTRIBUTION] D%d: n_AL=%d, n_IF=%d; total shift = %+.3f d ",
    "(95%% CI %+.3f to %+.3f; one-sided p = %.4f)\n",
    "[ATTRIBUTION] Focus (%s): %+.3f d (95%% CI %+.3f to %+.3f); ",
    "%.1f%% of total (95%% CI %.1f to %.1f)\n",
    "[ATTRIBUTION] Clock weights fixed; B1 training uncertainty not propagated.\n"),
    day, nA, nI, total_obs, total_ci[1], total_ci[2], summ$p_boot_total_younger,
    paste(focus, collapse = " + "), focus_obs, focus_ci[1], focus_ci[2],
    summ$focus_share_pct, focus_share_ci[1], focus_share_ci[2]))
  print(res[, c("display", "contribution_days", "ci_low", "ci_high",
                "share_of_total_pct", "p_boot_BH", "reliable")],
        row.names = FALSE, digits = 3)
  
  ## ---- Supplementary figure ----
  n_zero <- sum(!res$in_clock)
  p_att <- attr_plot(res, show_pct = TRUE)
  
  fn <- next_fig()
  save_fig(p_att, fn, fig_dir, "large", h = 130)
  register_legend(fn, sprintf(
    paste0("Trait contributions to the difference in clock age between Batch 2 ",
           "IF and AL worms at day %d (n = %d AL, %d IF worms). Contribution = ",
           "Batch 1 clock weight x (mean IF - mean AL trait value); negative ",
           "values indicate a younger clock age in IF. Total shift = %+.2f d ",
           "(95%% bootstrap CI %+.2f to %+.2f; one-sided bootstrap P = %.4f). ",
           "Focus traits (%s) contribute %+.2f d (95%% CI %+.2f to %+.2f), ",
           "i.e. %.1f%% of the total shift (95%% CI %.1f to %.1f%%). Whiskers: ",
           "95%% bootstrap CI over worms (n = %s resamples; AL and IF ",
           "resampled separately); Batch 1 clock weights held fixed and ",
           "training uncertainty not propagated. Asterisk: BH-adjusted ",
           "two-sided bootstrap P < 0.05. Percentages: share of the total shift.%s"),
    day, nA, nI, total_obs, total_ci[1], total_ci[2],
    summ$p_boot_total_younger,
    paste(unname(feat_display[focus]), collapse = " + "),
    focus_obs, focus_ci[1], focus_ci[2], summ$focus_share_pct,
    focus_share_ci[1], focus_share_ci[2], format(n_boot, big.mark = ","),
    if (n_zero > 0L) sprintf(" %d zero-weight trait(s) omitted.", n_zero)
    else ""), fig_dir)
  
  invisible(list(table = res, summary = summ, boot = boot,
                 boot_total = boot_total, fig_num = fn, n_zero = n_zero))
}

## ====
## I: TRAIT TESTS AND EXACT CLOCK CONTRIBUTIONS
## ====

compare_traits_B2 <- function(pw_AL, pw_IF, clock_obj, day) {
  al <- pw_AL[pw_AL$day == day, , drop = FALSE]
  fi <- pw_IF[pw_IF$day == day, , drop = FALSE]
  stopifnot(nrow(al) > 0L, nrow(fi) > 0L)
  
  cliff <- function(x, y) {
    if (!length(x) || !length(y)) return(NA_real_)
    mean(sign(outer(x, y, "-")))
  }
  
  res <- do.call(rbind, lapply(feature_cols, function(f) {
    a <- al[[f]]; a <- a[is.finite(a)]
    x <- fi[[f]]; x <- x[is.finite(x)]
    delta <- finite_mean(x) - finite_mean(a)
    sa <- finite_sd(a)
    p <- if (length(a) >= 3L && length(x) >= 3L) {
      suppressWarnings(wilcox.test(x, a, exact = FALSE)$p.value)
    } else NA_real_
    data.frame(trait = f, n_AL = length(a), n_IF = length(x),
               median_AL = finite_median(a), median_IF = finite_median(x),
               diff_mean = delta,
               std_diff = if (is.finite(sa) && sa > 0) delta / sa else NA_real_,
               cliff_delta = cliff(x, a), p_wilcox = p)
  }))
  
  res <- clean_numeric_nonfinite(res)
  res$p_BH <- p.adjust(res$p_wilcox, method = "BH")
  
  XA <- clock_apply(clock_matrix(al), clock_obj$prep)
  XI <- clock_apply(clock_matrix(fi), clock_obj$prep)
  delta_clock <- colMeans(XI) - colMeans(XA)
  
  res$clock_coef <- unname(clock_obj$coefficient[res$trait])
  res$diff_mean_clock <- unname(delta_clock[res$trait])
  res$contrib_days <- res$clock_coef * res$diff_mean_clock
  
  actual_shift <- mean(predict_clock_missing(clock_obj, clock_matrix(fi))) -
    mean(predict_clock_missing(clock_obj, clock_matrix(al)))
  stopifnot(isTRUE(all.equal(sum(res$contrib_days), actual_shift,
                             tolerance = 1e-8)))
  
  res <- res[order(res$p_wilcox, na.last = TRUE), , drop = FALSE]
  write.csv(res, file.path(B2_OUT_DIR,
                           sprintf("Batch2_IF_vs_AL_trait_wilcoxon_D%d.csv", day)),
            row.names = FALSE)
  
  cat(sprintf("\n[TRAITS] D%d: n_AL=%d, n_IF=%d; summed contributions=%+.3f d\n",
              day, nrow(al), nrow(fi), sum(res$contrib_days)))
  print(res, row.names = FALSE, digits = 3)
  
  pc <- res
  pc$channel <- unname(feature_channel[pc$trait])
  pc$display <- unname(feat_display[pc$trait])
  pc$sig <- is.finite(pc$p_BH) & pc$p_BH < 0.05
  pc$lab <- paste0(pc$display, ifelse(pc$sig, " *", ""))
  pc <- pc[order(pc$contrib_days), , drop = FALSE]
  pc$lab <- factor(pc$lab, levels = pc$lab)
  
  pB <- ggplot(pc, aes(contrib_days, lab, fill = channel)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = LW) +
    geom_col(width = 0.65, colour = "grey20", linewidth = LW) +
    scale_fill_manual(values = chan_pal) +
    labs(x = "Contribution to clock age shift (days)", y = NULL) +
    pub_theme + trait_bold +
    theme(axis.text.x = element_text(face = "bold", colour = "black"),
          axis.text.y = element_text(face = "bold", colour = "black"),
          plot.margin = margin(5, 8, 5, 5),
          legend.position = "bottom")
  
  fn <- next_fig()
  save_fig(pB, fn, B2_FIG_DIR, "large", h = 130)
  register_legend(fn, sprintf(
    paste0("Clock-weighted trait contributions to the Batch 2 IF-versus-AL ",
           "difference in predicted age at day %d (n = %d AL, %d IF worms). ",
           "Contribution = Batch 1 coefficient x (Batch 1-imputed mean IF - ",
           "mean AL trait value); contributions sum exactly to the raw mean ",
           "prediction shift (%+.3f d). Asterisk: Benjamini-Hochberg-adjusted ",
           "Wilcoxon rank-sum P < 0.05 for the raw trait difference. ",
           "Bars are coloured by channel."),
    day, nrow(al), nrow(fi), sum(res$contrib_days)), B2_FIG_DIR)
  
  invisible(res)
}

## ====
## J: MAIN FIGURES 1-3 
## Figure 1: AL external validation | Figure 2: AL vs IF error at MAIN_DAY |
## Figure 3: trait contributions with % share of total shift
## ====

build_main_figures_B2 <- function(day, res_al, if_proj, attrib,
                                  fig_dir = B2_FIG_DIR, font = FIG_FONT) {
  stopifnot(!is.null(if_proj), !is.null(attrib))
  G_AL <- "AL (external)"; G_IF <- "IF (projected)"
  cols <- setNames(c("#4D4D4D", "#D7191C"), c(G_AL, G_IF))
  
  ## ---- Figure 1: AL external validation ----
  met <- res_al$metrics
  save_fig(res_al$plot, 1L, fig_dir, "large", h = 130)
  register_legend(1L, paste0(
    "External validation of the Batch 1 clock. Predicted versus chronological ",
    "age for Batch 2 ad libitum-fed (AL) worms (n = ", met$n, "), scored with ",
    "the Batch 1-trained clock and the Batch 1 image-noise reference. The ",
    "clock was not refit or tuned on Batch 2. The dashed line marks perfect ",
    "prediction. Spearman ", SYM_RHO, " = ", sprintf("%.2f", met$rho),
    ", Pearson r\u00b2 = ", sprintf("%.2f", met$r2_pearson),
    ", MAE = ", sprintf("%.2f", met$MAE), " d, RMSE = ",
    sprintf("%.2f", met$RMSE), " d."), fig_dir)
  
  ## ---- Figure 2: error at one day; day-specific linear model ----
  cb <- if_proj$comparison
  cb <- cb[cb$day == day, , drop = FALSE]
  cb$group <- factor(cb$group, levels = c(G_AL, G_IF))
  stopifnot(all(table(cb$group) >= 2L))
  
  fit_d <- lm(error ~ group, data = cb)
  cf_d <- summary(fit_d)$coefficients
  rn_d <- grep("^group", rownames(cf_d), value = TRUE)
  est_d <- unname(cf_d[rn_d, 1])
  p_less_d <- pt(unname(cf_d[rn_d, 3]), fit_d$df.residual)
  ci_d <- as.numeric(confint(fit_d)[rn_d, ])
  nA <- sum(cb$group == G_AL); nI <- sum(cb$group == G_IF)
  
  pB <- ggplot(cb, aes(group, error, fill = group, colour = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = LW) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.3,
                 linewidth = LW) +
    geom_jitter(width = 0.12, height = 0, size = PT_SZ, alpha = 0.75) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    scale_x_discrete(labels = c("AL", "IF")) +
    labs(x = NULL, y = "Prediction error (days)") +
    pub_theme + theme(legend.position = "none")
  
  save_fig(pB, 2L, fig_dir, "large", h = 130)
  register_legend(2L, paste0(
    "Prediction error at adult day ", day, " for Batch 2 ad libitum-fed ",
    "(AL; grey; external, n = ", nA, ") and intermittent-fasting (IF; red; ",
    "projected, n = ", nI, ") worms scored with the Batch 1 clock. IF worms ",
    "were assigned a younger age than AL worms (", mn(est_d), " d, 95% CI ",
    mn(ci_d[1]), " to ", mn(ci_d[2]), "; one-sided p = ",
    sprintf("%.3f", p_less_d), ", linear model). Boxes show the median and ",
    "IQR; whiskers extend to 1.5 \u00d7 IQR."), fig_dir)
  
  ## ---- Figure 3: trait contributions with % share ----
  sm <- attrib$summary
  pC <- attr_plot(attrib$table, show_pct = TRUE, font = font)
  
  save_fig(pC, 3L, fig_dir, "large", h = 130)
  register_legend(3L, paste0(
    "Trait contributions to the younger clock age of Batch 2 IF worms at day ",
    day, " (n = ", sm$n_AL, " AL, ", sm$n_IF, " IF), with Batch 1 clock ",
    "weights held fixed. Percentages show each trait's share of the total ",
    "shift (", mn(sm$total_shift_days), " d, 95% CI ", mn(sm$total_ci_low),
    " to ", mn(sm$total_ci_high), "); negative shares oppose the shift. ",
    "Whiskers show 95% bootstrap CIs over worms (",
    format(sm$n_boot, big.mark = ","), " resamples). ",
    if (attrib$n_zero > 0L)
      paste0(attrib$n_zero, " traits with zero clock weight are omitted. ")
    else "",
    "Asterisk: BH-adjusted bootstrap P < 0.05. ",
    "Bars are colored by channel: PGL-1 (red), MUT-16 (green), ",
    "both channels (purple)."), fig_dir)
  
  cat(sprintf("\n[FIGURES 1-3] Day %d: IF - AL = %+.3f d (95%% CI %+.3f to %+.3f); one-sided p = %.4f\n",
              day, est_d, ci_d[1], ci_d[2], p_less_d))
  invisible(list(fig1 = res_al$plot, fig2 = pB, fig3 = pC))
}

## ====
## RUN
## ====

cat("\n", strrep("*", 70),
    "\nBATCH 2 — external validation of a B1-only AL clock\n",
    strrep("*", 70), "\n")

## Freeze the B1-only model before processing B2.
clock_B1 <- train_B1_clock(B1_per_worm, seed = 19)

pairs_AL <- collect_group_pairs(B2_ROOT, "AL", AL_REGEX)
pairs_IF <- collect_group_pairs(B2_ROOT, "IF", IF_REGEX)

if (is.null(pairs_AL)) {
  stop("No B2 AL pairs found; check B2_ROOT / AL_REGEX / filenames.")
}

pw_AL_B2 <- process_group_B2(pairs_AL, "AL")
pw_IF_B2 <- if (!is.null(pairs_IF)) process_group_B2(pairs_IF, "IF") else NULL

al_tbl <- project_B2(clock_B1, pw_AL_B2, "AL")
if_tbl <- if (!is.null(pw_IF_B2)) project_B2(clock_B1, pw_IF_B2, "IF") else NULL

write.csv(al_tbl, file.path(B2_OUT_DIR, "Batch2_AL_projected.csv"),
          row.names = FALSE)
if (!is.null(if_tbl)) {
  write.csv(if_tbl, file.path(B2_OUT_DIR, "Batch2_IF_projected.csv"),
            row.names = FALSE)
}

## AL external validation (Figure 1 plot is built here, saved below)
res_B2_AL <- evaluate_B2_AL(al_tbl)

## Figures 4-5: IF vs AL projection, all shared days (supplementary)
if_projection_B2 <- project_IF_vs_AL_B2(al_tbl, if_tbl, tag = "_B1ref")

## Figures 6+: clock attribution and trait contributions, per shared day
attribution_B2 <- list()

if (!is.null(pw_IF_B2)) {
  shared_days <- sort(intersect(unique(pw_AL_B2$day), unique(pw_IF_B2$day)))
  
  for (d in shared_days) {
    attribution_B2[[paste0("D", d)]] <- tryCatch(
      clock_attribution_B2(pw_AL_B2, pw_IF_B2, clock_B1, day = d,
                           tag = "_B1ref", n_boot = 5000),
      error = function(e) {
        message(sprintf("[ATTRIBUTION] D%d skipped: %s", d, conditionMessage(e)))
        NULL
      }
    )
    compare_traits_B2(pw_AL_B2, pw_IF_B2, clock_B1, day = d)
  }
}

## Figures 1-3: main figures
att_main <- attribution_B2[[paste0("D", MAIN_DAY)]]
if (!is.null(if_projection_B2) && !is.null(att_main)) {
  main_figs_B2 <- build_main_figures_B2(MAIN_DAY, res_B2_AL,
                                        if_projection_B2, att_main)
} else {
  warning("Figures 2-3 not built: no IF data or attribution for day ", MAIN_DAY,
          ". Set MAIN_DAY to a shared day.")
  ## Figure 1 does not depend on IF data
  save_fig(res_B2_AL$plot, 1L, B2_FIG_DIR, "large", h = 130)
}

cat("\nDone. External-validation outputs in:\n", B2_OUT_DIR, "\n")
cat("Figures and Figure_legends.txt in:\n", B2_FIG_DIR, "\n")