# ============================================================
# Age-associated changes in the germ granule architecture (Batch 1 & 2)
# ============================================================
rm(list = ls()); options(warn = 1)

for (p in c("dplyr", "tidyr", "readr", "lme4", "lmerTest", "ggplot2", "scales", "ragg"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("EBImage", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("EBImage", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(EBImage); library(dplyr); library(tidyr); library(readr)
  library(lme4); library(lmerTest); library(ggplot2); library(scales)
})

# ============================================================
# FIGURES
# ============================================================
FONT     <- "Times New Roman"
FS_TXT   <- 10
FS_TICK  <- 8
FS_STRIP <- 7
LW       <- 0.4          # reference lines, error bars
LW_AXIS  <- 0.6         # axis lines and ticks (bold look)
PT_SZ    <- 1.3
ANN      <- FS_TICK / .pt
MM_SMALL <- 80
MM_LARGE <- 180
FIG_DPI  <- 600
FIG_MAX_MB <- 10

SYM_ARROW <- "\u2192"

pub_theme <- theme_classic(base_size = FS_TXT, base_family = FONT) +
  theme(
    text              = element_text(family = FONT, colour = "black"),
    plot.title        = element_blank(),
    plot.subtitle     = element_blank(),
    plot.caption      = element_blank(),
    axis.title        = element_text(face = "bold", size = FS_TXT),
    axis.text         = element_text(face = "bold", size = FS_TICK, colour = "black"),
    axis.line         = element_line(linewidth = LW_AXIS, colour = "black"),
    axis.ticks        = element_line(linewidth = LW_AXIS, colour = "black"),
    legend.title      = element_blank(),
    legend.text       = element_text(size = FS_TICK),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.size   = unit(3, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    strip.text        = element_text(face = "bold", size = FS_STRIP),
    strip.background  = element_blank(),
    panel.grid        = element_blank(),
    plot.margin       = margin(4, 4, 4, 4)
  )
theme_set(pub_theme)

## Bold trait labels (y axis of the forest plot)
trait_bold <- theme(axis.text.y = element_text(face = "bold", size = FS_TICK,
                                               colour = "black"))

update_geom_defaults("text", list(family = FONT, size = ANN, colour = "black"))

save_fig <- function(p, num, fig_dir, size = c("small", "large"), h = NULL) {
  size <- match.arg(size)
  w <- if (size == "small") MM_SMALL else MM_LARGE
  if (is.null(h)) h <- if (size == "small") 70 else 110
  stem <- file.path(fig_dir, paste0("Figure_", num))
  
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
         units = "mm", device = grDevices::cairo_pdf)
  
  tif_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_tiff else "tiff"
  tif <- paste0(stem, ".tiff")
  ggsave(tif, p, width = w, height = h, units = "mm", dpi = FIG_DPI,
         device = tif_dev, compression = "lzw", bg = "white")
  
  mb <- file.size(tif) / 1024^2
  if (mb > FIG_MAX_MB) {
    ggsave(tif, p, width = w, height = h, units = "mm", dpi = 300,
           device = tif_dev, compression = "lzw", bg = "white")
    warning(sprintf("Figure_%s.tiff was %.1f MB at %d dpi; re-saved at 300 dpi.",
                    num, mb, FIG_DPI))
  }
  invisible(NULL)
}

.legend_store <- new.env()
register_legend <- function(num, text, fig_dir) {
  prefix <- paste0(fig_dir, "||")
  assign(paste0(prefix, sprintf("%02d", num)),
         paste0("Figure ", num, ". ", text), envir = .legend_store)
  keys <- sort(ls(.legend_store))
  keys <- keys[startsWith(keys, prefix)]
  writeLines(unlist(mget(keys, envir = .legend_store)),
             file.path(fig_dir, "Figure_legends.txt"), useBytes = FALSE)
}

# ---------------- CONFIG ----------------
ROOT       <- "/Users/purbasha/Desktop/Germ granules/Age clock"
IMG_ROOT   <- file.path(ROOT, "Analysis_raw image")        # B1_D1..B1_D5, B2_D1..B2_D5
SEARCH     <- "/Users/purbasha/Desktop/Germ granules"      # where Batch_1.R / Batch_2.R live
NEW_OUT    <- file.path(ROOT, "reextracted_from_images")
RES_DIR    <- file.path(NEW_OUT, "age_architecture_results")
FIG_DIR    <- file.path(RES_DIR, "figures")
RERUN_EXTRACTION <- TRUE   
EXPECTED_ROWS <- c(Batch1 = 777, Batch2 = 807)   # warning only
for (d in c(NEW_OUT, RES_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

features <- c("red_n","red_size","red_elong","green_n","green_size","green_elong",
              "green_ring","red_ring","red_radial_cv","green_app_net",
              "red_app_net_norm","size_ratio")        # 12 prespecified; intensity excluded
log_features        <- c("green_app_net")             # right-skewed -> log1p
foci_dependent      <- c("green_size","green_elong","green_app_net")  # undefined if green_n == 0
min_animals_per_day <- 3
NA_STR <- c("", "NA", "nan", "NaN")

# ---- Feature names ----
feat_label <- c(
  red_n            = "PGL-1 granule count",
  red_size         = "PGL-1 granule area (px\u00B2)",
  red_elong        = "PGL-1 eccentricity",
  green_n          = "MUT-16 foci count",
  green_size       = "MUT-16 focus area (px\u00B2)",
  green_elong      = "MUT-16 eccentricity",
  green_ring       = "MUT-16 NN-spacing CV",
  red_ring         = "PGL-1 NN-spacing CV",
  red_radial_cv    = "PGL-1 radial dispersion CV",
  green_app_net    = paste("MUT-16", SYM_ARROW, "PGL-1 net association (%)"),
  red_app_net_norm = paste("PGL-1", SYM_ARROW, "MUT-16 net association (norm.)"),
  size_ratio       = "PGL-1 / MUT-16 size ratio")
stopifnot(setequal(names(feat_label), features))
write_csv(tibble(feature_code = names(feat_label), feature_name = unname(feat_label)),
          file.path(RES_DIR, "feature_name_map.csv"))
relabel <- function(x) unname(feat_label[x])

# ============================================================
# PART A. EXTRACTION FROM RAW IMAGES 
# ============================================================
cat("\n================ IMAGE SOURCE CHECK ================\n")
cat("Image root :", IMG_ROOT, "\n")
stopifnot("Image root does not exist" = dir.exists(IMG_ROOT))
stopifnot("Image root must be the RAW folder, not a merged one" =
            grepl("raw", basename(IMG_ROOT), ignore.case = TRUE) &&
            !grepl("merged", IMG_ROOT, ignore.case = TRUE))
cat("Batch/day folders found:\n"); print(sort(basename(list.dirs(IMG_ROOT, recursive = FALSE))))

find_script <- function(fname, label) {
  hits <- list.files(SEARCH, pattern = paste0("^", fname, "$"), recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (length(hits)) { cat(label, "->", hits[1], "\n"); return(hits[1]) }
  cat(label, "not found automatically. Select the file...\n"); file.choose()
}

load_batch_env <- function(script, stop_pattern) {
  E <- new.env(parent = globalenv())
  exprs <- parse(script, keep.source = FALSE)
  txt   <- vapply(exprs, function(e) deparse(e)[1], character(1))
  stop_i <- grep(stop_pattern, txt)[1]
  if (is.na(stop_i)) stop("No line matching '", stop_pattern, "' in ", script)
  for (i in seq_len(stop_i - 1)) {
    if (grepl("^rm\\(list", txt[i])) next
    eval(exprs[[i]], envir = E)
  }
  need <- c("extract_day","extract_worm_num","find_nuclei_flat",
            "compute_image_noise_scale","extract_features_real","LOCAL_CONTRAST_Z_MIN")
  miss <- need[!vapply(need, exists, logical(1), envir = E)]
  if (length(miss)) stop("Missing in ", basename(script), ": ", paste(miss, collapse = ", "))
  E$OUT_ROOT <- NEW_OUT
  E
}

process_all_images <- function(root_dir, bt, label) {
  day_dirs <- list.dirs(root_dir, recursive = FALSE)
  day_dirs <- day_dirs[grepl(paste0("^", bt, "_D[0-9]+$"), basename(day_dirs))]
  if (!length(day_dirs)) stop("No ", bt, "_D* folders under: ", root_dir)
  
  prs <- list()
  for (dd in day_dirs) {
    wfs <- list.dirs(dd, recursive = FALSE)
    wfs <- wfs[grepl(paste0("^AL_", bt, "_D[0-9]+_W[0-9]+$"), basename(wfs))]
    for (wf in wfs) {
      day <- extract_day(wf); worm_no <- extract_worm_num(wf)
      if (is.na(day) || is.na(worm_no)) next
      red <- find_nuclei_flat(wf, "R"); green <- find_nuclei_flat(wf, "G")
      ids <- intersect(red$nucleus, green$nucleus)
      if (!length(ids)) next
      prs[[length(prs) + 1]] <- data.frame(
        worm_id = paste0("adlib_D", day, "_W", worm_no), condition = "adlib",
        day = as.integer(day), worm = worm_no, nucleus = ids,
        base_id = paste0("AL_", bt, "_D", day, "_W", worm_no, "_n", ids),
        red_path = red$path[match(ids, red$nucleus)],
        green_path = green$path[match(ids, green$nucleus)],
        stringsAsFactors = FALSE)
    }
  }
  pair_tbl <- do.call(rbind, prs)
  stopifnot(!is.null(pair_tbl), nrow(pair_tbl) > 0)
  
  all_paths <- normalizePath(c(pair_tbl$red_path, pair_tbl$green_path), mustWork = FALSE)
  root_norm <- normalizePath(root_dir, mustWork = FALSE)
  in_root   <- startsWith(all_paths, root_norm)
  has_merge <- grepl("merged", all_paths, ignore.case = TRUE)
  cat(sprintf("\n[CONFIRM] %s: %d red/green nucleus pairs from %d animals\n",
              label, nrow(pair_tbl), length(unique(pair_tbl$worm_id))))
  cat(sprintf("[CONFIRM] %d/%d image paths are inside the RAW root; %d contain 'merged'\n",
              sum(in_root), length(in_root), sum(has_merge)))
  stopifnot("Some images are outside the raw folder" = all(in_root),
            "Some image paths contain 'merged'"       = !any(has_merge))
  cat("[CONFIRM] Images per day / animals per day:\n")
  print(pair_tbl %>% group_by(day) %>%
          summarise(animals = n_distinct(worm_id), nuclei = n(), .groups = "drop"))
  cat("[CONFIRM] Example paths:\n"); print(head(pair_tbl$red_path, 2)); print(head(pair_tbl$green_path, 2))
  
  ni <- compute_image_noise_scale(pair_tbl)
  pair_tbl <- merge(pair_tbl, ni$scale_tbl[, c("base_id", "noise_ratio")],
                    by = "base_id", all.x = TRUE)
  pair_tbl$noise_ratio[!is.finite(pair_tbl$noise_ratio)] <- 1
  pair_tbl$adaptive_z_min <- LOCAL_CONTRAST_Z_MIN * pair_tbl$noise_ratio
  
  ov <- file.path(NEW_OUT, label, "qc_overlays")
  dir.create(ov, showWarnings = FALSE, recursive = TRUE)
  feats <- vector("list", nrow(pair_tbl))
  for (i in seq_len(nrow(pair_tbl))) {
    if (i %% 20 == 0 || i == 1) cat(sprintf("  [%d/%d] %s\n", i, nrow(pair_tbl), pair_tbl$base_id[i]))
    meta <- as.list(pair_tbl[i, c("condition","day","worm_id","nucleus","base_id")])
    feats[[i]] <- tryCatch(
      extract_features_real(pair_tbl$red_path[i], pair_tbl$green_path[i], meta,
                            overlay_dir = ov, contrast_z_min = pair_tbl$adaptive_z_min[i]),
      error = function(e) { warning("SKIP ", pair_tbl$base_id[i], ": ", e$message); NULL })
  }
  per_image <- do.call(rbind, Filter(Negate(is.null), feats))
  write.csv(per_image, file.path(NEW_OUT, paste0(label, "_per_image_FROM_IMAGES.csv")), row.names = FALSE)
  cat(sprintf("[CONFIRM] %s: %d of %d pairs yielded features\n", label, nrow(per_image), nrow(pair_tbl)))
  per_image
}

run_in_env <- function(E, ...) { f <- process_all_images; environment(f) <- E; f(...) }

if (RERUN_EXTRACTION) {
  B1_SCRIPT <- find_script("Batch_1.R", "Batch 1 script")
  B2_SCRIPT <- find_script("Batch_2.R", "Batch 2 script")
  E1 <- load_batch_env(B1_SCRIPT, "^b1 <- process_batch1")
  run_in_env(E1, IMG_ROOT, "B1", "Batch1")
  E2 <- load_batch_env(B2_SCRIPT, "^b2 <- process_batch2")
  run_in_env(E2, IMG_ROOT, "B2", "Batch2")
}

# ============================================================
# PART B. NESTED AGE MODELS (nuclei within animals)
# ============================================================
files_raw <- c(Batch1 = file.path(NEW_OUT, "Batch1_per_image_FROM_IMAGES.csv"),
               Batch2 = file.path(NEW_OUT, "Batch2_per_image_FROM_IMAGES.csv"))
stopifnot(all(file.exists(files_raw)))

cat("\n================ MODEL INPUT CHECK ================\n")
for (b in names(files_raw)) {
  n <- nrow(read_csv(files_raw[[b]], na = NA_STR, show_col_types = FALSE))
  cat(sprintf("[CONFIRM] %s model input: %s\n          rows = %d | written %s\n", b,
              files_raw[[b]], n, format(file.info(files_raw[[b]])$mtime)))
  stopifnot("Model input is not a FROM_IMAGES (raw-folder) table" = grepl("FROM_IMAGES", files_raw[[b]]))
  if (n != EXPECTED_ROWS[[b]])
    warning(sprintf("%s has %d rows; expected %d from the earlier raw run", b, n, EXPECTED_ROWS[[b]]))
}

load_batch <- function(path, b, mask_no_foci = TRUE) {
  d <- read_csv(path, na = NA_STR, show_col_types = FALSE) %>%
    mutate(batch = b,
           animal_uid = paste(b, worm_id, sep = "_"),
           day_num = as.numeric(day))
  if (mask_no_foci)
    d <- d %>% mutate(across(all_of(foci_dependent),
                             ~ ifelse(!is.na(green_n) & green_n == 0, NA_real_, .x)))
  d %>% mutate(across(all_of(intersect(log_features, features)), ~ log1p(pmax(.x, 0))))
}

zscore <- function(d) d %>% mutate(across(all_of(features), ~ as.numeric(scale(.x))))

analyse <- function(d, adjust_batch = FALSE) {
  d <- zscore(d)
  out <- bind_rows(lapply(features, function(f) {
    x <- d %>% select(y = all_of(f), day_num, batch, animal_uid) %>% filter(is.finite(y))
    x$day_f <- factor(x$day_num)
    base <- tibble(feature = f, n_nuclei = nrow(x), n_animals = n_distinct(x$animal_uid),
                   slope_SD_per_day = NA_real_, SE = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
                   p_slope = NA_real_, p_day_categorical = NA_real_, ICC_animal = NA_real_,
                   rho_animal = NA_real_, p_animal = NA_real_, singular = NA, note = NA_character_)
    per_day <- x %>% group_by(day_num) %>% summarise(a = n_distinct(animal_uid), .groups = "drop")
    if (nrow(per_day) < 3 || any(per_day$a < min_animals_per_day)) {
      base$note <- "skipped: too few animals on one or more days"; return(base)
    }
    f_lin  <- if (adjust_batch) y ~ day_num + batch + (1 | animal_uid) else y ~ day_num + (1 | animal_uid)
    f_full <- if (adjust_batch) y ~ day_f   + batch + (1 | animal_uid) else y ~ day_f   + (1 | animal_uid)
    m_lin  <- suppressMessages(lmer(f_lin,  data = x))
    m_full <- suppressMessages(lmer(f_full, data = x))
    co <- summary(m_lin)$coefficients
    ci <- suppressMessages(confint(m_lin, parm = "day_num", method = "Wald"))
    vc <- as.data.frame(VarCorr(m_lin))
    va <- vc$vcov[vc$grp == "animal_uid"]; vr <- vc$vcov[vc$grp == "Residual"]
    am <- x %>% group_by(animal_uid, day_num) %>% summarise(y = mean(y), .groups = "drop")
    ct <- suppressWarnings(cor.test(am$y, am$day_num, method = "spearman", exact = FALSE))
    base$slope_SD_per_day <- co["day_num", "Estimate"]
    base$SE        <- co["day_num", "Std. Error"]
    base$CI_low    <- ci["day_num", 1]; base$CI_high <- ci["day_num", 2]
    base$p_slope   <- co["day_num", "Pr(>|t|)"]
    base$p_day_categorical <- anova(m_full, type = "III")["day_f", "Pr(>F)"]
    base$ICC_animal <- va / (va + vr)
    base$rho_animal <- unname(ct$estimate); base$p_animal <- ct$p.value
    base$singular   <- isSingular(m_lin)
    base
  }))
  out %>%
    mutate(padj_slope           = p.adjust(p_slope, "BH"),
           padj_day_categorical = p.adjust(p_day_categorical, "BH"),
           padj_animal          = p.adjust(p_animal, "BH"),
           direction = ifelse(is.na(slope_SD_per_day), NA_character_,
                              ifelse(slope_SD_per_day < 0, "decrease", "increase")),
           feature = relabel(feature))
}

run_all <- function(files, tag, mask_no_foci = TRUE) {
  dats <- lapply(names(files), function(b) load_batch(files[[b]], b, mask_no_foci))
  names(dats) <- names(files)
  for (b in names(dats)) {
    cat("\n=====", tag, b, "=====\n")
    print(dats[[b]] %>% group_by(day_num) %>%
            summarise(animals = n_distinct(animal_uid), nuclei = n(), .groups = "drop"))
  }
  res_b <- lapply(names(dats), function(b) analyse(dats[[b]]) %>% mutate(batch = b))
  names(res_b) <- names(dats)
  res_all <- analyse(bind_rows(dats), adjust_batch = TRUE) %>% arrange(p_slope)
  for (b in names(res_b)) write_csv(res_b[[b]], file.path(RES_DIR, paste0(tag, "_", b, "_age_mixed_model.csv")))
  write_csv(res_all, file.path(RES_DIR, paste0(tag, "_Combined_age_mixed_model.csv")))
  list(dats = dats, res_b = res_b, res_all = res_all)
}

R_raw  <- run_all(files_raw, "RAW",               mask_no_foci = TRUE)
R_zero <- run_all(files_raw, "RAW_zerosIncluded", mask_no_foci = FALSE)

cat("\n===== RAW-FOLDER COMBINED MODEL (primary: foci-conditional) =====\n")
print(R_raw$res_all, width = Inf)
for (b in names(R_raw$res_b)) { cat("\n==", b, "==\n"); print(R_raw$res_b[[b]], width = Inf) }

no_foci_tbl <- bind_rows(R_raw$dats) %>%
  group_by(batch, day_num) %>%
  summarise(n_images = n(), n_no_foci = sum(green_n == 0, na.rm = TRUE),
            pct_no_foci = round(100 * n_no_foci / n_images, 1), .groups = "drop")
write_csv(no_foci_tbl, file.path(RES_DIR, "no_foci_images_by_day.csv"))
cat("\n===== Images with no MUT-16 foci =====\n"); print(no_foci_tbl)

sens_tbl <- R_raw$res_all %>%
  select(feature, n_cond = n_nuclei, slope_cond = slope_SD_per_day, padj_cond = padj_slope) %>%
  inner_join(R_zero$res_all %>%
               select(feature, n_zero = n_nuclei, slope_zero = slope_SD_per_day, padj_zero = padj_slope),
             by = "feature") %>%
  mutate(changed_rows = n_zero != n_cond,
         same_direction = sign(slope_cond) == sign(slope_zero),
         same_significance = (padj_cond < 0.05) == (padj_zero < 0.05))
write_csv(sens_tbl, file.path(RES_DIR, "sensitivity_zero_vs_conditional.csv"))
cat("\n===== SENSITIVITY: foci-conditional vs zeros included (combined) =====\n")
print(sens_tbl %>% filter(changed_rows), width = Inf)

# ============================================================
# PART C1. FIGURE 1 - forest plot (primary analysis)
# * = combined model (batch as fixed effect) BH-adj. P < 0.05
# ============================================================
plot_df <- bind_rows(
  R_raw$res_b$Batch1 %>% mutate(set = "Batch 1"),
  R_raw$res_b$Batch2 %>% mutate(set = "Batch 2"),
  R_raw$res_all      %>% mutate(set = "Combined (batch-adjusted)")
) %>% filter(!is.na(slope_SD_per_day))

sig <- R_raw$res_all %>% filter(!is.na(padj_slope), padj_slope < 0.05) %>% pull(feature)
plot_df <- plot_df %>% mutate(age_assoc = feature %in% sig)

ord <- R_raw$res_all %>% filter(!is.na(slope_SD_per_day)) %>%
  arrange(slope_SD_per_day) %>% pull(feature)
lab <- setNames(ifelse(ord %in% sig, paste0(ord, " *"), ord), ord)
wrap_lab <- function(x) vapply(lab[x], function(s) paste(strwrap(s, 34), collapse = "\n"),
                               character(1), USE.NAMES = FALSE)

plot_df$feature <- factor(plot_df$feature, levels = ord)
plot_df$set <- factor(plot_df$set, levels = c("Batch 1", "Batch 2", "Combined (batch-adjusted)"))

p_all <- ggplot(plot_df, aes(slope_SD_per_day, feature, colour = set,
                             shape = padj_slope < 0.05, alpha = age_assoc)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = LW) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high), size = 0.3, linewidth = LW,
                  position = position_dodge(width = 0.7)) +
  scale_y_discrete(labels = wrap_lab) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH adj. P < 0.05",
                                `FALSE` = "BH adj. P \u2265 0.05")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "black")) +
  labs(x = "Age slope (SD per day, 95% CI)", y = NULL) +
  pub_theme + trait_bold +
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.y = element_text(face = "bold", size = FS_TICK,
                                   colour = "black", lineheight = 0.9))

save_fig(p_all, 1, FIG_DIR, "large", h = 130)
register_legend(1, sprintf(
  paste0("Age-associated changes in germ-granule architecture. Points: age slope ",
         "(SD units per day) from nested linear mixed models (nuclei within animals; ",
         "animal random intercept) with 95%% Wald CIs, shown for Batch 1, Batch 2 and ",
         "the combined data (batch as fixed effect). Unit of replication: animal ",
         "(Batch 1 + Batch 2: %d animals). Features are z-scored; MUT-16 focus area, ",
         "eccentricity and net association are undefined and set to missing in images ",
         "with no MUT-16 foci. Filled symbols: Benjamini-Hochberg adjusted P < 0.05 ",
         "within each model family; * = combined-model adjusted P < 0.05 (faded ",
         "features: no combined age effect). Dashed line: zero slope."),
  max(R_raw$res_all$n_animals, na.rm = TRUE)), FIG_DIR)

# ---- Do age slopes differ between batches? (day x batch; primary data) ----
dall <- bind_rows(R_raw$dats) %>% mutate(across(all_of(features), ~ as.numeric(scale(.x))))
int_tbl <- bind_rows(lapply(features, function(f) {
  x <- dall %>% transmute(y = .data[[f]], day_num, batch, animal_uid) %>% filter(is.finite(y))
  if (n_distinct(x$animal_uid[x$batch == "Batch1"]) < 10 ||
      n_distinct(x$animal_uid[x$batch == "Batch2"]) < 10)
    return(tibble(feature = f, p_day_x_batch = NA_real_))
  m0 <- suppressMessages(lmer(y ~ day_num + batch + (1 | animal_uid), data = x, REML = FALSE))
  m1 <- suppressMessages(lmer(y ~ day_num * batch + (1 | animal_uid), data = x, REML = FALSE))
  tibble(feature = f, p_day_x_batch = anova(m0, m1)$`Pr(>Chisq)`[2])
})) %>% mutate(padj_day_x_batch = p.adjust(p_day_x_batch, "BH"),
               feature = relabel(feature))
write_csv(int_tbl, file.path(RES_DIR, "batch_by_age_interaction.csv"))
print(int_tbl)

# ============================================================
# PART C2. FIGURE 2 - per-day plots, significant traits only
# ============================================================
stopifnot(length(sig) > 0)
back <- function(code, v) ifelse(code %in% log_features, expm1(v), v)

sig_tbl <- R_raw$res_all %>% filter(feature %in% sig) %>% arrange(padj_slope)
feat_lv <- sig_tbl$feature

traits_sig <- bind_rows(R_raw$dats) %>%
  select(batch, animal_uid, day_num, all_of(features)) %>%
  pivot_longer(all_of(features), names_to = "code", values_to = "value") %>%
  mutate(feature = relabel(code)) %>%
  filter(feature %in% sig, is.finite(value)) %>%
  mutate(value = back(code, value)) %>%
  group_by(batch, animal_uid, day_num, feature) %>%            # animal = unit
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(feature = factor(feature, levels = feat_lv))

day_sum <- traits_sig %>% group_by(feature, day_num) %>%
  summarise(m = mean(value), se = sd(value) / sqrt(n()), .groups = "drop")

rng_tbl <- traits_sig %>% group_by(feature) %>%
  summarise(nonneg = all(value >= 0), lo = min(value), hi = max(value), .groups = "drop")

fit_df <- traits_sig %>% group_by(feature) %>%
  group_modify(~ {
    m  <- lm(value ~ day_num, data = .x)
    nd <- data.frame(day_num = seq(1, 5, length.out = 50))
    cbind(nd, as.data.frame(predict(m, nd, interval = "confidence")))
  }) %>% ungroup() %>%
  left_join(rng_tbl, by = "feature") %>%
  mutate(lwr = ifelse(nonneg, pmax(lwr, 0), lwr))

lims <- fit_df %>% group_by(feature) %>%
  summarise(flo = min(lwr), fhi = max(upr), .groups = "drop") %>%
  left_join(rng_tbl, by = "feature") %>%
  mutate(lo_ = pmin(lo, flo), hi_ = pmax(hi, fhi), rng = hi_ - lo_,
         ylo = lo_ - 0.05 * rng, ylo = ifelse(nonneg, pmax(ylo, 0), ylo),
         yhi = hi_ + 0.05 * rng)

n_f    <- length(feat_lv)
ncol_f <- if (n_f == 4) 2 else min(3, n_f)
nrow_f <- ceiling(n_f / ncol_f)
spare  <- (n_f %% ncol_f) != 0

p_sig <- ggplot(traits_sig, aes(day_num, value)) +
  geom_blank(data = lims, aes(x = 3, y = ylo), inherit.aes = FALSE) +
  geom_blank(data = lims, aes(x = 3, y = yhi), inherit.aes = FALSE) +
  geom_jitter(aes(colour = batch), width = 0.12, alpha = 0.55, size = PT_SZ, stroke = 0) +
  geom_ribbon(data = fit_df, aes(y = fit, ymin = lwr, ymax = upr),
              fill = "grey70", alpha = 0.3) +
  geom_line(data = fit_df, aes(y = fit), colour = "black", linewidth = 0.5) +
  geom_errorbar(data = day_sum, aes(y = m, ymin = m - se, ymax = m + se),
                width = 0.15, linewidth = LW) +
  geom_point(data = day_sum, aes(y = m), shape = 21, fill = "white",
             colour = "black", size = PT_SZ + 0.4, stroke = 0.5) +
  facet_wrap(~ feature, scales = "free_y", ncol = ncol_f, axes = "all_x",
             labeller = labeller(feature = label_wrap_gen(width = 30))) +
  scale_colour_manual(values = c(Batch1 = "#1b9e77", Batch2 = "#d95f02"),
                      labels = c(Batch1 = "Batch 1", Batch2 = "Batch 2")) +
  scale_x_continuous(breaks = 1:5, labels = paste0("D", 1:5),
                     limits = c(0.6, 5.4), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Adult day", y = "Animal mean (original units)") +
  pub_theme +
  theme(legend.position = if (spare) "inside" else "bottom",
        legend.position.inside = c(0.84, 0.14),
        panel.spacing = unit(1.1, "lines"))

h_sig <- if (nrow_f <= 2) 130 else min(220, 62 * nrow_f + 15)
save_fig(p_sig, 2, FIG_DIR, "large", h = h_sig)
register_legend(2, paste0(
  "Per-day values of the germ-granule features with a significant combined age ",
  "association (BH-adjusted P < 0.05; see Figure 1). Small points: individual animals ",
  "(jittered; colored by batch); open circles: day mean +/- SE of animal means (unit of ",
  "replication: animal; nuclei averaged within animal). Line and band: linear fit ",
  "with 95% CI, a visual summary of the linear trend (some features are not strictly ",
  "linear across days)."),
  FIG_DIR)

print(p_all); print(p_sig)
cat("\nDone. Figures and legends in:\n", FIG_DIR, "\n")