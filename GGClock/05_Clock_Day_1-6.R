## ====
## GG Clock — Batch 1 (DAYS 1-6)
## ====

rm(list = ls())
options(warn = 1)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("EBImage", quietly = TRUE)) {
  BiocManager::install("EBImage", update = FALSE, ask = FALSE)
}
if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
  install.packages("ggbeeswarm")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

suppressMessages({
  library(EBImage)
  library(glmnet)
  library(ggplot2)
  library(ggbeeswarm)
  library(patchwork)
})

set.seed(19)

## ====
## FIGURES
## ====
FONT   <- "Helvetica"   # built-in PDF font
FS_TXT <- 11            # axis titles, legend titles, annotations (pt)
FS_TICK<- 10            # tick labels, legend text (pt)
FS_TTL <- 12            # plot title (pt, bold)
FS_SUB <- 10            # subtitle (pt)
FS_STRIP <- 10          # facet strip labels (pt)
LW     <- 0.5           # axes, ticks, bar outlines, reference lines
PT_SZ  <- 2.4           # all scatter points
ANN    <- FS_TXT / .pt  # geom_text / annotate size

## 
SYM_RHO   <- "rho"
SYM_ARROW <- "->"
SYM_R2    <- "r^2"

pub_theme <- theme_classic(base_size = FS_TXT, base_family = FONT) +
  theme(
    text              = element_text(family = FONT, colour = "black"),
    plot.title        = element_text(face = "bold", size = FS_TTL, hjust = 0.5),
    plot.subtitle     = element_text(size = FS_SUB, hjust = 0.5, colour = "grey30"),
    axis.title        = element_text(face = "bold", size = FS_TXT),
    axis.text         = element_text(size = FS_TICK, colour = "black"),
    axis.line         = element_line(linewidth = LW),
    axis.ticks        = element_line(linewidth = LW),
    legend.title      = element_text(face = "bold", size = FS_TXT),
    legend.text       = element_text(size = FS_TICK),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    strip.text        = element_text(face = "bold", size = FS_STRIP),
    strip.background  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(12, 12, 12, 12)
  )
theme_set(pub_theme)

## Makes geom_text() and annotate("text") use the same font and size
update_geom_defaults("text", list(family = FONT, size = ANN, colour = "black"))

## Day colours (day 6 added)
col_day <- c(
  "1" = "#4393C3",
  "2" = "#74C476",
  "3" = "#FD8D3C",
  "4" = "#9E66AB",
  "5" = "#D6604D",
  "6" = "#8C6D31"
)

save_fig <- function(p, name, fig_dir, w = 6.5, h = 5) {   # one default size
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p,
         width = w, height = h, device = "pdf")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  invisible(NULL)
}

## ====
## DEFINED PATHS AND CONFIGURATION
## ====

B1_ROOT <- paste0(
  "/Users/purbasha/Desktop/Germ granules/",
  "Age clock/Analysis_raw image"
)

## Days present in the AL data.
EXPECTED_DAYS <- 1:6

##Change RUN_TAG for every new run.
RUN_TAG  <- "run3"
OUT_ROOT <- file.path(
  B1_ROOT,
  paste0("results_v6_B1_D1to6_", RUN_TAG)
)

## Safety guard: never write into a folder that already has files.
if (dir.exists(OUT_ROOT) &&
    length(list.files(OUT_ROOT, recursive = TRUE, all.files = TRUE)) > 0L) {
  stop(
    "Output folder already exists and is not empty:\n  ", OUT_ROOT,
    "\nChange RUN_TAG so that no existing results are overwritten.",
    call. = FALSE
  )
}

BLEED_COEFF <- 0

GREEN_THR_MULT        <- 2.0
GREEN_GLOBAL_MULT     <- 1.3
GREEN_MIN_AREA        <- 20
GREEN_MAX_AREA        <- 1000
GREEN_MAX_ECCENTRICITY <- 0.90
LOCAL_BG_RING_PX      <- 4
LOCAL_CONTRAST_Z_MIN  <- 1.0
GREEN_N_THRESHOLD    <- 4

## effective_z_min = LOCAL_CONTRAST_Z_MIN * noise_ratio
## noise_ratio = image median local_bg_mad / batch reference
CLOCK_SEED <- 19L
CLOCK_ALPHA <- 0.5
CLOCK_INNER_FOLDS <- 10L
CLOCK_STABILITY_SEEDS <- 1:30

CLOCK_LAMBDA <- 10^seq(2, -5, length.out = 100)

RUN_OLD_REF_COMPARISON <- FALSE

## ====
## 
## ====

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}

finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2L) sd(x) else NA_real_
}

finite_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x)) {
    unname(quantile(x, probs = prob))
  } else {
    NA_real_
  }
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2L ||
      length(unique(x)) < 2L ||
      length(unique(y)) < 2L) {
    return(NA_real_)
  }
  
  ans <- suppressWarnings(cor(x, y, method = method))
  if (is.finite(ans)) ans else NA_real_
}

clean_numeric_nonfinite <- function(d) {
  for (nm in names(d)) {
    if (is.numeric(d[[nm]])) {
      d[[nm]][!is.finite(d[[nm]])] <- NA_real_
    }
  }
  d
}

validate_noise_reference <- function(x, label) {
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(
      label,
      " is not a finite positive number. ",
      "Inspect the image-noise measurements before proceeding.",
      call. = FALSE
    )
  }
  invisible(x)
}

## ====
## SECTION A: IMAGE AND SEGMENTATION HELPERS
## ====

as_mat <- function(x) {
  a <- imageData(x)
  matrix(as.double(a), nrow = dim(a)[1], ncol = dim(a)[2])
}

as_lab <- function(x) {
  m <- as_mat(x)
  storage.mode(m) <- "integer"
  m
}

eb <- function(m) Image(m, colormode = "Grayscale")

read_single_tiff <- function(path, channel_name) {
  if (!grepl("\\.(tif|tiff)$", path, ignore.case = TRUE)) {
    stop(channel_name, ": select a TIFF file.", call. = FALSE)
  }
  
  im <- readImage(path)
  
  if (colorMode(im) != Grayscale) {
    stop(
      channel_name,
      ": this is a color/RGB image. Export the original grayscale channel.",
      call. = FALSE
    )
  }
  
  a <- imageData(im)
  d <- dim(a)
  
  if (length(d) < 2L ||
      (length(d) > 2L && any(d[-c(1L, 2L)] != 1L))) {
    stop(
      channel_name,
      ": multiple planes/channels detected. ",
      "Export one grayscale optical section.",
      call. = FALSE
    )
  }
  
  m <- matrix(as.double(a), nrow = d[1], ncol = d[2])
  
  if (any(!is.finite(m))) {
    stop(channel_name, ": non-finite pixels detected.", call. = FALSE)
  }
  if (any(m < 0)) {
    stop(channel_name, ": negative input intensities detected.", call. = FALSE)
  }
  if (min(dim(m)) < 25L) {
    stop(
      channel_name,
      ": image is smaller than 25 pixels on one axis.",
      call. = FALSE
    )
  }
  if (diff(range(m)) == 0) {
    warning(channel_name, ": image has constant intensity.")
  }
  
  cat(sprintf(
    "%s: %s\n  Array dimensions: %d x %d | Loaded range: %.6g to %.6g\n",
    channel_name, basename(path), nrow(m), ncol(m), min(m), max(m)
  ))
  
  m
}

safe_gblur_m <- function(m, sigma) {
  min_dim <- min(nrow(m), ncol(m))
  max_sigma <- max(0.3, floor((min_dim - 3) / 2) / 3)
  sigma <- min(sigma, max_sigma)
  
  if (sigma < 0.3 || min_dim < 7) return(m)
  
  as_mat(gblur(eb(m), sigma = sigma))
}

safe_brush <- function(size, min_dim, shape = "disc") {
  size <- min(as.integer(size), min_dim - 2L)
  if (size < 1L) size <- 1L
  if (size %% 2L == 0L) size <- size - 1L
  if (size < 1L) size <- 1L
  makeBrush(size, shape)
}

safe_opening_m <- function(m, brush) {
  as_mat(opening(eb(m), brush))
}

safe_closing_m <- function(m, brush) {
  as_mat(closing(eb(m), brush))
}

safe_fillhull_m <- function(m) {
  as_mat(fillHull(eb(m)))
}

safe_bwlabel_m <- function(m) {
  as_lab(bwlabel(eb(m)))
}

safe_rmobj <- function(lab_m, bad_ids) {
  if (!length(bad_ids)) return(lab_m)
  as_lab(rmObjects(eb(lab_m), bad_ids, reenumerate = TRUE))
}

safe_erode_m <- function(m, brush) {
  as_mat(erode(eb(m), brush))
}

safe_dilate_m <- function(m, brush) {
  as_mat(dilate(eb(m), brush))
}

segment_red_unified <- function(channel, brush_size = 51, thr_mult = 1.0,
                                min_area = 10, max_area = 3000) {
  min_dim <- min(nrow(channel), ncol(channel))
  sm <- safe_gblur_m(channel, sigma = 1)
  bg <- safe_opening_m(sm, safe_brush(brush_size, min_dim))
  th <- pmax(sm - bg, 0)
  pos <- th[th > 0]
  
  if (length(pos) < 10) {
    return(matrix(0L, nrow(channel), ncol(channel)))
  }
  
  thr <- mean(pos) + thr_mult * sd(pos)
  bw <- safe_fillhull_m((th > thr) * 1.0)
  bw <- safe_closing_m(bw, safe_brush(3, min_dim))
  lab <- safe_bwlabel_m(bw)
  
  if (max(lab) == 0) return(lab)
  
  sh <- computeFeatures.shape(eb(lab))
  bad_ids <- which(
    sh[, "s.area"] < min_area |
      sh[, "s.area"] > max_area
  )
  
  safe_rmobj(lab, bad_ids)
}

local_contrast_stats <- function(lab, channel, ring_px = LOCAL_BG_RING_PX) {
  nobj <- max(lab)
  if (nobj == 0) return(NULL)
  
  min_dim <- min(nrow(lab), ncol(lab))
  brush <- safe_brush(2 * ring_px + 1, min_dim)
  
  out <- data.frame(
    candidate = seq_len(nobj),
    local_bg_median = NA_real_,
    local_bg_mad = NA_real_,
    contrast_z = NA_real_
  )
  
  for (i in seq_len(nobj)) {
    obj_mask <- (lab == i) * 1.0
    ring_mask <- (safe_dilate_m(obj_mask, brush) > 0) &
      (obj_mask == 0)
    
    bg_vals <- channel[ring_mask]
    bg_vals <- bg_vals[is.finite(bg_vals)]
    cand_mean <- mean(channel[lab == i])
    
    if (length(bg_vals) >= 5) {
      bg_med <- median(bg_vals)
      bg_mad <- mad(bg_vals, constant = 1.4826)
      
      if (!is.finite(bg_mad) || bg_mad < 1e-9) {
        bg_mad <- 1e-9
      }
      
      out$local_bg_median[i] <- bg_med
      out$local_bg_mad[i] <- bg_mad
      out$contrast_z[i] <- (cand_mean - bg_med) / bg_mad
    }
  }
  
  out
}

shape_quality_stats <- function(lab, channel) {
  nobj <- max(lab)
  if (nobj == 0) return(NULL)
  
  out <- data.frame(
    candidate = seq_len(nobj),
    peakiness = NA_real_,
    solidity = NA_real_
  )
  
  for (i in seq_len(nobj)) {
    mask <- lab == i
    vals <- channel[mask]
    vals <- vals[is.finite(vals)]
    
    if (length(vals) >= 1 && mean(vals) > 1e-12) {
      out$peakiness[i] <- max(vals) / mean(vals)
    }
    
    hull_area <- tryCatch({
      idx <- which(mask, arr.ind = TRUE)
      
      if (nrow(idx) < 3) {
        NA_real_
      } else {
        ch <- chull(idx[, 1], idx[, 2])
        pts <- idx[ch, , drop = FALSE]
        x <- pts[, 1]
        y <- pts[, 2]
        
        0.5 * abs(sum(
          x * c(y[-1], y[1]) -
            c(x[-1], x[1]) * y
        ))
      }
    }, error = function(e) NA_real_)
    
    area_i <- sum(mask)
    
    if (is.finite(hull_area) && hull_area > 1e-9) {
      out$solidity[i] <- area_i / hull_area
    }
  }
  
  ## 
  out
}

empty_candidate_table <- function() {
  data.frame(
    candidate = integer(),
    area_px = numeric(),
    mean_green = numeric(),
    eccentricity = numeric(),
    cx = numeric(),
    cy = numeric(),
    local_bg_median = numeric(),
    local_bg_mad = numeric(),
    contrast_z = numeric(),
    global_floor = numeric(),
    floor_gap = numeric(),
    peakiness = numeric(),
    solidity = numeric(),
    adaptive_z_min = numeric(),
    pass_area = logical(),
    pass_global_gate = logical(),
    pass_eccentricity = logical(),
    pass_local_contrast = logical(),
    kept = logical(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

diagnose_green <- function(result, contrast_z_min = LOCAL_CONTRAST_Z_MIN) {
  if (length(contrast_z_min) != 1L ||
      !is.finite(contrast_z_min) ||
      contrast_z_min < 0) {
    stop("contrast_z_min must be one finite non-negative number.")
  }
  
  channel <- result$G_seg
  nr <- nrow(channel)
  nc <- ncol(channel)
  min_dim <- min(nr, nc)
  
  fine <- safe_gblur_m(channel, sigma = 1)
  coarse <- safe_gblur_m(channel, sigma = min(8, min_dim / 6))
  dog <- pmax(fine - coarse, 0)
  pos <- dog[dog > 0]
  
  lab <- matrix(0L, nr, nc)
  thr <- NA_real_
  
  if (length(pos) >= 10L) {
    thr <- mean(pos) + GREEN_THR_MULT * sd(pos)
    bw <- safe_fillhull_m((dog > thr) * 1.0)
    bw <- safe_closing_m(bw, safe_brush(3, min_dim))
    lab <- safe_bwlabel_m(bw)
  }
  
  global_floor <- mean(channel) + GREEN_GLOBAL_MULT * sd(channel)
  checks <- empty_candidate_table()
  final_lab <- lab
  
  if (max(lab) > 0L) {
    sh <- computeFeatures.shape(eb(lab))
    mom <- computeFeatures.moment(eb(lab))
    ints <- computeFeatures.basic(eb(lab), ref = eb(channel))
    
    lc <- local_contrast_stats(
      lab, channel, ring_px = LOCAL_BG_RING_PX
    )
    sq <- shape_quality_stats(lab, channel)
    
    checks <- data.frame(
      candidate = seq_len(max(lab)),
      area_px = sh[, "s.area"],
      mean_green = ints[, "b.mean"],
      eccentricity = mom[, "m.eccentricity"],
      cx = mom[, "m.cx"],
      cy = mom[, "m.cy"],
      local_bg_median = lc$local_bg_median,
      local_bg_mad = lc$local_bg_mad,
      contrast_z = lc$contrast_z,
      global_floor = global_floor,
      floor_gap = ints[, "b.mean"] - global_floor,
      peakiness = sq$peakiness,
      solidity = sq$solidity,
      adaptive_z_min = contrast_z_min,
      pass_area =
        is.finite(sh[, "s.area"]) &
        sh[, "s.area"] >= GREEN_MIN_AREA &
        sh[, "s.area"] <= GREEN_MAX_AREA,
      pass_global_gate =
        is.finite(ints[, "b.mean"]) &
        ints[, "b.mean"] >= global_floor,
      pass_eccentricity =
        is.finite(mom[, "m.eccentricity"]) &
        mom[, "m.eccentricity"] <= GREEN_MAX_ECCENTRICITY,
      pass_local_contrast =
        is.finite(lc$contrast_z) &
        lc$contrast_z >= contrast_z_min,
      stringsAsFactors = FALSE
    )
    
    checks$kept <- with(
      checks,
      pass_area &
        pass_global_gate &
        pass_eccentricity &
        pass_local_contrast
    )
    
    checks$reason <- vapply(seq_len(nrow(checks)), function(i) {
      failures <- c(
        if (!checks$pass_area[i]) "area",
        if (!checks$pass_global_gate[i]) "global gate",
        if (!checks$pass_eccentricity[i]) "elongated/streak-like",
        if (!checks$pass_local_contrast[i]) {
          "not significantly brighter than local background (adaptive)"
        }
      )
      
      if (!length(failures)) {
        "KEPT"
      } else {
        paste(failures, collapse = "; ")
      }
    }, character(1))
    
    final_lab <- safe_rmobj(
      lab, checks$candidate[!checks$kept]
    )
  }
  
  list(
    candidate_lab = lab,
    final_lab = final_lab,
    table = checks,
    dog = dog,
    threshold = thr,
    global_floor = global_floor
  )
}

## ---- Segmentation sanity check ----

.dummy_run_check <- function() {
  dummy <- matrix(runif(40 * 40, 0.01, 0.02), 40, 40)
  dummy[15:20, 15:20] <- 0.08
  
  res <- diagnose_green(
    list(G_seg = dummy),
    contrast_z_min = LOCAL_CONTRAST_Z_MIN
  )
  
  needed_cols <- c(
    "pass_local_contrast", "global_floor", "floor_gap",
    "peakiness", "solidity", "adaptive_z_min"
  )
  
  present <- needed_cols %in% names(res$table)
  names(present) <- needed_cols
  
  list(res = res, present = present)
}

.chk <- .dummy_run_check()

cat("\n[SANITY CHECK] diagnose_green() column check:\n")
print(.chk$present)

if (!all(.chk$present)) {
  stop("diagnose_green() is missing required diagnostic columns.")
}
if ("pass_edge_margin" %in% names(.chk$res$table)) {
  stop("Edge-margin gate is unexpectedly present.")
}

cat(
  "[SANITY CHECK] PASSED: adaptive contrast active; ",
  "edge gate removed; peakiness/solidity diagnostic-only.\n\n",
  sep = ""
)

rm(.chk, .dummy_run_check)

run_pair <- function(red_path, green_path, bleed_coeff = 0,
                     contrast_z_min = LOCAL_CONTRAST_Z_MIN) {
  R <- read_single_tiff(red_path, "RED / PGL-1")
  G <- read_single_tiff(green_path, "GREEN / MUT-16")
  
  if (!identical(dim(R), dim(G))) {
    stop(
      "Red and green dimensions differ: ",
      red_path, " / ", green_path,
      call. = FALSE
    )
  }
  
  R_blurred <- safe_gblur_m(R, sigma = 3)
  G_seg <- pmax(G - bleed_coeff * R_blurred, 0)
  rlab <- segment_red_unified(R)
  
  green_diagnostic <- diagnose_green(
    list(G_seg = G_seg, bleed_coeff = bleed_coeff),
    contrast_z_min = contrast_z_min
  )
  
  glab <- green_diagnostic$final_lab
  
  list(
    R = R, G = G, G_seg = G_seg,
    rlab = rlab, glab = glab,
    green_diagnostic = green_diagnostic
  )
}

mask_edge <- function(lab) {
  mask <- (lab > 0L) * 1.0
  brush <- safe_brush(3, min(dim(mask)))
  (mask - safe_erode_m(mask, brush)) > 0
}

display_scale <- function(m) {
  limits <- range(m)
  
  if (limits[2] <= limits[1]) {
    return(matrix(0, nrow(m), ncol(m)))
  }
  
  (m - limits[1]) / (limits[2] - limits[1])
}

make_gray_overlay <- function(view_m, lab, edge_color) {
  edge <- mask_edge(lab)
  channels <- list(view_m, view_m, view_m)
  
  for (k in 1:3) {
    channels[[k]][edge] <- edge_color[k]
  }
  
  rgbImage(
    eb(channels[[1]]),
    eb(channels[[2]]),
    eb(channels[[3]])
  )
}

save_qc_overlay <- function(R, G, rlab, glab,
                            out_path_red, out_path_green) {
  red_view <- display_scale(R)
  green_view <- display_scale(G)
  
  red_overlay <- make_gray_overlay(
    red_view, rlab, edge_color = c(0, 1, 1)
  )
  green_overlay <- make_gray_overlay(
    green_view, glab, edge_color = c(1, 0, 1)
  )
  
  writeImage(red_overlay, out_path_red, quality = 100)
  writeImage(green_overlay, out_path_green, quality = 100)
}

## ====
## SECTION B: FEATURES, MASKING AND AGGREGATION
## ====

nn_cv <- function(cen) {
  if (is.null(cen) || nrow(cen) < 4) return(NA_real_)
  
  d <- as.matrix(dist(cen))
  diag(d) <- NA_real_
  nn <- apply(d, 1, min, na.rm = TRUE)
  mn <- mean(nn)
  
  if (!is.finite(mn) || mn <= 0) return(NA_real_)
  
  ans <- sd(nn) / mn
  if (is.finite(ans)) ans else NA_real_
}

radial_cv <- function(cen, img_nr, img_nc) {
  if (is.null(cen) || nrow(cen) < 3) return(NA_real_)
  
  dists <- sqrt(
    (cen[, 1] - img_nc / 2)^2 +
      (cen[, 2] - img_nr / 2)^2
  )
  
  mn <- mean(dists)
  if (!is.finite(mn) || mn < 1e-6) return(NA_real_)
  
  ans <- sd(dists) / mn
  if (is.finite(ans)) ans else NA_real_
}

peak_coords <- function(lab_m, ref_m) {
  nobj <- max(lab_m)
  if (nobj == 0) return(NULL)
  
  do.call(rbind, lapply(seq_len(nobj), function(i) {
    mask <- ref_m
    mask[lab_m != as.integer(i)] <- 0
    idx <- which.max(mask)
    
    c(
      m.cx = ((idx - 1L) %/% nrow(ref_m)) + 1L,
      m.cy = ((idx - 1L) %% nrow(ref_m)) + 1L
    )
  }))
}

feature_cols <- c(
  "red_n", "red_size", "red_elong",
  "green_n", "green_size", "green_elong", "green_ring",
  "red_ring", "red_radial_cv",
  "green_app_net", "red_app_net_norm", "size_ratio"
)

count_cols <- c("red_n", "green_n")
shape_cols <- setdiff(feature_cols, count_cols)

foci_dependent <- c(
  "green_size", "green_elong", "green_app_net"
)

mask_clock_features <- function(d) {
  stopifnot(all(feature_cols %in% names(d)))
  
  for (f in feature_cols) {
    if (!is.numeric(d[[f]])) {
      stop("Non-numeric feature: ", f)
    }
    d[[f]][!is.finite(d[[f]])] <- NA_real_
  }
  
  no_foci <- !is.na(d$green_n) & d$green_n == 0
  
  for (f in foci_dependent) {
    d[[f]][no_foci] <- NA_real_
  }
  
  ## Counts and all rows are retained.
  d
}

aggregate_per_worm <- function(per_image) {
  if (is.null(per_image) || nrow(per_image) == 0L) {
    stop("Cannot aggregate an empty per-image table.")
  }
  if (anyNA(per_image$worm_id)) {
    stop("Missing worm_id in per-image table.")
  }
  
  out <- do.call(
    rbind,
    lapply(split(per_image, per_image$worm_id), function(df) {
      if (length(unique(df$day)) != 1L ||
          length(unique(df$condition)) != 1L) {
        stop("Inconsistent day or condition within worm: ", df$worm_id[1])
      }
      
      cnt_vals <- vapply(
        count_cols,
        function(col) finite_mean(df[[col]]),
        numeric(1)
      )
      
      shp_vals <- vapply(
        shape_cols,
        function(col) finite_median(df[[col]]),
        numeric(1)
      )
      
      agg <- c(cnt_vals, shp_vals)[feature_cols]
      
      data.frame(
        worm_id = df$worm_id[1],
        condition = df$condition[1],
        day = df$day[1],
        n_nuclei = nrow(df),
        as.data.frame(t(agg)),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  clean_numeric_nonfinite(out)
}

write_feature_tables <- function(per_image_raw, out_dir,
                                 prefix = "Batch1", tag = "") {
  if (is.null(per_image_raw) || nrow(per_image_raw) == 0L) {
    stop("No successfully extracted images to write.")
  }
  
  ## Preserve original zero-coded conditional measurements.
  ## Only non-finite numeric entries are normalized to NA.
  raw <- clean_numeric_nonfinite(per_image_raw)
  raw_worm <- aggregate_per_worm(raw)
  
  primary <- mask_clock_features(raw)
  primary_worm <- aggregate_per_worm(primary)
  
  write.csv(
    raw,
    file.path(out_dir, paste0(prefix, "_per_image_zerosIncluded", tag, ".csv")),
    row.names = FALSE, na = "NA"
  )
  write.csv(
    raw_worm,
    file.path(out_dir, paste0(prefix, "_per_worm_zerosIncluded", tag, ".csv")),
    row.names = FALSE, na = "NA"
  )
  write.csv(
    primary,
    file.path(out_dir, paste0(prefix, "_per_image", tag, ".csv")),
    row.names = FALSE, na = "NA"
  )
  write.csv(
    primary_worm,
    file.path(out_dir, paste0(prefix, "_per_worm", tag, ".csv")),
    row.names = FALSE, na = "NA"
  )
  
  no_foci <- !is.na(primary$green_n) & primary$green_n == 0
  
  stopifnot(
    identical(raw$green_n, primary$green_n),
    nrow(raw) == nrow(primary),
    all(is.na(as.matrix(
      primary[no_foci, foci_dependent, drop = FALSE]
    )))
  )
  
  missingness <- data.frame(
    trait = feature_cols,
    n_worms = nrow(primary_worm),
    n_missing = vapply(
      primary_worm[, feature_cols, drop = FALSE],
      function(x) sum(!is.finite(x)),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  missingness$pct_missing <- 100 *
    missingness$n_missing / missingness$n_worms
  
  write.csv(
    missingness,
    file.path(out_dir, paste0(prefix, "_per_worm_missingness", tag, ".csv")),
    row.names = FALSE
  )
  
  cat(sprintf(
    "\n[%s] %d images; %d worms; %d no-foci images retained.\n",
    prefix, nrow(primary), nrow(primary_worm), sum(no_foci)
  ))
  
  list(
    per_image = primary,
    per_worm = primary_worm,
    per_image_zerosIncluded = raw,
    per_worm_zerosIncluded = raw_worm
  )
}

## ====
## SECTION B2: TRAINING-ONLY CLOCK PREPROCESSING
## ====

clock_matrix <- function(d) {
  stopifnot(all(feature_cols %in% names(d)))
  
  if (!all(vapply(
    d[, feature_cols, drop = FALSE],
    is.numeric, logical(1)
  ))) {
    stop("All clock features must be numeric.")
  }
  
  X <- as.matrix(d[, feature_cols, drop = FALSE])
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- NA_real_
  X
}

clock_preprocess <- function(X) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- NA_real_
  
  med <- apply(X, 2, finite_median)
  all_missing <- !is.finite(med)
  
  ## 
  med[all_missing] <- 0
  
  for (j in seq_len(ncol(X))) {
    X[is.na(X[, j]), j] <- med[j]
  }
  
  sds <- apply(X, 2, finite_sd)
  disabled <- all_missing | !is.finite(sds) | sds == 0
  
  X[, disabled] <- 0
  sds[disabled] <- 0
  
  stopifnot(all(is.finite(X)))
  
  list(
    X = X,
    med = med,
    sds = sds,
    all_missing = all_missing,
    disabled = disabled
  )
}

clock_apply <- function(X, prep) {
  X <- as.matrix(X)
  
  if (!all(names(prep$med) %in% colnames(X))) {
    stop("Prediction matrix is missing required feature columns.")
  }
  
  X <- X[, names(prep$med), drop = FALSE]
  X[!is.finite(X)] <- NA_real_
  
  for (j in seq_len(ncol(X))) {
    X[is.na(X[, j]), j] <- prep$med[j]
  }
  
  ## A feature unavailable or constant during training cannot
  ## acquire an effect because it varies in held-out data.
  X[, prep$disabled] <- 0
  
  stopifnot(all(is.finite(X)))
  X
}

clock_path <- function(X, y, lambda) {
  has_variation <- any(apply(X, 2, finite_sd) > 0, na.rm = TRUE)
  
  if (length(unique(y)) < 2L || !has_variation) {
    return(NULL)
  }
  
  glmnet::glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = CLOCK_ALPHA,
    lambda = lambda,
    standardize = TRUE
  )
}

clock_path_predict <- function(fit, X, lambda, y_train) {
  if (is.null(fit)) {
    return(matrix(
      mean(y_train),
      nrow = nrow(X),
      ncol = length(lambda)
    ))
  }
  
  as.matrix(predict(fit, newx = X, s = lambda))
}

fit_clock_missing <- function(X, y, nfolds = CLOCK_INNER_FOLDS) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  
  stopifnot(
    length(y) == n,
    all(is.finite(y)),
    n >= 4L,
    !is.null(colnames(X))
  )
  
  k <- min(as.integer(nfolds), n)
  stopifnot(k >= 3L)
  
  foldid <- sample(rep(seq_len(k), length.out = n))
  loss <- matrix(NA_real_, n, length(CLOCK_LAMBDA))
  
  ## 
  for (fold in seq_len(k)) {
    test <- which(foldid == fold)
    train <- which(foldid != fold)
    
    prep <- clock_preprocess(X[train, , drop = FALSE])
    X_test <- clock_apply(X[test, , drop = FALSE], prep)
    
    fit <- clock_path(prep$X, y[train], CLOCK_LAMBDA)
    pred <- clock_path_predict(
      fit, X_test, CLOCK_LAMBDA, y[train]
    )
    
    loss[test, ] <- sweep(pred, 1, y[test], "-")^2
  }
  
  stopifnot(all(is.finite(loss)))
  
  cvm <- colMeans(loss)
  best <- which.min(cvm)
  
  ## 
  prep <- clock_preprocess(X)
  fit <- clock_path(prep$X, y, CLOCK_LAMBDA)
  
  beta <- setNames(rep(0, ncol(X)), colnames(X))
  intercept <- mean(y)
  
  if (!is.null(fit)) {
    cc <- as.matrix(coef(fit, s = CLOCK_LAMBDA[best]))
    beta[] <- cc[colnames(X), 1]
    intercept <- cc["(Intercept)", 1]
  }
  
  stopifnot(
    all(is.finite(beta)),
    is.finite(intercept)
  )
  
  list(
    fit = fit,
    prep = prep,
    coefficient = beta,
    intercept = intercept,
    lambda.min = CLOCK_LAMBDA[best],
    lambda = CLOCK_LAMBDA,
    cvm = cvm,
    foldid = foldid,
    boundary = best %in% c(1L, length(CLOCK_LAMBDA))
  )
}

predict_clock_missing <- function(object, X) {
  X <- clock_apply(X, object$prep)
  out <- as.numeric(object$intercept + X %*% object$coefficient)
  stopifnot(all(is.finite(out)))
  out
}

clock_loo <- function(per_worm, seed = CLOCK_SEED) {
  X <- clock_matrix(per_worm)
  y <- as.numeric(per_worm$day)
  n <- nrow(X)
  
  stopifnot(
    n >= 5L,
    all(is.finite(y)),
    !anyNA(per_worm$worm_id),
    !anyDuplicated(per_worm$worm_id)
  )
  
  set.seed(seed)
  
  preds <- rep(NA_real_, n)
  coefficients <- matrix(
    0,
    nrow = n,
    ncol = length(feature_cols),
    dimnames = list(per_worm$worm_id, feature_cols)
  )
  
  selected_lambda <- rep(NA_real_, n)
  boundary <- rep(FALSE, n)
  
  for (i in seq_len(n)) {
    fit <- fit_clock_missing(
      X[-i, , drop = FALSE],
      y[-i]
    )
    
    preds[i] <- predict_clock_missing(
      fit, X[i, , drop = FALSE]
    )
    
    coefficients[i, ] <- fit$coefficient[feature_cols]
    selected_lambda[i] <- fit$lambda.min
    boundary[i] <- fit$boundary
  }
  
  ## Never silently drop failed predictions from performance.
  stopifnot(all(is.finite(preds)))
  
  if (any(boundary)) {
    warning(
      sum(boundary),
      " outer fits selected a lambda-grid boundary (seed ",
      seed, "). Inspect the tuning range."
    )
  }
  
  list(
    predictions = preds,
    coefficients = coefficients,
    selected_lambda = selected_lambda,
    boundary = boundary
  )
}

## ---- Checks executed before image processing ----

local({
  d <- as.data.frame(matrix(1, 2, length(feature_cols)))
  names(d) <- feature_cols
  d$green_n <- c(0, 1)
  
  z <- mask_clock_features(d)
  
  stopifnot(
    z$green_n[1] == 0,
    all(is.na(as.matrix(z[1, foci_dependent, drop = FALSE]))),
    z$green_size[2] == 1
  )
  
  train <- cbind(
    a = c(1, NA, 3),
    b = c(NA_real_, NaN, NA_real_),
    c = c(2, 2, 2)
  )
  
  prep <- clock_preprocess(train)
  test <- clock_apply(
    cbind(a = NA_real_, b = 100, c = 200),
    prep
  )
  
  stopifnot(
    test[1, "a"] == 2,
    test[1, "b"] == 0,
    test[1, "c"] == 0
  )
  
  cat("[SANITY CHECK] Conditional masking and imputation checks passed.\n")
})

## ====
## SECTION C: FEATURE EXTRACTION + QC OVERLAYS
## ====

extract_features_real <- function(red_path, green_path, meta,
                                  overlay_dir = NULL,
                                  contrast_z_min = LOCAL_CONTRAST_Z_MIN) {
  p <- tryCatch(
    run_pair(
      red_path, green_path,
      bleed_coeff = BLEED_COEFF,
      contrast_z_min = contrast_z_min
    ),
    error = function(e) {
      warning(e$message)
      NULL
    }
  )
  
  if (is.null(p)) return(NULL)
  
  R <- p$R
  G <- p$G
  rlab <- p$rlab
  glab <- p$glab
  
  nr <- nrow(R)
  nc <- ncol(R)
  npix <- nr * nc
  rN <- as.integer(max(rlab))
  gN <- as.integer(max(glab))
  
  if (!is.null(overlay_dir)) {
    tryCatch({
      save_qc_overlay(
        R, G, rlab, glab,
        out_path_red = file.path(
          overlay_dir, paste0(meta$base_id, "_red.png")
        ),
        out_path_green = file.path(
          overlay_dir, paste0(meta$base_id, "_green.png")
        )
      )
    }, error = function(e) {
      warning(
        "Overlay save failed for ", meta$base_id, ": ", e$message
      )
    })
  }
  
  if (rN > 0) {
    sh_r <- tryCatch(
      computeFeatures.shape(eb(rlab)),
      error = function(e) NULL
    )
    mom_r <- tryCatch(
      computeFeatures.moment(eb(rlab)),
      error = function(e) NULL
    )
    
    red_n <- rN
    red_size <- if (!is.null(sh_r)) {
      median(sh_r[, "s.area"])
    } else {
      NA_real_
    }
    red_cov <- sum(rlab > 0L) / npix * 100
    red_elong <- if (!is.null(mom_r)) {
      median(mom_r[, "m.eccentricity"])
    } else {
      NA_real_
    }
    
    pk_r <- peak_coords(rlab, R)
    red_ring <- nn_cv(pk_r)
    red_radial_cv <- radial_cv(pk_r, nr, nc)
  } else {
    red_n <- 0L
    red_size <- 0
    red_cov <- 0
    red_elong <- 0
    pk_r <- NULL
    red_ring <- NA_real_
    red_radial_cv <- NA_real_
  }
  
  if (gN > 0) {
    sh_g <- tryCatch(
      computeFeatures.shape(eb(glab)),
      error = function(e) NULL
    )
    mom_g <- tryCatch(
      computeFeatures.moment(eb(glab)),
      error = function(e) NULL
    )
    
    green_n <- gN
    green_size <- if (!is.null(sh_g)) {
      median(sh_g[, "s.area"])
    } else {
      NA_real_
    }
    green_elong <- if (!is.null(mom_g)) {
      median(mom_g[, "m.eccentricity"])
    } else {
      NA_real_
    }
    
    pk_g <- peak_coords(glab, p$G_seg)
    green_ring <- nn_cv(pk_g)
  } else {
    ## Retained here for the raw zerosIncluded tables.
    ## 
    green_n <- 0L
    green_size <- 0
    green_elong <- 0
    pk_g <- NULL
    green_ring <- NA_real_
  }
  
  green_app_net <- 0
  red_app_net_norm <- NA_real_
  coloc_nucleus <- NA_real_
  
  if (rN > 0 && gN > 0) {
    min_dim <- min(nr, nc)
    br7 <- safe_brush(7, min_dim)
    br5 <- safe_brush(5, min_dim)
    
    rdil <- as_mat(dilate(eb((rlab > 0L) * 1.0), br7)) > 0
    gdil <- as_mat(dilate(eb((glab > 0L) * 1.0), br7)) > 0
    
    tissue <- (
      (R > mean(R) + 0.5 * sd(R)) |
        (G > mean(G) + 0.5 * sd(G))
    )
    tissue <- as_mat(dilate(eb(tissue * 1.0), br5)) > 0
    
    chance_rg <- sum(rdil & tissue) / max(sum(tissue), 1) * 100
    chance_gr <- sum(gdil & tissue) / max(sum(tissue), 1) * 100
    
    touch_g <- if (!is.null(pk_g)) {
      sum(sapply(seq_len(nrow(pk_g)), function(i) {
        rdil[pk_g[i, "m.cy"], pk_g[i, "m.cx"]]
      }))
    } else {
      0L
    }
    
    green_app_net <- touch_g / gN * 100 - chance_rg
    
    touch_r <- if (!is.null(pk_r)) {
      sum(sapply(seq_len(nrow(pk_r)), function(i) {
        gdil[pk_r[i, "m.cy"], pk_r[i, "m.cx"]]
      }))
    } else {
      0L
    }
    
    red_app_net_norm <- (
      touch_r / rN * 100 - chance_gr
    ) / max(gN, 1)
    
    r_thresh <- mean(R) + 0.5 * sd(R)
    g_thresh <- mean(G) + 0.5 * sd(G)
    
    both_high <- (R > r_thresh) & (G > g_thresh)
    either_hi <- (R > r_thresh) | (G > g_thresh)
    
    coloc_nucleus <- sum(both_high) / max(sum(either_hi), 1) * 100
  }
  
  size_ratio <- if (
    rN > 0 && gN > 0 &&
    is.finite(green_size) && green_size > 0
  ) {
    red_size / green_size
  } else {
    NA_real_
  }
  
  out <- data.frame(
    file = meta$base_id,
    condition = meta$condition,
    day = meta$day,
    worm_id = meta$worm_id,
    nucleus = meta$nucleus,
    red_n = red_n,
    red_size = red_size,
    red_cov = red_cov,
    red_elong = red_elong,
    green_n = green_n,
    green_size = green_size,
    green_elong = green_elong,
    green_ring = green_ring,
    red_ring = red_ring,
    red_radial_cv = red_radial_cv,
    green_app_net = green_app_net,
    red_app_net_norm = red_app_net_norm,
    coloc_nucleus = coloc_nucleus,
    size_ratio = size_ratio,
    stringsAsFactors = FALSE
  )
  
  clean_numeric_nonfinite(out)
}

## ====
## SECTION D: FOLDER DISCOVERY
## ====

find_nuclei_flat <- function(folder, channel_letter) {
  files <- list.files(
    folder, pattern = "\\.(tif|tiff)$",
    full.names = TRUE, ignore.case = TRUE, recursive = FALSE
  )
  files <- files[!dir.exists(files)]
  
  if (!length(files)) {
    return(data.frame(nucleus = numeric(0), path = character(0)))
  }
  
  pattern <- paste0("_", channel_letter, "n?([0-9]+)\\.(tif|tiff)$")
  
  ids <- vapply(basename(files), function(name) {
    m <- regmatches(
      name, regexec(pattern, name, ignore.case = TRUE)
    )[[1]]
    
    if (length(m) < 2) return(NA_real_)
    as.numeric(m[2])
  }, numeric(1))
  
  valid <- is.finite(ids) & ids >= 1
  
  data.frame(
    nucleus = ids[valid],
    path = files[valid],
    stringsAsFactors = FALSE
  )
}

extract_day <- function(worm_folder) {
  hit <- regmatches(
    worm_folder,
    regexpr("D([0-9]+)", worm_folder, ignore.case = TRUE)
  )
  
  if (length(hit) == 0 || hit == "") return(NA_integer_)
  as.integer(sub("D", "", hit, ignore.case = TRUE))
}

extract_worm_num <- function(worm_folder) {
  wname <- basename(worm_folder)
  m <- regmatches(
    wname,
    regexpr("W([0-9]+)$", wname, ignore.case = TRUE)
  )
  
  if (length(m) == 0 || m == "") return(NA_character_)
  sub("W", "", m, ignore.case = TRUE)
}

collect_worm_pairs <- function(worm_folder) {
  day <- extract_day(worm_folder)
  worm_no <- extract_worm_num(worm_folder)
  
  if (is.na(day) || is.na(worm_no)) return(NULL)
  
  red <- find_nuclei_flat(worm_folder, "R")
  green <- find_nuclei_flat(worm_folder, "G")
  ids <- intersect(red$nucleus, green$nucleus)
  
  if (length(ids) == 0) return(NULL)
  
  data.frame(
    worm_id = paste0("adlib_D", day, "_W", worm_no),
    condition = "adlib",
    day = as.integer(day),
    worm = worm_no,
    nucleus = ids,
    base_id = paste0("AL_B1_D", day, "_W", worm_no, "_n", ids),
    red_path = red$path[match(ids, red$nucleus)],
    green_path = green$path[match(ids, green$nucleus)],
    stringsAsFactors = FALSE
  )
}

find_nuclei_flat_IF <- function(folder, channel_letter) {
  files <- list.files(
    folder, pattern = "\\.(tif|tiff)$",
    full.names = TRUE, ignore.case = TRUE, recursive = FALSE
  )
  files <- files[!dir.exists(files)]
  
  if (!length(files)) {
    return(data.frame(nucleus = numeric(0), path = character(0)))
  }
  
  pattern <- paste0(channel_letter, "n([0-9]+)\\.(tif|tiff)$")
  
  ids <- vapply(basename(files), function(name) {
    m <- regmatches(
      name, regexec(pattern, name, ignore.case = TRUE)
    )[[1]]
    
    if (length(m) < 2) return(NA_real_)
    as.numeric(m[2])
  }, numeric(1))
  
  valid <- is.finite(ids) & ids >= 1
  
  data.frame(
    nucleus = ids[valid],
    path = files[valid],
    stringsAsFactors = FALSE
  )
}

collect_worm_pairs_IF <- function(worm_folder) {
  wname <- basename(worm_folder)
  
  if (!grepl("^IF_B1_D[0-9]+_W[0-9]+$", wname)) {
    return(NULL)
  }
  
  day <- extract_day(worm_folder)
  worm_no <- extract_worm_num(worm_folder)
  
  if (is.na(day) || is.na(worm_no)) return(NULL)
  
  red <- find_nuclei_flat_IF(worm_folder, "R")
  green <- find_nuclei_flat_IF(worm_folder, "G")
  ids <- intersect(red$nucleus, green$nucleus)
  
  if (length(ids) == 0) return(NULL)
  
  data.frame(
    worm_id = paste0("IF_D", day, "_W", worm_no),
    condition = "IF",
    day = as.integer(day),
    worm = worm_no,
    nucleus = ids,
    base_id = paste0("IF_B1_D", day, "_W", worm_no, "_n", ids),
    red_path = red$path[match(ids, red$nucleus)],
    green_path = green$path[match(ids, green$nucleus)],
    stringsAsFactors = FALSE
  )
}

## ====
## SECTION D0: ADAPTIVE NOISE SCALE
## ====

compute_image_noise_scale <- function(pair_tbl) {
  cat(
    "\n[ADAPTIVE] Measuring per-image local_bg_mad for",
    nrow(pair_tbl), "images...\n"
  )
  
  scales <- vector("list", nrow(pair_tbl))
  
  for (i in seq_len(nrow(pair_tbl))) {
    scales[[i]] <- tryCatch({
      p <- run_pair(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        bleed_coeff = BLEED_COEFF,
        contrast_z_min = LOCAL_CONTRAST_Z_MIN
      )
      
      tbl <- p$green_diagnostic$table
      med_mad <- finite_median(tbl$local_bg_mad)
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        img_bg_mad = med_mad,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("Noise measurement failed: ", pair_tbl$base_id[i], ": ", e$message)
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        img_bg_mad = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf(
        "  [%d/%d] %s\n",
        i, nrow(pair_tbl), pair_tbl$base_id[i]
      ))
    }
  }
  
  scale_tbl <- do.call(rbind, scales)
  batch_ref <- finite_median(scale_tbl$img_bg_mad)
  
  validate_noise_reference(batch_ref, "AL batch noise reference")
  
  scale_tbl$noise_ratio <- scale_tbl$img_bg_mad / batch_ref
  unavailable <- !is.finite(scale_tbl$noise_ratio)
  scale_tbl$noise_ratio[unavailable] <- 1.0
  
  if (any(unavailable)) {
    warning(
      sum(unavailable),
      " images have no usable local noise estimate; ",
      "using the original noise_ratio = 1 fallback."
    )
  }
  
  cat(sprintf(
    "[ADAPTIVE] Batch reference local_bg_mad = %.6g\n",
    batch_ref
  ))
  cat(sprintf(
    "[ADAPTIVE] noise_ratio range: %.3f to %.3f (median %.3f)\n",
    min(scale_tbl$noise_ratio),
    max(scale_tbl$noise_ratio),
    median(scale_tbl$noise_ratio)
  ))
  
  list(scale_tbl = scale_tbl, batch_ref = batch_ref)
}

compute_image_noise_scale_IF <- function(pair_tbl, batch_ref = NULL) {
  cat(
    "\n[ADAPTIVE-IF] Measuring per-image background noise for",
    nrow(pair_tbl), "IF images...\n"
  )
  
  scales <- vector("list", nrow(pair_tbl))
  
  for (i in seq_len(nrow(pair_tbl))) {
    scales[[i]] <- tryCatch({
      p <- run_pair(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        bleed_coeff = BLEED_COEFF,
        contrast_z_min = LOCAL_CONTRAST_Z_MIN
      )
      
      tbl <- p$green_diagnostic$table
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        img_bg_mad = finite_median(tbl$local_bg_mad),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning(
        "IF noise measurement failed: ",
        pair_tbl$base_id[i], ": ", e$message
      )
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        img_bg_mad = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf(
        "  [%d/%d] %s\n",
        i, nrow(pair_tbl), pair_tbl$base_id[i]
      ))
    }
  }
  
  scale_tbl <- do.call(rbind, scales)
  own_ref <- finite_median(scale_tbl$img_bg_mad)
  used_ref <- if (is.null(batch_ref)) own_ref else batch_ref
  
  validate_noise_reference(used_ref, "IF segmentation noise reference")
  
  scale_tbl$noise_ratio <- scale_tbl$img_bg_mad / used_ref
  unavailable <- !is.finite(scale_tbl$noise_ratio)
  scale_tbl$noise_ratio[unavailable] <- 1.0
  
  if (any(unavailable)) {
    warning(
      sum(unavailable),
      " IF images have no usable noise estimate; using noise_ratio = 1."
    )
  }
  
  cat(sprintf(
    "[ADAPTIVE-IF] Reference used = %.6g (%s) | IF-own median = %.6g\n",
    used_ref,
    if (is.null(batch_ref)) "IF-own median" else "AL reference",
    own_ref
  ))
  cat(sprintf(
    "[ADAPTIVE-IF] noise_ratio range: %.3f to %.3f (median %.3f)\n",
    min(scale_tbl$noise_ratio),
    max(scale_tbl$noise_ratio),
    median(scale_tbl$noise_ratio)
  ))
  
  list(scale_tbl = scale_tbl, batch_ref = used_ref)
}

## ====
## SECTION D2: IMAGE-LEVEL GREEN QC
## ====

diagnose_green_detection <- function(pair_tbl, out_dir,
                                     green_n_threshold = GREEN_N_THRESHOLD) {
  cat(
    "\n[QC] Running green-detection diagnostics on",
    nrow(pair_tbl), "images...\n"
  )
  
  diag_rows <- vector("list", nrow(pair_tbl))
  
  for (i in seq_len(nrow(pair_tbl))) {
    z_min_i <- pair_tbl$adaptive_z_min[i]
    
    res <- tryCatch({
      p <- run_pair(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        bleed_coeff = BLEED_COEFF,
        contrast_z_min = z_min_i
      )
      
      gN <- as.integer(max(p$glab))
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        worm_id = pair_tbl$worm_id[i],
        day = pair_tbl$day[i],
        adaptive_z_min = z_min_i,
        green_n = gN,
        below_threshold = gN < green_n_threshold,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("QC failed: ", pair_tbl$base_id[i], ": ", e$message)
      
      data.frame(
        base_id = pair_tbl$base_id[i],
        worm_id = pair_tbl$worm_id[i],
        day = pair_tbl$day[i],
        adaptive_z_min = z_min_i,
        green_n = NA_integer_,
        below_threshold = NA,
        stringsAsFactors = FALSE
      )
    })
    
    diag_rows[[i]] <- res
    
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf(
        "  [%d/%d] %s green_n=%s (z_min=%.3f)\n",
        i, nrow(pair_tbl), pair_tbl$base_id[i], res$green_n, z_min_i
      ))
    }
  }
  
  diag_tbl <- do.call(rbind, diag_rows)
  
  write.csv(
    diag_tbl,
    file.path(out_dir, "Batch1_green_detection_QC.csv"),
    row.names = FALSE
  )
  
  n_total <- nrow(diag_tbl)
  n_flagged <- sum(diag_tbl$below_threshold, na.rm = TRUE)
  n_zero <- sum(diag_tbl$green_n == 0, na.rm = TRUE)
  n_failed <- sum(is.na(diag_tbl$green_n))
  
  cat(sprintf(
    "\n[QC SUMMARY] %d / %d images have green_n < %d.\n",
    n_flagged, n_total, green_n_threshold
  ))
  cat(sprintf(
    "  %d images have no detected green foci; %d segmentation failures.\n",
    n_zero, n_failed
  ))
  
  by_day <- do.call(rbind, lapply(
    split(diag_tbl, diag_tbl$day),
    function(df) {
      valid <- !is.na(df$below_threshold)
      
      data.frame(
        day = df$day[1],
        n_images = nrow(df),
        n_evaluated = sum(valid),
        n_failed = sum(!valid),
        pct_flagged = if (any(valid)) {
          100 * mean(df$below_threshold[valid])
        } else {
          NA_real_
        }
      )
    }
  ))
  
  worm_agg <- do.call(rbind, lapply(
    split(diag_tbl, diag_tbl$worm_id),
    function(df) {
      valid <- !is.na(df$below_threshold)
      n_flag <- sum(df$below_threshold, na.rm = TRUE)
      
      data.frame(
        worm_id = df$worm_id[1],
        n_images = nrow(df),
        n_evaluated = sum(valid),
        n_failed = sum(!valid),
        n_flagged = n_flag,
        pct_flagged = if (any(valid)) {
          round(100 * n_flag / sum(valid), 1)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  ))
  
  worm_agg <- worm_agg[
    order(-worm_agg$pct_flagged, na.last = TRUE),
    ,
    drop = FALSE
  ]
  
  write.csv(
    worm_agg,
    file.path(out_dir, "Batch1_green_detection_by_worm.csv"),
    row.names = FALSE
  )
  
  cat("\n[QC] Percent flagged among evaluated images by day:\n")
  print(by_day)
  
  cat("\n[QC] Most affected worms:\n")
  print(head(worm_agg, 10))
  
  list(
    diag_tbl = diag_tbl,
    by_day = by_day,
    by_worm = worm_agg
  )
}

## ====
## SECTION D3: CANDIDATE-LEVEL GREEN QC
## ====

dump_green_candidate_qc <- function(pair_tbl, out_dir) {
  cat(
    "\n[QC] Dumping per-candidate diagnostics for",
    nrow(pair_tbl), "images...\n"
  )
  
  ## lapply keeps NULL entries without altering list indices.
  ## A zero-candidate image must not return from this whole function.
  rows <- lapply(seq_len(nrow(pair_tbl)), function(i) {
    tryCatch({
      p <- run_pair(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        bleed_coeff = BLEED_COEFF,
        contrast_z_min = pair_tbl$adaptive_z_min[i]
      )
      
      tbl <- p$green_diagnostic$table
      nr <- nrow(tbl)
      
      tbl$base_id <- rep(pair_tbl$base_id[i], nr)
      tbl$worm_id <- rep(pair_tbl$worm_id[i], nr)
      tbl$day <- rep(pair_tbl$day[i], nr)
      tbl$noise_ratio <- rep(pair_tbl$noise_ratio[i], nr)
      tbl
    }, error = function(e) {
      warning(
        "Candidate QC failed: ",
        pair_tbl$base_id[i], ": ", e$message
      )
      NULL
    })
  })
  
  valid_rows <- Filter(Negate(is.null), rows)
  
  if (length(valid_rows)) {
    full <- do.call(rbind, valid_rows)
  } else {
    full <- empty_candidate_table()
    full$base_id <- character(0)
    full$worm_id <- character(0)
    full$day <- integer(0)
    full$noise_ratio <- numeric(0)
  }
  
  full <- clean_numeric_nonfinite(full)
  
  write.csv(
    full,
    file.path(out_dir, "Batch1_green_candidate_QC.csv"),
    row.names = FALSE
  )
  
  if (nrow(full) > 0L) {
    cat(sprintf(
      paste0(
        "\n[QC] %d candidates | %d kept | ",
        "%d eccentricity rejections | %d local-contrast rejections\n"
      ),
      nrow(full),
      sum(full$kept, na.rm = TRUE),
      sum(!full$pass_eccentricity, na.rm = TRUE),
      sum(!full$pass_local_contrast, na.rm = TRUE)
    ))
    
    kept <- which(full$kept)
    
    cat(sprintf(
      "[QC] Peakiness, kept: median=%.3f; IQR=[%.3f, %.3f]\n",
      finite_median(full$peakiness[kept]),
      finite_quantile(full$peakiness[kept], 0.25),
      finite_quantile(full$peakiness[kept], 0.75)
    ))
    cat(sprintf(
      "[QC] Solidity, kept: median=%.3f; IQR=[%.3f, %.3f]\n",
      finite_median(full$solidity[kept]),
      finite_quantile(full$solidity[kept], 0.25),
      finite_quantile(full$solidity[kept], 0.75)
    ))
    
    cat("\n[QC] Kept candidates by day:\n")
    print(aggregate(kept ~ day, data = full, FUN = sum))
  } else {
    cat("[QC] Candidate table is empty; column headers were saved.\n")
  }
  
  full
}

## ====
## SECTION E: AL BATCH PROCESSING
## ====

process_batch1 <- function(root_dir) {
  out_dir <- file.path(OUT_ROOT, "Batch1_results_D1to6")
  fig_dir <- file.path(out_dir, "figures")
  overlay_dir <- file.path(out_dir, "qc_overlays")
  
  for (d in c(out_dir, fig_dir, overlay_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  day_dirs <- list.dirs(root_dir, recursive = FALSE)
  day_dirs <- day_dirs[grepl("^B1_D[0-9]+$", basename(day_dirs))]
  
  if (!length(day_dirs)) {
    stop("No B1_D* folders found under: ", root_dir)
  }
  
  cat("Day folders found:", paste(basename(day_dirs), collapse = ", "), "\n")
  
  all_pairs <- list()
  
  for (dd in day_dirs) {
    worm_folders <- list.dirs(dd, recursive = FALSE)
    worm_folders <- worm_folders[
      grepl("^AL_B1_D[0-9]+_W[0-9]+$", basename(worm_folders))
    ]
    
    for (wf in worm_folders) {
      pr <- collect_worm_pairs(wf)
      
      if (!is.null(pr)) {
        all_pairs[[length(all_pairs) + 1L]] <- pr
      }
    }
  }
  
  pair_tbl <- do.call(rbind, all_pairs)
  
  if (is.null(pair_tbl) || nrow(pair_tbl) == 0L) {
    stop("No matched AL nucleus pairs found under: ", root_dir)
  }
  
  if (anyDuplicated(pair_tbl$base_id)) {
    stop("Duplicate AL base_id values detected.")
  }
  
  ## Every expected day (including D6) must have AL worms.
  days_found <- sort(unique(pair_tbl$day))
  days_missing <- setdiff(EXPECTED_DAYS, days_found)
  
  if (length(days_missing)) {
    stop(
      "No AL worms found for expected day(s): ",
      paste(days_missing, collapse = ", "),
      ". Check folder names (B1_D<day>/AL_B1_D<day>_W<n>)."
    )
  }
  
  cat(sprintf(
    "Matched %d nucleus pairs across %d AL worms (days: %s).\n",
    nrow(pair_tbl), length(unique(pair_tbl$worm_id)),
    paste(days_found, collapse = ", ")
  ))
  print(table(
    day = pair_tbl$day[!duplicated(pair_tbl$worm_id)]
  ))
  
  noise_info <- compute_image_noise_scale(pair_tbl)
  
  write.csv(
    noise_info$scale_tbl,
    file.path(out_dir, "Batch1_image_noise_scale.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(batch_ref = noise_info$batch_ref),
    file.path(out_dir, "Batch1_noise_reference.csv"),
    row.names = FALSE
  )
  
  pair_tbl <- merge(
    pair_tbl,
    noise_info$scale_tbl[, c("base_id", "noise_ratio")],
    by = "base_id", all.x = TRUE
  )
  
  pair_tbl$noise_ratio[!is.finite(pair_tbl$noise_ratio)] <- 1.0
  pair_tbl$adaptive_z_min <- LOCAL_CONTRAST_Z_MIN * pair_tbl$noise_ratio
  
  cat("[ADAPTIVE] Median effective z-min by day:\n")
  print(aggregate(
    adaptive_z_min ~ day, data = pair_tbl, FUN = median
  ))
  
  green_qc <- diagnose_green_detection(
    pair_tbl, out_dir,
    green_n_threshold = GREEN_N_THRESHOLD
  )
  candidate_qc <- dump_green_candidate_qc(pair_tbl, out_dir)
  
  cat("\n[Step] Extracting AL features and saving QC overlays...\n")
  
  feat_list <- lapply(seq_len(nrow(pair_tbl)), function(i) {
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf(
        "  [%d/%d] %s\n",
        i, nrow(pair_tbl), pair_tbl$base_id[i]
      ))
    }
    
    meta <- as.list(pair_tbl[
      i, c("condition", "day", "worm_id", "nucleus", "base_id")
    ])
    
    tryCatch(
      extract_features_real(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        meta,
        overlay_dir = overlay_dir,
        contrast_z_min = pair_tbl$adaptive_z_min[i]
      ),
      error = function(e) {
        warning("SKIP: ", pair_tbl$base_id[i], " - ", e$message)
        NULL
      }
    )
  })
  
  per_image_raw <- do.call(
    rbind, Filter(Negate(is.null), feat_list)
  )
  
  tables <- write_feature_tables(
    per_image_raw, out_dir, prefix = "Batch1"
  )
  
  list(
    per_image = tables$per_image,
    per_worm = tables$per_worm,
    per_image_zerosIncluded = tables$per_image_zerosIncluded,
    per_worm_zerosIncluded = tables$per_worm_zerosIncluded,
    out_dir = out_dir,
    fig_dir = fig_dir,
    overlay_dir = overlay_dir,
    green_qc = green_qc,
    candidate_qc = candidate_qc,
    pair_tbl = pair_tbl,
    batch_ref = noise_info$batch_ref
  )
}

## ====
## SECTION F: LOO-CV CLOCK + PERFORMANCE PLOTS
## ====

run_clock_LOOCV <- function(per_worm, out_dir, fig_dir) {
  y <- as.numeric(per_worm$day)
  loo <- clock_loo(per_worm, seed = CLOCK_SEED)
  preds <- loo$predictions
  
  result <- data.frame(
    worm_id = per_worm$worm_id,
    actual_day = y,
    predicted_day = preds,
    stringsAsFactors = FALSE
  )
  
  result$abs_err <- abs(result$predicted_day - result$actual_day)
  result$sq_err <- (result$predicted_day - result$actual_day)^2
  result$err <- result$predicted_day - result$actual_day
  
  write.csv(
    result,
    file.path(out_dir, "Batch1_LOOCV_predictions.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      worm_id = per_worm$worm_id,
      selected_lambda = loo$selected_lambda,
      lambda_boundary = loo$boundary
    ),
    file.path(out_dir, "Batch1_LOOCV_tuning.csv"),
    row.names = FALSE
  )
  
  rho <- safe_cor(y, preds, method = "spearman")
  
  ## Retains the original metric definition:
  ## squared Pearson correlation
  r2 <- safe_cor(y, preds)^2
  
  overall_mae <- mean(result$abs_err)
  overall_rmse <- sqrt(mean(result$sq_err))
  
  mae_by_day <- aggregate(
    cbind(MAE = abs_err, RMSE = sq_err, bias = err) ~ actual_day,
    data = result,
    FUN = mean
  )
  mae_by_day$RMSE <- sqrt(mae_by_day$RMSE)
  names(mae_by_day)[1] <- "day"
  
  write.csv(
    mae_by_day,
    file.path(out_dir, "Batch1_MAE_by_day.csv"),
    row.names = FALSE
  )
  
  cat(sprintf(
    paste0(
      "\nLOO-CV: n=%d | Spearman rho=%.3f | ",
      "Pearson r²=%.3f | MAE=%.3f d | RMSE=%.3f d\n"
    ),
    nrow(result), rho, r2, overall_mae, overall_rmse
  ))
  print(mae_by_day)
  
  day_levels <- sort(unique(result$actual_day))
  day_labels <- paste0("D", day_levels)
  result$actual_day_f <- factor(
    result$actual_day, levels = day_levels
  )
  
  p <- ggplot(
    result,
    aes(x = actual_day, y = predicted_day, colour = actual_day_f)
  ) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = "grey50", linewidth = LW
    ) +
    geom_jitter(
      width = 0.08, height = 0,
      size = PT_SZ, alpha = 0.85, stroke = 0.3
    ) +
    scale_colour_manual(values = col_day) +
    scale_x_continuous(
      breaks = day_levels,
      labels = day_labels,
      limits = c(min(day_levels) - 0.5, max(day_levels) + 0.5)
    ) +
    scale_y_continuous(
      breaks = day_levels,
      labels = day_labels,
      limits = c(min(day_levels) - 0.5, max(day_levels) + 0.5)
    ) +
    annotate(
      "text",
      x = min(day_levels),
      y = max(day_levels) + 0.3,
      hjust = 0, vjust = 1,
      label = sprintf(
        "MAE = %.3f d\nRMSE = %.3f d\nPearson %s = %.3f\n%s = %.3f",
        overall_mae, overall_rmse, SYM_R2, r2, SYM_RHO, rho
      )
    ) +
    labs(
      title = "GG-Clock — (LOO-CV)",
      x = "Chronological age (days)",
      y = "Predicted age (days)"
    ) +
    pub_theme +
    theme(legend.position = "none")
  
  save_fig(p, "Batch1_LOOCV_scatter", fig_dir, w = 5.5, h = 5.5)
  
  p_mae <- ggplot(
    mae_by_day,
    aes(
      x = factor(day, levels = day_levels, labels = day_labels),
      y = MAE,
      fill = factor(day, levels = day_levels)
    )
  ) +
    geom_col(width = 0.65, alpha = 0.9, colour = "grey20") +
    geom_hline(
      yintercept = overall_mae, linetype = "dashed",
      color = "firebrick", linewidth = LW
    ) +
    annotate(
      "text", x = Inf, y = overall_mae,
      label = sprintf("Overall MAE = %.2f", overall_mae),
      hjust = 1.05, vjust = -0.5,
      color = "firebrick"
    ) +
    scale_fill_manual(values = col_day, guide = "none") +
    labs(
      title = "LOO-CV Mean Absolute Error by Day",
      x = "Chronological age (days)",
      y = "MAE (days)"
    ) +
    pub_theme
  
  save_fig(p_mae, "Batch1_MAE_by_day", fig_dir, w = 5.8, h = 4.8)
  
  p_rmse <- ggplot(
    mae_by_day,
    aes(
      x = factor(day, levels = day_levels, labels = day_labels),
      y = RMSE,
      fill = factor(day, levels = day_levels)
    )
  ) +
    geom_col(width = 0.65, alpha = 0.9, colour = "grey20") +
    geom_hline(
      yintercept = overall_rmse, linetype = "dashed",
      color = "firebrick", linewidth = LW
    ) +
    annotate(
      "text", x = Inf, y = overall_rmse,
      label = sprintf("Overall RMSE = %.2f", overall_rmse),
      hjust = 1.05, vjust = -0.5,
      color = "firebrick"
    ) +
    scale_fill_manual(values = col_day, guide = "none") +
    labs(
      title = "LOO-CV Root Mean Squared Error by Day",
      x = "Chronological age (days)",
      y = "RMSE (days)"
    ) +
    pub_theme
  
  save_fig(p_rmse, "Batch1_RMSE_by_day", fig_dir, w = 5.8, h = 4.8)
  
  df_err <- data.frame(
    day = result$actual_day_f,
    error = result$err
  )
  
  p_err <- ggplot(
    df_err, aes(x = day, y = error, fill = day, colour = day)
  ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      colour = "grey50", linewidth = LW
    ) +
    geom_boxplot(
      width = 0.45, outlier.shape = NA,
      alpha = 0.35, linewidth = LW
    ) +
    geom_beeswarm(
      cex = 2.2, size = PT_SZ, alpha = 0.75, stroke = 0.2
    ) +
    scale_x_discrete(labels = day_labels) +
    scale_fill_manual(values = col_day, guide = "none") +
    scale_colour_manual(values = col_day, guide = "none") +
    labs(
      title = "Clock prediction error by age",
      subtitle = "LOO-CV | dashed line = perfect prediction",
      x = "Chronological age (days)",
      y = "Prediction error (days)"
    ) +
    pub_theme
  
  save_fig(
    p_err, "Batch1_prediction_error_by_day",
    fig_dir, w = 5.5, h = 4.5
  )
  
  list(
    predictions = result,
    rho = rho,
    r2 = r2,
    mae = overall_mae,
    rmse = overall_rmse,
    mae_by_day = mae_by_day,
    loo = loo
  )
}

## ====
## SECTION G: TRAIT IMPORTANCE
## ====

feature_channel <- c(
  red_n = "RED",
  red_size = "RED",
  red_elong = "RED",
  green_n = "GREEN",
  green_size = "GREEN",
  green_elong = "GREEN",
  green_ring = "GREEN",
  red_ring = "RED",
  red_radial_cv = "RED",
  green_app_net = "BOTH",
  red_app_net_norm = "BOTH",
  size_ratio = "BOTH"
)

feat_display <- c(
  red_n = "PGL-1 granule count",
  red_size = "PGL-1 granule size (px)",
  red_elong = "PGL-1 eccentricity",
  green_n = "MUT-16 foci count",
  green_size = "MUT-16 focus area (px)",
  green_elong = "MUT-16 eccentricity",
  green_ring = "MUT-16 NN-spacing CV",
  red_ring = "PGL-1 NN-spacing CV",
  red_radial_cv = "PGL-1 radial dispersion CV",
  green_app_net = paste("MUT-16", SYM_ARROW, "PGL-1 net association (%)"),
  red_app_net_norm = paste("PGL-1", SYM_ARROW, "MUT-16 net association (norm.)"),
  size_ratio = "PGL-1 / MUT-16 size ratio"
)

chan_pal <- c(
  RED = "#C0392B",
  GREEN = "#27AE60",
  BOTH = "#8E44AD"
)

trait_importance <- function(per_worm, out_dir, fig_dir,
                             final_model = NULL, loo = NULL) {
  X <- clock_matrix(per_worm)
  y <- as.numeric(per_worm$day)
  
  if (is.null(final_model)) {
    set.seed(CLOCK_SEED)
    final_model <- fit_clock_missing(X, y)
  }
  
  if (is.null(loo)) {
    loo <- clock_loo(per_worm, seed = CLOCK_SEED)
  }
  
  ## Standardized weights are beta * training SD.
  coef_df <- data.frame(
    trait = feature_cols,
    coefficient = unname(final_model$coefficient[feature_cols]),
    std_coefficient = unname(
      final_model$coefficient[feature_cols] *
        final_model$prep$sds[feature_cols]
    ),
    stringsAsFactors = FALSE
  )
  
  selection_pct <- 100 * colMeans(loo$coefficients != 0)
  coef_df$loocv_selection_pct <- round(
    unname(selection_pct[coef_df$trait]), 1
  )
  coef_df$selected <- coef_df$coefficient != 0
  coef_df <- coef_df[
    order(-abs(coef_df$std_coefficient)),
    ,
    drop = FALSE
  ]
  
  attr(coef_df, "clock_model") <- final_model
  
  write.csv(
    coef_df,
    file.path(out_dir, "Batch1_trait_importance.csv"),
    row.names = FALSE
  )
  
  cat("\nTrait importance: raw and standardized weights from one model.\n")
  print(coef_df, row.names = FALSE)
  
  abs_w <- abs(coef_df$std_coefficient)
  pct_w <- if (sum(abs_w) > 0) {
    100 * abs_w / sum(abs_w)
  } else {
    rep(0, length(abs_w))
  }
  
  df_coef <- data.frame(
    feature = coef_df$trait,
    coef = coef_df$std_coefficient,
    abs_w = abs_w,
    pct = pct_w,
    channel = unname(feature_channel[coef_df$trait]),
    active = coef_df$selected,
    stringsAsFactors = FALSE
  )
  
  df_coef$display <- unname(feat_display[df_coef$feature])
  df_coef <- df_coef[order(df_coef$coef), , drop = FALSE]
  df_coef$display <- factor(
    df_coef$display, levels = df_coef$display
  )
  
  p10 <- ggplot(
    df_coef,
    aes(x = coef, y = display, fill = channel, alpha = active)
  ) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = LW) +
    geom_col(width = 0.65, colour = "grey20", linewidth = LW) +
    geom_text(
      aes(
        label = ifelse(
          active,
          sprintf("%+.3f (%.1f%%)", coef, pct),
          "zeroed"
        ),
        x = ifelse(coef >= 0, coef + 0.01, coef - 0.01),
        hjust = ifelse(coef >= 0, 0, 1)
      ),
      colour = "grey20"
    ) +
    scale_fill_manual(values = chan_pal, name = "Channel") +
    scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.35),
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.45, 0.40))) +
    labs(
      title = "Elastic-net feature weights — GG-Clock",
      subtitle = "Standardized coefficients; percentages = share of absolute weight",
      x = "Standardized coefficient",
      y = NULL
    ) +
    pub_theme +
    theme(legend.position = c(0.88, 0.15))
  
  save_fig(
    p10, "Batch1_feature_weights",
    fig_dir, w = 8, h = 5.5
  )
  
  plot_df <- coef_df[coef_df$selected, , drop = FALSE]
  
  if (nrow(plot_df) > 0L) {
    plot_df$display <- unname(feat_display[plot_df$trait])
    
    p <- ggplot(
      plot_df,
      aes(x = reorder(display, abs(coefficient)), y = coefficient)
    ) +
      geom_col(
        aes(fill = coefficient > 0),
        width = 0.7, colour = "grey20", linewidth = LW
      ) +
      coord_flip() +
      scale_fill_manual(
        values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
        guide = "none"
      ) +
      labs(
        title = "Trait coefficients — Germ-granule Aging Clock (Batch 1, D1-6)",
        subtitle = "Elastic net, alpha = 0.5, lambda.min (full-data fit)",
        x = NULL,
        y = "Coefficient (effect on predicted day)"
      ) +
      pub_theme
    
    save_fig(
      p, "Batch1_trait_importance",
      fig_dir, w = 8.2, h = 5.4
    )
    
    p_std <- ggplot(
      plot_df,
      aes(
        x = reorder(display, abs(std_coefficient)),
        y = std_coefficient
      )
    ) +
      geom_col(
        aes(fill = std_coefficient > 0),
        width = 0.7, colour = "grey20", linewidth = LW
      ) +
      coord_flip() +
      scale_fill_manual(
        values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
        guide = "none"
      ) +
      labs(
        title = "Standardized trait coefficients — Germ-granule Aging Clock",
        subtitle = "Same fitted model; raw coefficient multiplied by feature SD",
        x = NULL,
        y = "Standardized coefficient"
      ) +
      pub_theme
    
    save_fig(
      p_std, "Batch1_trait_importance_standardized",
      fig_dir, w = 8.2, h = 5.4
    )
  } else {
    cat("All trait coefficients are zero at the selected lambda.\n")
  }
  
  coef_df
}

## ====
## SECTION H: MULTI-SEED LOO-CV STABILITY
## ====

check_seed_stability <- function(per_worm, out_dir,
                                 seeds = CLOCK_STABILITY_SEEDS) {
  y <- as.numeric(per_worm$day)
  
  out <- data.frame(
    seed = seeds,
    rho = NA_real_,
    r2 = NA_real_,
    mae = NA_real_
  )
  
  for (s in seq_along(seeds)) {
    cat(sprintf(
      "\n[SEED STABILITY] Seed %d (%d/%d)\n",
      seeds[s], s, length(seeds)
    ))
    
    loo <- clock_loo(per_worm, seed = seeds[s])
    preds <- loo$predictions
    
    out$rho[s] <- safe_cor(y, preds, method = "spearman")
    out$r2[s] <- safe_cor(y, preds)^2
    out$mae[s] <- mean(abs(preds - y))
  }
  
  cat("\n[SEED STABILITY] Summary:\n")
  
  for (metric in c("rho", "r2", "mae")) {
    v <- out[[metric]]
    v <- v[is.finite(v)]
    
    if (!length(v)) {
      cat(metric, ": undefined for all seeds\n")
    } else {
      cat(sprintf(
        "%s: mean=%.3f, SD=%.3f, range=%.3f-%.3f\n",
        metric, mean(v), finite_sd(v), min(v), max(v)
      ))
    }
  }
  
  write.csv(
    out,
    file.path(out_dir, "Batch1_seed_stability.csv"),
    row.names = FALSE
  )
  
  out
}

## ====
## CLOCK ATTRIBUTION
## ====

clock_attribution_IF_vs_AL <- function(per_worm_AL, per_worm_IF,
                                       final_model, day = 4,
                                       out_dir, fig_dir, tag = "",
                                       n_boot = 5000, seed = CLOCK_SEED,
                                       focus_traits = c("green_n", "red_radial_cv"),
                                       conf = 0.95,
                                       min_per_group = 3L) {
  if (is.null(final_model)) stop("final_model is required.")
  if (is.null(per_worm_IF) || nrow(per_worm_IF) == 0L) {
    stop("No IF per-worm data available.")
  }
  if (!all(focus_traits %in% feature_cols)) {
    stop("focus_traits must be clock feature names.")
  }
  
  al4 <- per_worm_AL[per_worm_AL$day == day, , drop = FALSE]
  fi4 <- per_worm_IF[per_worm_IF$day == day, , drop = FALSE]
  nA <- nrow(al4)
  nF <- nrow(fi4)
  
  if (nA < min_per_group || nF < min_per_group) {
    stop(
      "Need at least ", min_per_group,
      " worms per group at day ", day, "."
    )
  }
  
  cat(sprintf(
    "\n[ATTRIBUTION] Day %d: n_AL=%d; n_IF=%d; bootstrap n=%d\n",
    day, nA, nF, n_boot
  ))
  
  ## Fixed weights and AL-derived preprocessing.
  beta <- final_model$coefficient[feature_cols]
  sds <- final_model$prep$sds[feature_cols]
  
  X_al <- clock_apply(clock_matrix(al4), final_model$prep)[
    , feature_cols, drop = FALSE
  ]
  X_if <- clock_apply(clock_matrix(fi4), final_model$prep)[
    , feature_cols, drop = FALSE
  ]
  
  contrib_fun <- function(Xa, Xf) {
    beta * (colMeans(Xf) - colMeans(Xa))
  }
  
  obs <- contrib_fun(X_al, X_if)
  total <- sum(obs)
  
  ## 
  pred_diff <- mean(predict_clock_missing(final_model, clock_matrix(fi4))) -
    mean(predict_clock_missing(final_model, clock_matrix(al4)))
  
  stopifnot(abs(total - pred_diff) < 1e-8)
  
  ## ---- Worm-level bootstrap ----
  set.seed(seed)
  
  boot <- matrix(
    NA_real_, nrow = n_boot, ncol = length(feature_cols),
    dimnames = list(NULL, feature_cols)
  )
  
  for (b in seq_len(n_boot)) {
    ia <- sample.int(nA, nA, replace = TRUE)
    ib <- sample.int(nF, nF, replace = TRUE)
    boot[b, ] <- contrib_fun(
      X_al[ia, , drop = FALSE],
      X_if[ib, , drop = FALSE]
    )
  }
  
  boot_total <- rowSums(boot)
  a <- (1 - conf) / 2
  
  ## ---- Per-trait table ----
  res <- data.frame(
    trait = feature_cols,
    display = unname(feat_display[feature_cols]),
    channel = unname(feature_channel[feature_cols]),
    in_clock = unname(beta != 0),
    clock_beta = unname(beta),
    training_SD = unname(sds),
    mean_AL = unname(colMeans(X_al)),
    mean_IF = unname(colMeans(X_if)),
    shift_raw = unname(colMeans(X_if) - colMeans(X_al)),
    shift_training_SD = unname(
      ifelse(sds > 0, (colMeans(X_if) - colMeans(X_al)) / sds, NA_real_)
    ),
    contribution_days = unname(obs),
    ci_low = apply(boot, 2, quantile, probs = a, names = FALSE),
    ci_high = apply(boot, 2, quantile, probs = 1 - a, names = FALSE),
    share_of_total_pct = if (abs(total) > 1e-12) {
      unname(100 * obs / total)
    } else {
      NA_real_
    },
    ## One-sided bootstrap p for IF (contribution < 0).
    p_boot_younger = vapply(
      feature_cols,
      function(f) (1 + sum(boot[, f] >= 0)) / (n_boot + 1),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  res$reliable <- res$in_clock &
    (res$ci_low > 0 | res$ci_high < 0)
  res$direction <- ifelse(
    !res$in_clock, "not in clock",
    ifelse(res$contribution_days < 0, "younger", "older")
  )
  
  res <- res[order(res$contribution_days), , drop = FALSE]
  
  ## ---- Total and focus-trait share ----
  total_ci <- quantile(boot_total, c(a, 1 - a), names = FALSE)
  
  focus_obs <- sum(obs[focus_traits])
  focus_boot <- rowSums(boot[, focus_traits, drop = FALSE])
  focus_ci <- quantile(focus_boot, c(a, 1 - a), names = FALSE)
  
  ## Share = focus / total, only meaningful where the bootstrap
  ## total keeps the observed sign (otherwise ratio explodes).
  same_sign <- sign(boot_total) == sign(total) & abs(boot_total) > 1e-12
  focus_share_obs <- if (abs(total) > 1e-12) 100 * focus_obs / total else NA_real_
  focus_share_boot <- 100 * focus_boot[same_sign] / boot_total[same_sign]
  focus_share_ci <- if (length(focus_share_boot) >= 100L) {
    quantile(focus_share_boot, c(a, 1 - a), names = FALSE)
  } else {
    c(NA_real_, NA_real_)
  }
  
  summary_tbl <- data.frame(
    day = day,
    n_AL = nA,
    n_IF = nF,
    n_boot = n_boot,
    total_shift_days = total,
    total_ci_low = total_ci[1],
    total_ci_high = total_ci[2],
    p_boot_total_younger = (1 + sum(boot_total >= 0)) / (n_boot + 1),
    focus_traits = paste(focus_traits, collapse = "+"),
    focus_contribution_days = focus_obs,
    focus_ci_low = focus_ci[1],
    focus_ci_high = focus_ci[2],
    focus_share_pct = focus_share_obs,
    focus_share_ci_low = focus_share_ci[1],
    focus_share_ci_high = focus_share_ci[2],
    pct_boot_same_sign_total = 100 * mean(same_sign),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    res,
    file.path(out_dir, paste0(
      "Batch1_IF_vs_AL_clock_attribution_D", day, tag, ".csv"
    )),
    row.names = FALSE
  )
  write.csv(
    summary_tbl,
    file.path(out_dir, paste0(
      "Batch1_IF_vs_AL_clock_attribution_summary_D", day, tag, ".csv"
    )),
    row.names = FALSE
  )
  write.csv(
    data.frame(boot_id = seq_len(n_boot), boot, total = boot_total),
    file.path(out_dir, paste0(
      "Batch1_IF_vs_AL_clock_attribution_boot_D", day, tag, ".csv"
    )),
    row.names = FALSE
  )
  
  ## ---- Console report ----
  cat(sprintf(
    "\n[ATTRIBUTION] Total IF - AL clock shift = %+.3f d (%.0f%% CI %+.3f to %+.3f)\n",
    total, 100 * conf, total_ci[1], total_ci[2]
  ))
  cat(sprintf(
    "[ATTRIBUTION] Bootstrap one-sided p (IF younger) = %.4f\n",
    summary_tbl$p_boot_total_younger
  ))
  cat("\n--- Per-trait contributions (days; negative = younger) ---\n")
  print(
    res[, c("display", "contribution_days", "ci_low", "ci_high",
            "share_of_total_pct", "p_boot_younger", "reliable",
            "direction")],
    row.names = FALSE, digits = 3
  )
  cat(sprintf(
    paste0(
      "\n[ATTRIBUTION] %s: %+.3f d (%.0f%% CI %+.3f to %+.3f) = ",
      "%.1f%% of the total shift (bootstrap %.0f%% CI %.1f to %.1f%%).\n",
      "[ATTRIBUTION] Share CI uses the %.1f%% of resamples whose total ",
      "has the observed sign.\n"
    ),
    paste(unname(feat_display[focus_traits]), collapse = " + "),
    focus_obs, 100 * conf, focus_ci[1], focus_ci[2],
    focus_share_obs, 100 * conf,
    focus_share_ci[1], focus_share_ci[2],
    100 * mean(same_sign)
  ))
  cat(paste0(
    "[ATTRIBUTION] Clock weights fixed; clock-training uncertainty is ",
    "not propagated. AL day predictions are in-sample.\n"
  ))
  
  ## ---- Figure 1: contributions with bootstrap CI ----
  pl <- res
  pl$display <- factor(pl$display, levels = pl$display)
  pl$focus <- pl$trait %in% focus_traits
  
  pl <- res[res$in_clock, , drop = FALSE]          # drop zero-weight traits, as in Batch 2
  pl$display <- factor(pl$display, levels = pl$display[order(
    pl$contribution_days, decreasing = TRUE)])
  pl$channel <- factor(
    pl$channel, c("RED", "GREEN", "BOTH"),
    c("PGL-1 (red)", "MUT-16 (green)", "Both channels")
  )
  pl$lab <- sprintf("%.0f%%", pl$share_of_total_pct)
  n_zero <- sum(!res$in_clock)
  
  p_c <- ggplot(pl, aes(x = contribution_days, y = display, fill = channel)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = LW) +
    geom_col(width = 0.65, colour = "grey20", linewidth = LW) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y", width = 0.25, linewidth = LW, colour = "grey20"
    ) +
    geom_text(aes(x = pmax(ci_high, 0), label = lab), hjust = -0.3) +
    scale_fill_manual(
      values = c("PGL-1 (red)" = "#C0392B",
                 "MUT-16 (green)" = "#3FA66B",
                 "Both channels" = "#7D3C98"),
      name = NULL, drop = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = "Trait contributions to the younger clock age in IF worms",
      subtitle = sprintf(
        "Batch 1 (discovery), day %d | total shift %+.2f d (%.0f%% CI %+.2f to %+.2f)",
        day, total, 100 * conf, total_ci[1], total_ci[2]
      ),
      x = "Contribution to clock age shift (days; negative = younger)",
      y = NULL,
      caption = sprintf(
        "Whiskers: %.0f%% bootstrap CI over worms (n = %s). Clock weights fixed.%s",
        100 * conf, format(n_boot, big.mark = ","),
        if (n_zero > 0L) sprintf(" %d zero-weight trait(s) omitted.", n_zero) else ""
      )
    ) +
    pub_theme +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, colour = "grey30")
    )
  
  save_fig(
    p_c,
    paste0("Batch1_clock_attribution_D", day, tag),
    fig_dir, w = 8.5, h = 6
  )
  
  ## ---- Figure 2: bootstrap distribution of focus-trait share ----
  if (length(focus_share_boot) >= 100L) {
    p_s <- ggplot(data.frame(s = focus_share_boot), aes(x = s)) +
      geom_histogram(bins = 50, fill = "grey75", colour = "grey30") +
      geom_vline(
        xintercept = focus_share_obs,
        colour = "#D7191C", linewidth = LW
      ) +
      labs(
        title = "Share of clock shift from focus traits",
        subtitle = sprintf(
          "%s | observed = %.1f%% (%.0f%% CI %.1f to %.1f%%)",
          paste(unname(feat_display[focus_traits]), collapse = " + "),
          focus_share_obs, 100 * conf,
          focus_share_ci[1], focus_share_ci[2]
        ),
        x = "Share of total shift (%)", y = "Bootstrap count"
      ) +
      pub_theme
    
    save_fig(
      p_s,
      paste0("Batch1_clock_attribution_focus_share_D", day, tag),
      fig_dir, w = 6.5, h = 4.8
    )
  }
  
  invisible(list(
    table = res,
    summary = summary_tbl,
    boot = boot,
    boot_total = boot_total
  ))
}

## ====
## SECTION I: IF EXTRACTION
## ====

get_AL_batch_ref <- function(b1) {
  ref <- b1$batch_ref
  validate_noise_reference(ref, "Stored AL batch noise reference")
  
  cat(sprintf(
    "[IF] Using stored AL batch reference local_bg_mad = %.6g\n",
    ref
  ))
  
  ref
}

process_batch1_IF <- function(root_dir, out_dir, overlay_dir,
                              batch_ref = NULL, tag = "") {
  day_dirs <- list.dirs(root_dir, recursive = FALSE)
  day_dirs <- day_dirs[grepl("^B1_D[0-9]+$", basename(day_dirs))]
  
  if (!length(day_dirs)) {
    stop("No B1_D* folders found under: ", root_dir)
  }
  
  all_pairs <- list()
  
  for (dd in day_dirs) {
    worm_folders <- list.dirs(dd, recursive = FALSE)
    worm_folders <- worm_folders[
      grepl("^IF_B1_D[0-9]+_W[0-9]+$", basename(worm_folders))
    ]
    
    for (wf in worm_folders) {
      pr <- collect_worm_pairs_IF(wf)
      
      if (!is.null(pr)) {
        all_pairs[[length(all_pairs) + 1L]] <- pr
      }
    }
  }
  
  pair_tbl <- do.call(rbind, all_pairs)
  
  if (is.null(pair_tbl) || nrow(pair_tbl) == 0L) {
    cat("\n[IF] No matched IF nucleus pairs found; skipping.\n")
    return(NULL)
  }
  
  if (anyDuplicated(pair_tbl$base_id)) {
    stop("Duplicate IF base_id values detected.")
  }
  
  cat(sprintf(
    "\n[IF] Matched %d nucleus pairs across %d IF worms.\n",
    nrow(pair_tbl), length(unique(pair_tbl$worm_id))
  ))
  
  ## Separate sensitivity overlays so an optional rerun does
  ## not overwrite the main IF segmentation overlays.
  if_overlay_dir <- file.path(overlay_dir, paste0("IF", tag))
  dir.create(if_overlay_dir, recursive = TRUE, showWarnings = FALSE)
  
  noise_info <- compute_image_noise_scale_IF(
    pair_tbl, batch_ref = batch_ref
  )
  
  write.csv(
    noise_info$scale_tbl,
    file.path(
      out_dir, paste0("Batch1_IF_image_noise_scale", tag, ".csv")
    ),
    row.names = FALSE
  )
  
  pair_tbl <- merge(
    pair_tbl,
    noise_info$scale_tbl[, c("base_id", "noise_ratio")],
    by = "base_id", all.x = TRUE
  )
  
  pair_tbl$noise_ratio[!is.finite(pair_tbl$noise_ratio)] <- 1.0
  pair_tbl$adaptive_z_min <- LOCAL_CONTRAST_Z_MIN * pair_tbl$noise_ratio
  
  cat(sprintf(
    "[ADAPTIVE-IF] Effective z-min range: %.3f to %.3f (median %.3f)\n",
    min(pair_tbl$adaptive_z_min),
    max(pair_tbl$adaptive_z_min),
    median(pair_tbl$adaptive_z_min)
  ))
  
  feat_list <- lapply(seq_len(nrow(pair_tbl)), function(i) {
    if (i %% 20 == 0 || i == 1) {
      cat(sprintf(
        "  [IF %d/%d] %s\n",
        i, nrow(pair_tbl), pair_tbl$base_id[i]
      ))
    }
    
    meta <- as.list(pair_tbl[
      i, c("condition", "day", "worm_id", "nucleus", "base_id")
    ])
    
    tryCatch(
      extract_features_real(
        pair_tbl$red_path[i],
        pair_tbl$green_path[i],
        meta,
        overlay_dir = if_overlay_dir,
        contrast_z_min = pair_tbl$adaptive_z_min[i]
      ),
      error = function(e) {
        warning("SKIP: ", pair_tbl$base_id[i], " - ", e$message)
        NULL
      }
    )
  })
  
  per_image_IF <- do.call(
    rbind, Filter(Negate(is.null), feat_list)
  )
  
  if (is.null(per_image_IF) || nrow(per_image_IF) == 0L) {
    cat("[IF] No successfully processed IF images; skipping.\n")
    return(NULL)
  }
  
  tables <- write_feature_tables(
    per_image_IF,
    out_dir,
    prefix = "Batch1_IF",
    tag = tag
  )
  
  tables$per_worm
}

## ====
## IF vs AL COMPARISON
## ====

project_IF_vs_AL <- function(per_worm_AL, per_worm_IF,
                             al_loocv, out_dir, fig_dir,
                             tag = "", final_model = NULL) {
  if (is.null(per_worm_IF) || nrow(per_worm_IF) == 0L) {
    cat("\n[IF] Skipping comparison: no IF worm data.\n")
    return(NULL)
  }
  
  X_al <- clock_matrix(per_worm_AL)
  X_if <- clock_matrix(per_worm_IF)
  y_al <- as.numeric(per_worm_AL$day)
  
  if (is.null(final_model)) {
    set.seed(CLOCK_SEED)
    final_model <- fit_clock_missing(X_al, y_al)
  }
  
  cv_al <- final_model
  if_pred <- predict_clock_missing(cv_al, X_if)
  
  G_AL <- "AL (LOO-CV)"
  G_IF <- "IF (projected)"
  
  if_tbl <- data.frame(
    worm_id = per_worm_IF$worm_id,
    group = G_IF,
    day = per_worm_IF$day,
    predicted = if_pred,
    stringsAsFactors = FALSE
  )
  
  al_tbl <- data.frame(
    worm_id = al_loocv$worm_id,
    group = G_AL,
    day = al_loocv$actual_day,
    predicted = al_loocv$predicted_day,
    stringsAsFactors = FALSE
  )
  
  if (!all(is.finite(al_tbl$predicted))) {
    stop("AL LOO predictions contain non-finite values.")
  }
  if (!all(is.finite(if_tbl$day))) {
    stop("IF chronological ages contain non-finite values.")
  }
  
  if_days <- sort(unique(if_tbl$day))
  
  if (!all(if_days %in% al_tbl$day)) {
    stop(
      "Cannot calculate same-day AL correction for IF day(s): ",
      paste(setdiff(if_days, al_tbl$day), collapse = ", ")
    )
  }
  
  al_tbl <- al_tbl[al_tbl$day %in% if_days, , drop = FALSE]
  
  comb <- rbind(al_tbl, if_tbl)
  comb$error <- comb$predicted - comb$day
  comb$group <- factor(comb$group, levels = c(G_AL, G_IF))
  
  al_bias <- tapply(
    comb$error[comb$group == G_AL],
    comb$day[comb$group == G_AL],
    mean
  )
  
  ## 
  raw_accel <- if_pred - per_worm_IF$day
  same_day_bias <- as.numeric(
    al_bias[as.character(per_worm_IF$day)]
  )
  
  if_out <- data.frame(
    worm_id = per_worm_IF$worm_id,
    n_nuclei = per_worm_IF$n_nuclei,
    true_day = per_worm_IF$day,
    predicted_age = if_pred,
    age_accel_raw = raw_accel,
    al_bias_same_day = same_day_bias,
    age_accel_corrected = raw_accel - same_day_bias,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    if_out,
    file.path(
      out_dir, paste0("Batch1_IF_projected_onto_AL_clock", tag, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    comb,
    file.path(out_dir, paste0("Batch1_IF_vs_AL_combined", tag, ".csv")),
    row.names = FALSE
  )
  
  cat("\n==== IF PROJECTED ONTO AL CLOCK", tag, "====\n")
  print(if_out, row.names = FALSE, digits = 4)
  
  per_day <- do.call(rbind, lapply(if_days, function(d) {
    a <- comb$error[comb$group == G_AL & comb$day == d]
    f <- comb$error[comb$group == G_IF & comb$day == d]
    
    p <- if (length(a) >= 2L && length(f) >= 2L) {
      suppressWarnings(
        wilcox.test(
          f, a, alternative = "less", exact = FALSE
        )$p.value
      )
    } else {
      NA_real_
    }
    
    if (!is.finite(p)) p <- NA_real_
    
    data.frame(
      day = d,
      n_AL = length(a),
      n_IF = length(f),
      mean_err_AL = mean(a),
      mean_err_IF = mean(f),
      IF_minus_AL = mean(f) - mean(a),
      p_one_sided = p
    )
  }))
  
  write.csv(
    per_day,
    file.path(out_dir, paste0("Batch1_IF_vs_AL_by_day", tag, ".csv")),
    row.names = FALSE
  )
  
  cat("\n--- Same-day prediction errors ---\n")
  print(per_day, row.names = FALSE, digits = 3)
  
  fit <- if (length(unique(comb$day)) > 1L) {
    lm(error ~ factor(day) + group, data = comb)
  } else {
    lm(error ~ group, data = comb)
  }
  
  cf <- summary(fit)$coefficients
  rn <- grep("^group", rownames(cf), value = TRUE)
  
  est <- NA_real_
  ci <- c(NA_real_, NA_real_)
  p_two <- NA_real_
  p_less <- NA_real_
  
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
  
  cat(sprintf(
    paste0(
      "\n[IF vs AL%s] Adjusted IF effect: %+.3f d ",
      "(95%% CI %+.3f to %+.3f)\n",
      "Two-sided p=%.4f; one-sided IF-younger p=%.4f\n"
    ),
    tag, est, ci[1], ci[2], p_two, p_less
  ))
  
  plab <- function(p) {
    if (!is.finite(p)) {
      "p unavailable"
    } else if (p < 0.0001) {
      "p < 0.0001"
    } else {
      sprintf("p = %.3f", p)
    }
  }
  
  day_levels <- sort(unique(c(y_al, if_days)))
  day_labels <- paste0("D", day_levels)
  cols <- setNames(c("#4D4D4D", "#D7191C"), c(G_AL, G_IF))
  
  set.seed(1)
  off <- ifelse(comb$group == G_IF, 0.16, -0.16)
  comb$x_plot <- comb$day + off + runif(nrow(comb), -0.05, 0.05)
  
  sm <- do.call(rbind, lapply(
    split(comb, list(comb$group, comb$day), drop = TRUE),
    function(df) {
      data.frame(
        group = df$group[1],
        day = df$day[1],
        xc = df$day[1] +
          ifelse(df$group[1] == G_IF, 0.16, -0.16),
        m = mean(df$predicted),
        s = finite_sd(df$predicted)
      )
    }
  ))
  
  ## Missing SD stays NA
  p_sc <- ggplot(
    comb,
    aes(x = x_plot, y = predicted, colour = group, shape = group)
  ) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = "grey50", linewidth = LW
    ) +
    geom_point(size = PT_SZ, alpha = 0.7) +
    geom_crossbar(
      data = sm,
      aes(x = xc, y = m, ymin = m - s, ymax = m + s, colour = group),
      inherit.aes = FALSE,
      width = 0.26, linewidth = LW,
      fill = NA, na.rm = TRUE
    ) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_x_continuous(breaks = day_levels, labels = day_labels) +
    scale_y_continuous(breaks = day_levels, labels = day_labels) +
    coord_cartesian(
      xlim = range(day_levels) + c(-0.6, 0.6),
      ylim = range(c(day_levels, comb$predicted)) + c(-0.3, 0.3)
    ) +
    annotate(
      "text",
      x = min(day_levels) - 0.5,
      y = max(day_levels) + 0.3,
      hjust = 0, vjust = 1,
      label = sprintf(
        "IF vs AL (same day): %+.2f d\n%s (one-sided)",
        est, plab(p_less)
      )
    ) +
    labs(
      title = "GG-Clock: AL vs IF projection",
      subtitle = "Bars = mean +/- SD | dashed line = perfect prediction",
      x = "Chronological age (days)",
      y = "Predicted age (days)"
    ) +
    pub_theme +
    theme(legend.position = "bottom")
  
  save_fig(
    p_sc, paste0("Batch1_IF_vs_AL_scatter", tag),
    fig_dir, w = 6.5, h = 6
  )
  
  comb$day_f <- factor(comb$day, levels = day_levels)
  
  p_err <- ggplot(
    comb, aes(x = day_f, y = error, fill = group, colour = group)
  ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      colour = "grey50", linewidth = LW
    ) +
    geom_boxplot(
      width = 0.6, outlier.shape = NA,
      alpha = 0.3, linewidth = LW,
      position = position_dodge(width = 0.75)
    ) +
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.12, dodge.width = 0.75
      ),
      size = PT_SZ, alpha = 0.75
    ) +
    scale_x_discrete(
      labels = setNames(day_labels, as.character(day_levels))
    ) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(
      title = "Prediction error by age: AL vs IF",
      subtitle = sprintf(
        "Same-day adjusted IF effect = %+.2f d (%s, one-sided)",
        est, plab(p_less)
      ),
      x = "Chronological age (days)",
      y = "Prediction error (days)"
    ) +
    pub_theme +
    theme(legend.position = "bottom")
  
  save_fig(
    p_err, paste0("Batch1_IF_vs_AL_error_by_day", tag),
    fig_dir, w = 6.5, h = 5
  )
  
  list(
    if_table = if_out,
    comparison = comb,
    per_day = per_day,
    lm_fit = fit,
    effect = est,
    ci = ci,
    p_two = p_two,
    p_less = p_less,
    final_cv = cv_al
  )
}

## ====
## RUN — AL BATCH 1 (DAYS 1-6)
## ====

cat("\n", strrep("*", 70), "\n")
cat("  GERM-GRANULE AGING CLOCK — BATCH 1 (DAYS 1-6)\n")
cat("  Adaptive noise scaling; edge-margin gate removed\n")
cat("  Conditional feature masking; nested training-only imputation\n")
cat(
  "  Green thresholds: DoG =", GREEN_THR_MULT,
  "| global =", GREEN_GLOBAL_MULT,
  "| max eccentricity =", GREEN_MAX_ECCENTRICITY,
  "| baseline local z-min =", LOCAL_CONTRAST_Z_MIN, "\n"
)
cat("  Output folder (new):", OUT_ROOT, "\n")
cat(strrep("*", 70), "\n\n")

b1 <- process_batch1(B1_ROOT)

required_cols <- c(
  "pass_local_contrast", "global_floor", "floor_gap",
  "peakiness", "solidity", "adaptive_z_min"
)

missing_cols <- setdiff(required_cols, names(b1$candidate_qc))

if (length(missing_cols)) {
  stop(
    "Candidate QC table is missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}
if ("pass_edge_margin" %in% names(b1$candidate_qc)) {
  stop("Unexpected edge-margin gate in candidate QC.")
}

clock <- run_clock_LOOCV(
  b1$per_worm, b1$out_dir, b1$fig_dir
)

## Fit the final AL model ONCE.
## Reuse it for importance, IF projection and attribution.
set.seed(CLOCK_SEED)

final_AL_clock <- fit_clock_missing(
  clock_matrix(b1$per_worm),
  as.numeric(b1$per_worm$day)
)

if (final_AL_clock$boundary) {
  warning(
    "Final AL model selected a lambda-grid boundary. ",
    "Inspect Batch1_final_clock_tuning.csv."
  )
}

write.csv(
  data.frame(
    lambda = final_AL_clock$lambda,
    inner_CV_MSE = final_AL_clock$cvm,
    selected = final_AL_clock$lambda == final_AL_clock$lambda.min
  ),
  file.path(b1$out_dir, "Batch1_final_clock_tuning.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    trait = feature_cols,
    imputation_median = unname(final_AL_clock$prep$med[feature_cols]),
    training_SD = unname(final_AL_clock$prep$sds[feature_cols]),
    all_missing = unname(final_AL_clock$prep$all_missing[feature_cols]),
    disabled = unname(final_AL_clock$prep$disabled[feature_cols])
  ),
  file.path(b1$out_dir, "Batch1_final_clock_preprocessing.csv"),
  row.names = FALSE
)

saveRDS(
  final_AL_clock,
  file.path(b1$out_dir, "Batch1_final_AL_clock.rds")
)

importance <- trait_importance(
  b1$per_worm,
  b1$out_dir,
  b1$fig_dir,
  final_model = final_AL_clock,
  loo = clock$loo
)

seed_check <- check_seed_stability(
  b1$per_worm,
  b1$out_dir,
  seeds = CLOCK_STABILITY_SEEDS
)

## ====
## RUN — IF WITH SHARED AL NOISE REFERENCE
## ====

AL_REF <- get_AL_batch_ref(b1)

per_worm_IF <- process_batch1_IF(
  B1_ROOT,
  b1$out_dir,
  b1$overlay_dir,
  batch_ref = AL_REF,
  tag = "_ALref"
)

if_projection <- project_IF_vs_AL(
  b1$per_worm,
  per_worm_IF,
  clock$predictions,
  b1$out_dir,
  b1$fig_dir,
  tag = "_ALref",
  final_model = final_AL_clock
)

if (!is.null(if_projection)) {
  cat("\nIF results written to:\n")
  
  for (f in c(
    "Batch1_IF_per_image_zerosIncluded_ALref.csv",
    "Batch1_IF_per_image_ALref.csv",
    "Batch1_IF_per_worm_zerosIncluded_ALref.csv",
    "Batch1_IF_per_worm_ALref.csv",
    "Batch1_IF_projected_onto_AL_clock_ALref.csv",
    "Batch1_IF_vs_AL_combined_ALref.csv",
    "Batch1_IF_vs_AL_by_day_ALref.csv"
  )) {
    cat("  ", file.path(b1$out_dir, f), "\n")
  }
} else {
  cat("\nIF comparison skipped: no usable IF data.\n")
}

## ====
## IF vs AL: CLOCK ATTRIBUTION 
## IF comparison day stays at day 4.
## ====

if (!is.null(if_projection)) {
  attribution_D4 <- clock_attribution_IF_vs_AL(
    per_worm_AL = b1$per_worm,
    per_worm_IF = per_worm_IF,
    final_model = final_AL_clock,
    day = 4,
    out_dir = b1$out_dir,
    fig_dir = b1$fig_dir,
    tag = "_ALref",
    n_boot = 5000,
    focus_traits = c("green_n", "red_radial_cv")
  )
}

## ====
## IF-OWN-NOISE-REFERENCE SENSITIVITY
## ====

if (RUN_OLD_REF_COMPARISON && !is.null(if_projection)) {
  per_worm_IF_old <- process_batch1_IF(
    B1_ROOT,
    b1$out_dir,
    b1$overlay_dir,
    batch_ref = NULL,
    tag = "_IFownref"
  )
  
  if_projection_old <- project_IF_vs_AL(
    b1$per_worm,
    per_worm_IF_old,
    clock$predictions,
    b1$out_dir,
    b1$fig_dir,
    tag = "_IFownref",
    final_model = final_AL_clock
  )
  
  if (!is.null(if_projection_old)) {
    cat(sprintf(
      "\n[SENSITIVITY] AL reference: %+.3f d (two-sided p=%.4f)\n",
      if_projection$effect,
      if_projection$p_two
    ))
    cat(sprintf(
      "[SENSITIVITY] IF-own reference: %+.3f d (two-sided p=%.4f)\n",
      if_projection_old$effect,
      if_projection_old$p_two
    ))
  }
}

cat("\nDone. Outputs in:\n", b1$out_dir, "\n")
cat("Figures in:\n", b1$fig_dir, "\n")
cat("QC overlays in:\n", b1$overlay_dir, "\n")