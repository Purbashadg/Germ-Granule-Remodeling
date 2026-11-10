# ============================================================
# GGProfiler - v1.0
# ============================================================

# ---- Dependency check and optional one-click installation -------------------
CRAN_PKGS <- c("shiny", "ggplot2", "DT")   # required, from CRAN
BIOC_PKGS <- c("EBImage")                  # required, from Bioconductor
OPTIONAL_PKGS <- c("zip")                  # optional: faster, safer ZIP export

is_installed <- function(pkgs) {
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
}

install_missing_packages <- function(include_optional = TRUE) {
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos[["CRAN"]]), "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  
  cran_needed <- CRAN_PKGS[!is_installed(CRAN_PKGS)]
  if (include_optional) {
    cran_needed <- c(cran_needed, OPTIONAL_PKGS[!is_installed(OPTIONAL_PKGS)])
  }
  bioc_needed <- BIOC_PKGS[!is_installed(BIOC_PKGS)]
  
  if (length(cran_needed)) {
    message("Installing from CRAN: ", paste(cran_needed, collapse = ", "))
    install.packages(cran_needed)
  }
  
  if (length(bioc_needed)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    message("Installing from Bioconductor: ", paste(bioc_needed, collapse = ", "))
    BiocManager::install(bioc_needed, update = FALSE, ask = FALSE)
  }
  invisible(TRUE)
}

required <- c(CRAN_PKGS, BIOC_PKGS)
missing <- required[!is_installed(required)]

if (length(missing)) {
  if (interactive()) {
    message(
      "GGProfiler needs the following missing package(s):\n  ",
      paste(missing, collapse = ", ")
    )
    answer <- tolower(trimws(readline(
      "Download and install them now? [y/n]: "
    )))
    if (answer %in% c("y", "yes")) {
      tryCatch(
        install_missing_packages(include_optional = TRUE),
        error = function(e) {
          stop("Package installation failed: ", conditionMessage(e), call. = FALSE)
        }
      )
      missing <- required[!is_installed(required)]
    }
  }
  
  if (length(missing)) {
    stop(
      "Install these packages before running: ",
      paste(missing, collapse = ", "),
      "\nCRAN: install.packages(c(",
      paste0('"', intersect(missing, CRAN_PKGS), '"', collapse = ", "), "))",
      if (any(missing %in% BIOC_PKGS)) paste0(
        "\nBioconductor: BiocManager::install(c(",
        paste0('"', intersect(missing, BIOC_PKGS), '"', collapse = ", "), "))"
      ) else "",
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(shiny)
  library(EBImage)
  library(ggplot2)
  library(DT)
})

suppressPackageStartupMessages({
  library(shiny)
  library(EBImage)
  library(ggplot2)
  library(DT)
})

options(shiny.maxRequestSize = 2 * 1024^3) # 2 GiB upload limit

APP_VERSION <- "GGProfiler reliability refactor 1 + nuclear crop module"
MAX_CONDS <- 20L
MAX_CHANNELS <- 4L

`%||%` <- function(a, b) if (is.null(a)) b else a

COLORS <- c(
  Cyan = "#00BCD4",
  Magenta = "#CC79A7",
  Yellow = "#F0E442",
  Green = "#009E73",
  Red = "#D55E00",
  Blue = "#0072B2",
  Orange = "#E69F00",
  White = "#FFFFFF"
)

METHODS <- c(
  "Top-hat: diffuse / larger objects" = "tophat",
  "DoG: punctate objects" = "dog"
)

SOURCES <- c(
  "R / plane 1" = 1,
  "G / plane 2" = 2,
  "B / plane 3" = 3
)

OVERLAYS <- c(
  "Peak positions" = "peaks",
  "Nearest-neighbor links" = "nnd",
  "Radial spokes: image center" = "radial",
  "B peaks near A objects" = "green_app_net",
  "A peaks near B objects" = "red_app_net_norm"
)

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

cv_value <- function(x, minimum = 3L) {
  x <- x[is.finite(x)]
  if (length(x) < minimum || mean(x) <= 1e-12) return(NA_real_)
  sd(x) / mean(x)
}

bind_rows <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- Image operations -------------------------------------------------------

eb <- function(m) EBImage::Image(m, colormode = "Grayscale")

matrix_data <- function(x) {
  a <- EBImage::imageData(x)
  if (length(dim(a)) != 2L) stop("Expected a two-dimensional image.")
  matrix(as.numeric(a), nrow = dim(a)[1], ncol = dim(a)[2])
}

read_channel <- function(path, src) {
  im <- EBImage::readImage(path)
  a <- EBImage::imageData(im)
  d <- dim(a)
  
  if (length(d) == 2L) {
    # A genuine grayscale image has only one plane; use it regardless
    # of the RGB selector, matching the intended original behavior.
    m <- matrix(as.numeric(a), d[1], d[2])
  } else if (
    length(d) == 3L &&
    EBImage::colorMode(im) == EBImage::Color &&
    d[3] >= 3L
  ) {
    if (src %in% 1:3) {
      m <- matrix(as.numeric(a[, , src]), d[1], d[2])
    } else if (src == 4L) {
      m <- (a[, , 1] + a[, , 2] + a[, , 3]) / 3
    } else {
      stop("Invalid RGB source-plane selection.")
    }
  } else {
    stop(
      "Unsupported image dimensions/type: ",
      paste(d, collapse = " x "),
      ". Use a single 2D grayscale or RGB image, not a stack."
    )
  }
  
  if (min(dim(m)) < 7L) stop("Image must be at least 7 x 7 pixels.")
  if (any(!is.finite(m))) stop("Image contains non-finite pixels.")
  m
}

brush <- function(size, dims) {
  size <- as.integer(size)
  if (!is.finite(size) || size < 1L) stop("Invalid brush size.")
  limit <- min(dims) - 2L
  size <- min(size, limit)
  if (size %% 2L == 0L) size <- size - 1L
  EBImage::makeBrush(max(size, 1L), shape = "disc")
}

morph <- function(fn, m, size) {
  matrix_data(fn(eb(m * 1), brush(size, dim(m))))
}

blur <- function(m, sigma) {
  if (!is.finite(sigma) || sigma <= 0) stop("Sigma must be positive.")
  matrix_data(EBImage::gblur(eb(m), sigma = sigma))
}

label_mask <- function(mask) {
  out <- matrix_data(EBImage::bwlabel(eb(mask * 1)))
  storage.mode(out) <- "integer"
  out
}

relabel_keep <- function(lab, keep) {
  out <- matrix(0L, nrow(lab), ncol(lab))
  for (j in seq_along(keep)) out[lab == keep[j]] <- as.integer(j)
  out
}

candidate_labels <- function(response, multiplier) {
  pos <- response[response > 0]
  if (length(pos) < 10L) {
    return(matrix(0L, nrow(response), ncol(response)))
  }
  
  threshold <- mean(pos) + multiplier * sd(pos)
  bw <- matrix_data(EBImage::fillHull(eb((response > threshold) * 1)))
  bw <- morph(EBImage::closing, bw, 3L)
  label_mask(bw > 0)
}

dog_response <- function(m, p) {
  coarse_sigma <- min(p$sigma2, min(dim(m)) / 6)
  if (coarse_sigma <= p$sigma1) {
    stop(
      "Effective coarse sigma must exceed fine sigma. ",
      "Reduce fine sigma or use a larger image."
    )
  }
  pmax(blur(m, p$sigma1) - blur(m, coarse_sigma), 0)
}

local_contrast <- function(lab, m, ring_px = 4L) {
  n <- max(lab)
  out <- data.frame(
    candidate = seq_len(n),
    local_bg_median = rep(NA_real_, n),
    local_bg_mad = rep(NA_real_, n),
    contrast_z = rep(NA_real_, n)
  )
  
  for (j in seq_len(n)) {
    object <- lab == j
    ring <- morph(
      EBImage::dilate, object, 2L * ring_px + 1L
    ) > 0 & !object
    
    # Legacy behavior retained:
    # neighboring objects are not excluded from the background ring.
    bg <- m[ring]
    if (length(bg) >= 5L) {
      med <- median(bg)
      scale <- max(mad(bg, constant = 1.4826), 1e-9)
      out$local_bg_median[j] <- med
      out$local_bg_mad[j] <- scale
      out$contrast_z[j] <- (mean(m[object]) - med) / scale
    }
  }
  out
}

noise_scale <- function(m, p) {
  lab <- candidate_labels(dog_response(m, p), 2)
  lc <- local_contrast(lab, m)
  scale <- finite_median(lc$local_bg_mad)
  
  if (!is.finite(scale)) {
    warning(
      "Noise estimate unavailable; legacy fallback scale of 1 used.",
      call. = FALSE
    )
    return(1)
  }
  scale
}

segment_channel <- function(m, cfg) {
  p <- cfg$params
  
  if (cfg$method == "tophat") {
    sm <- blur(m, 1)
    bg <- morph(EBImage::opening, sm, p$brush_size)
    lab <- candidate_labels(pmax(sm - bg, 0), p$thr_mult)
  } else if (cfg$method == "dog") {
    lab <- candidate_labels(dog_response(m, p), p$thr_mult)
  } else {
    stop("Unknown segmentation method: ", cfg$method)
  }
  
  n <- max(lab)
  if (!n) return(lab)
  
  areas <- tabulate(lab, nbins = n)
  keep <- areas >= p$min_area & areas <= p$max_area
  
  if (cfg$method == "dog") {
    moments <- EBImage::computeFeatures.moment(eb(lab))
    if (!"m.eccentricity" %in% colnames(moments)) {
      stop("Eccentricity feature unavailable.")
    }
    
    lc <- local_contrast(lab, m)
    means <- vapply(seq_len(n), function(j) mean(m[lab == j]), numeric(1))
    ecc <- moments[, "m.eccentricity"]
    
    if (length(ecc) != n) stop("Unexpected moment-feature dimensions.")
    
    keep <- keep &
      means >= mean(m) + p$global_mult * sd(m) &
      is.finite(ecc) & ecc <= p$max_eccentricity &
      is.finite(lc$contrast_z) & lc$contrast_z >= p$contrast_z_min
  }
  
  relabel_keep(lab, which(keep))
}

# ---- Coordinate-safe measurements ------------------------------------------

peak_coords <- function(lab, m) {
  n <- max(lab)
  out <- matrix(
    numeric(0), nrow = 0L, ncol = 2L,
    dimnames = list(NULL, c("m.cx", "m.cy"))
  )
  if (!n) return(out)
  
  out <- t(vapply(seq_len(n), function(j) {
    pixels <- which(lab == j)
    best <- pixels[which.max(m[pixels])]
    xy <- arrayInd(best, dim(m))
    c(m.cx = xy[1, 1], m.cy = xy[1, 2])
  }, numeric(2)))
  colnames(out) <- c("m.cx", "m.cy")
  out
}

inside_mask <- function(mask, pts) {
  if (!nrow(pts)) return(logical(0))
  xy <- round(pts[, c("m.cx", "m.cy"), drop = FALSE])
  valid <- xy[, 1] >= 1 & xy[, 1] <= nrow(mask) &
    xy[, 2] >= 1 & xy[, 2] <= ncol(mask)
  
  out <- rep(FALSE, nrow(xy))
  out[valid] <- mask[xy[valid, , drop = FALSE]] > 0
  out
}

nearest_info <- function(pts) {
  n <- nrow(pts)
  if (n < 2L) {
    return(list(distance = numeric(0), to = integer(0)))
  }
  
  # Row-wise calculation avoids a dense n-by-n distance matrix.
  distance <- numeric(n)
  to <- integer(n)
  
  for (j in seq_len(n)) {
    ds <- (pts[, 1] - pts[j, 1])^2 + (pts[, 2] - pts[j, 2])^2
    ds[j] <- Inf
    to[j] <- which.min(ds)
    distance[j] <- sqrt(ds[to[j]])
  }
  
  list(distance = distance, to = to)
}

measure_channel <- function(m, lab, cfg) {
  n <- max(lab)
  pk <- peak_coords(lab, m)
  nn <- nearest_info(pk)
  
  radial <- if (n) {
    sqrt(
      (pk[, 1] - (nrow(m) + 1) / 2)^2 +
        (pk[, 2] - (ncol(m) + 1) / 2)^2
    )
  } else numeric(0)
  
  if (n) {
    mom <- EBImage::computeFeatures.moment(eb(lab))
    lc <- local_contrast(lab, m)
    areas <- tabulate(lab, nbins = n)
    
    objects <- data.frame(
      channel = cfg$name,
      object_id = seq_len(n),
      area_px = areas,
      peak_x = pk[, 1],
      peak_y = pk[, 2],
      centroid_x = mom[, "m.cx"],
      centroid_y = mom[, "m.cy"],
      eccentricity = mom[, "m.eccentricity"],
      mean_intensity = vapply(
        seq_len(n), function(j) mean(m[lab == j]), numeric(1)
      ),
      peak_intensity = vapply(
        seq_len(n), function(j) max(m[lab == j]), numeric(1)
      ),
      local_bg_median = lc$local_bg_median,
      local_bg_mad = lc$local_bg_mad,
      contrast_z = lc$contrast_z,
      nnd_px = if (n >= 2L) nn$distance else NA_real_,
      radial_px = radial,
      stringsAsFactors = FALSE
    )
  } else {
    objects <- data.frame()
    areas <- numeric(0)
    mom <- NULL
  }
  
  list(
    name = cfg$name,
    color = cfg$color,
    mat = m,
    lab = lab,
    n = n,
    size = finite_median(areas),
    elong = if (n) finite_median(mom[, "m.eccentricity"]) else NA_real_,
    ring = cv_value(nn$distance, minimum = 4L),
    radcv = cv_value(radial, minimum = 3L),
    pk = pk,
    nn = nn,
    radial = radial,
    objects = objects
  )
}
# ---- MatCol-inspired segmented-object overlap -------------------------------
# Uses existing labeled masks without dilation.
# Background = 0; positive integers identify individual objects.
# Shared pixels are counted separately for each A-B object pair.

validate_overlap_min <- function(x) {
  if (!is.numeric(x) || length(x) != 1L ||
      !is.finite(x) || x < 1 || x != floor(x) ||
      x > .Machine$integer.max) {
    stop("Minimum shared pixels must be a positive whole number.")
  }
  as.integer(x)
}

segmented_overlap <- function(lab_a, lab_b, min_shared_px = 3L) {
  min_shared_px <- validate_overlap_min(min_shared_px)
  
  validate_labels <- function(lab, channel) {
    if (!is.matrix(lab) || !is.numeric(lab) ||
        !length(lab) || any(!is.finite(lab)) ||
        any(lab < 0) || any(lab != floor(lab))) {
      stop(channel, ": expected a non-negative integer label matrix.")
    }
  }
  
  validate_labels(lab_a, "A")
  validate_labels(lab_b, "B")
  
  if (!identical(dim(lab_a), dim(lab_b))) {
    stop("Segmented-overlap masks must have matching dimensions.")
  }
  
  ids_a <- sort(unique(as.vector(lab_a[lab_a > 0])))
  ids_b <- sort(unique(as.vector(lab_b[lab_b > 0])))
  
  n_a <- length(ids_a)
  n_b <- length(ids_b)
  
  # All object pairs sharing at least one pixel, including pairs
  # subsequently rejected by the minimum-overlap threshold.
  pair_table <- data.frame(
    a_object_id = numeric(),
    b_object_id = numeric(),
    shared_px = integer(),
    qualifies = logical()
  )
  
  qualifying_mask <- matrix(FALSE, nrow(lab_a), ncol(lab_a))
  shared_index <- which(lab_a > 0 & lab_b > 0)
  
  if (length(shared_index)) {
    a_id <- lab_a[shared_index]
    b_id <- lab_b[shared_index]
    
    # Group actual intersecting pairs rather than allocating a dense
    # all-objects-by-all-objects contingency table.
    pair_key <- paste(a_id, b_id, sep = ":")
    unique_keys <- unique(pair_key)
    group <- match(pair_key, unique_keys)
    first <- match(unique_keys, pair_key)
    
    pair_table <- data.frame(
      a_object_id = a_id[first],
      b_object_id = b_id[first],
      shared_px = tabulate(group, nbins = length(unique_keys))
    )
    
    pair_table$qualifies <-
      pair_table$shared_px >= min_shared_px
    
    # Highlight only pixels belonging to qualifying object pairs.
    qualifying_mask[shared_index] <-
      pair_table$qualifies[group]
  }
  
  qualifying_pairs <- pair_table[
    pair_table$qualifies, , drop = FALSE
  ]
  
  overlapping_a <- unique(qualifying_pairs$a_object_id)
  overlapping_b <- unique(qualifying_pairs$b_object_id)
  
  list(
    segcoloc_pairs_n = nrow(qualifying_pairs),
    
    segcoloc_b_to_a_pct = if (n_b > 0L) {
      100 * length(overlapping_b) / n_b
    } else {
      NA_real_
    },
    
    segcoloc_a_to_b_pct = if (n_a > 0L) {
      100 * length(overlapping_a) / n_a
    } else {
      NA_real_
    },
    
    min_shared_px = min_shared_px,
    pair_table = pair_table,
    qualifying_mask = qualifying_mask,
    overlapping_a_ids = overlapping_a,
    overlapping_b_ids = overlapping_b,
    total_a = n_a,
    total_b = n_b
  )
}

pair_measurements <- function(a, b, tissue) {
  adil <- morph(EBImage::dilate, a$lab > 0, 7L) > 0
  bdil <- morph(EBImage::dilate, b$lab > 0, 7L) > 0
  
  # Legacy domain and empty-domain denominator retained.
  t_area <- max(sum(tissue), 1L)
  chance_a <- 100 * sum(adil & tissue) / t_area
  chance_b <- 100 * sum(bdil & tissue) / t_area
  
  b_inside <- inside_mask(adil, b$pk)
  a_inside <- inside_mask(bdil, a$pk)
  
  gross_ba <- if (b$n) 100 * mean(b_inside) else NA_real_
  gross_ab <- if (a$n) 100 * mean(a_inside) else NA_real_
  
  valid <- a$n > 0L && b$n > 0L
  
  # Preserve the original requirement that BOTH channels have objects
  # before exporting pairwise measurements.
  list(
    green_app_net = if (valid) gross_ba - chance_a else NA_real_,
    red_app_net_norm = if (valid) {
      (gross_ab - chance_b) / b$n
    } else NA_real_,
    sizeratio = if (valid && b$size > 0) a$size / b$size else NA_real_,
    gross_ba = gross_ba,
    gross_ab = gross_ab,
    chance_a = chance_a,
    chance_b = chance_b,
    a_inside = a_inside,
    b_inside = b_inside,
    adil = adil,
    bdil = bdil
  )
}

extract_sample <- function(
    sample, configs, reference = NULL, min_shared_px = 3L
) {
  min_shared_px <- validate_overlap_min(min_shared_px)
  mats <- lapply(seq_along(configs), function(i) {
    tryCatch(
      read_channel(sample$paths[i], configs[[i]]$src),
      error = function(e) {
        stop(configs[[i]]$name, ": ", conditionMessage(e), call. = FALSE)
      }
    )
  })
  
  dims <- lapply(mats, dim)
  if (!all(vapply(dims, identical, logical(1), dims[[1]]))) {
    stop("Channel dimensions differ; images must have matching crops.")
  }
  
  effective <- configs
  
  for (i in seq_along(configs)) {
    if (configs[[i]]$method == "dog" && !is.null(reference)) {
      ref <- reference[i]
      if (!is.finite(ref) || ref <= 0) {
        stop("Valid calibration reference missing for ", configs[[i]]$name)
      }
      
      scale <- noise_scale(mats[[i]], configs[[i]]$params)
      effective[[i]]$params$contrast_z_min <-
        configs[[i]]$params$contrast_z_min * scale / ref
    }
  }
  
  channels <- lapply(seq_along(mats), function(i) {
    tryCatch({
      lab <- segment_channel(mats[[i]], effective[[i]])
      measure_channel(mats[[i]], lab, effective[[i]])
    }, error = function(e) {
      stop(
        effective[[i]]$name, ": ", conditionMessage(e),
        call. = FALSE
      )
    })
  })
  
  feats <- list()
  for (cf in channels) {
    prefix <- make.names(cf$name)
    for (field in c("n", "size", "elong", "ring", "radcv")) {
      feats[[paste0(prefix, "_", field)]] <- cf[[field]]
    }
  }
  
  tissue <- Reduce("|", lapply(channels, function(cf) {
    cf$mat > mean(cf$mat) + 0.5 * sd(cf$mat)
  }))
  
  pairs <- list()
  if (length(channels) >= 2L) {
    indices <- combn(seq_along(channels), 2, simplify = FALSE)
    for (ij in indices) {
      i <- ij[1]
      j <- ij[2]
      a <- channels[[i]]
      b <- channels[[j]]
      # Existing proximity measuremnets:
      pm <- pair_measurements(a, b, tissue)
      
      # Additional segmented-object overlap, without dilation:
      overlap <- segmented_overlap(
        a$lab, b$lab,
        min_shared_px = min_shared_px
      )
      
      pm$segcoloc_pairs_n <- overlap$segcoloc_pairs_n
      pm$segcoloc_b_to_a_pct <- overlap$segcoloc_b_to_a_pct
      pm$segcoloc_a_to_b_pct <- overlap$segcoloc_a_to_b_pct
      pm$segcoloc <- overlap
      
      key <- paste(i, j, sep = "_")
      pairs[[key]] <- pm
      
      nm_a <- make.names(a$name)
      nm_b <- make.names(b$name)
      prefix <- paste0(nm_a, "_vs_", nm_b)
      
      # Marker-name-based feature names (flexible for any channel setup):
      #   <B>_near_<A>_net       (was _green_app_net)
      #   <A>_near_<B>_net_norm  (was _red_app_net_norm)
      feats[[paste0(prefix, "_", nm_b, "_near_", nm_a, "_net")]] <-
        pm$green_app_net
      feats[[paste0(prefix, "_", nm_a, "_near_", nm_b, "_net_norm")]] <-
        pm$red_app_net_norm
      
      for (field in c(
        "sizeratio",
        "segcoloc_pairs_n",
        "segcoloc_b_to_a_pct",
        "segcoloc_a_to_b_pct"
      )) {
        feats[[paste0(prefix, "_", field)]] <- pm[[field]]
      }
    }
  }
  
  list(
    feats = as.data.frame(feats, check.names = FALSE),
    channels = channels,
    pairs = pairs,
    width = nrow(mats[[1]]),
    height = ncol(mats[[1]]),
    effective = effective,
    objects = bind_rows(lapply(channels, `[[`, "objects"))
  )
}

capture_analysis <- function(expr) {
  warnings <- character(0)
  error <- NULL
  
  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }
  )
  
  list(
    result = result,
    error = error,
    warnings = paste(unique(warnings), collapse = " | ")
  )
}

# ---- Pairing manifest -------------------------------------------------------

trailing_id <- function(filename) {
  base <- tools::file_path_sans_ext(basename(filename))
  hit <- regmatches(base, regexpr("[0-9]+$", base))
  if (!length(hit) || !nzchar(hit)) return(NA_character_)
  # Character IDs avoid integer overflow and normalise leading zeros.
  hit <- sub("^0+", "", hit)
  if (!nzchar(hit)) "0" else hit
}

build_manifest <- function(conditions, configs) {
  rows <- list()
  samples <- list()
  nch <- length(configs)
  
  for (cond in conditions) {
    tables <- lapply(seq_len(nch), function(i) {
      u <- cond$uploads[[i]]
      if (is.null(u)) {
        return(data.frame(
          name = character(), path = character(), id = character()
        ))
      }
      u$id <- vapply(u$name, trailing_id, character(1))
      u
    })
    
    # Show every unpairable filename rather than silently dropping it.
    for (i in seq_len(nch)) {
      invalid <- which(is.na(tables[[i]]$id))
      for (j in invalid) {
        rows[[length(rows) + 1L]] <- data.frame(
          condition = cond$name,
          sample = NA_character_,
          channel = configs[[i]]$name,
          filename = tables[[i]]$name[j],
          status = "invalid",
          detail = "No trailing numeric sample ID",
          stringsAsFactors = FALSE
        )
      }
    }
    
    ids <- unique(unlist(lapply(tables, function(x) x$id)))
    ids <- sort(ids[!is.na(ids)])
    
    for (id in ids) {
      hits <- lapply(tables, function(x) which(!is.na(x$id) & x$id == id))
      counts <- lengths(hits)
      valid <- all(counts == 1L)
      
      reason <- if (any(counts > 1L)) {
        "Duplicate numeric ID within channel"
      } else if (any(counts == 0L)) {
        "Missing channel file"
      } else {
        "Ready"
      }
      
      paths <- character(nch)
      filenames <- character(nch)
      
      for (i in seq_len(nch)) {
        h <- hits[[i]]
        names_here <- if (length(h)) tables[[i]]$name[h] else ""
        
        rows[[length(rows) + 1L]] <- data.frame(
          condition = cond$name,
          sample = paste0("n", id),
          channel = configs[[i]]$name,
          filename = paste(names_here, collapse = " | "),
          status = if (valid) "ready" else "excluded",
          detail = reason,
          stringsAsFactors = FALSE
        )
        
        if (valid) {
          paths[i] <- tables[[i]]$path[h]
          filenames[i] <- tables[[i]]$name[h]
        }
      }
      
      if (valid) {
        samples[[length(samples) + 1L]] <- list(
          condition = cond$name,
          condition_id = cond$id,
          stem = paste0("n", id),
          paths = paths,
          filenames = filenames
        )
      }
    }
  }
  
  list(table = bind_rows(rows), samples = samples)
}

# ---- Plotting: convert [x,y] to raster [row,column] ---------------------------

normalise_display <- function(m) {
  r <- range(m)
  if (diff(r) <= 1e-12) return(m * 0)
  (m - r[1]) / diff(r)
}

make_raster <- function(res, outlines = FALSE) {
  # Raster rows correspond to y, raster columns to x.
  arr <- array(0, c(res$height, res$width, 3L))
  
  for (cf in res$channels) {
    intensity <- t(normalise_display(cf$mat))
    rgb <- as.numeric(grDevices::col2rgb(cf$color)) / 255
    
    for (k in 1:3) {
      arr[, , k] <- pmin(
        1, arr[, , k] + intensity * rgb[k] *
          if (outlines) 0.35 else 1
      )
    }
  }
  
  if (outlines) {
    for (cf in res$channels) {
      mask <- cf$lab > 0
      edge <- t(mask & !(morph(EBImage::erode, mask, 3L) > 0))
      rgb <- as.numeric(grDevices::col2rgb(cf$color)) / 255
      
      for (k in 1:3) {
        plane <- arr[, , k]
        plane[edge] <- rgb[k]
        arr[, , k] <- plane
      }
    }
  }
  
  grDevices::as.raster(arr)
}

empty_plot <- function(message = "Run a preview or select an extracted image.") {
  plot.new()
  text(0.5, 0.5, message, col = "#555555", cex = 0.9)
}

plot_mask <- function(mask, color, alpha = 0.25) {
  xy <- which(mask, arr.ind = TRUE)
  if (!nrow(xy)) return(invisible(NULL))
  
  # Deterministic thinning: no random redraw changes.
  if (nrow(xy) > 6000L) {
    xy <- xy[
      unique(round(seq(1, nrow(xy), length.out = 6000L))),
      , drop = FALSE
    ]
  }
  points(
    xy[, 1], xy[, 2], pch = 15, cex = 0.35,
    col = grDevices::adjustcolor(color, alpha.f = alpha)
  )
}

draw_image <- function(res, mode = "raw", pair = c(1L, 2L), channel = 1L) {
  par(mar = c(2, 0.5, 2.5, 0.5), bg = "white")
  
  if (is.null(res)) {
    empty_plot()
    return(invisible(NULL))
  }
  
  plot(
    NA, xlim = c(0.5, res$width + 0.5),
    ylim = c(res$height + 0.5, 0.5),
    type = "n", axes = FALSE, xlab = "", ylab = "",
    asp = 1, xaxs = "i", yaxs = "i"
  )
  
  # With reversed y limits, the raster bottom is at height + 0.5.
  rasterImage(
    make_raster(res, outlines = mode == "seg"),
    0.5, res$height + 0.5, res$width + 0.5, 0.5,
    interpolate = FALSE
  )
  
  if (mode %in% c("raw", "seg")) {
    title(if (mode == "raw") "Display-normalized composite" else "Segmentation")
    return(invisible(NULL))
  }
  
  channel <- min(max(channel, 1L), length(res$channels))
  cf <- res$channels[[channel]]
  
  if (mode == "nnd" && length(cf$nn$to)) {
    points_from <- cf$pk
    points_to <- cf$pk[cf$nn$to, , drop = FALSE]
    segments(
      points_from[, 1], points_from[, 2],
      points_to[, 1], points_to[, 2], col = cf$color
    )
  }
  
  if (mode == "radial") {
    cx <- (res$width + 1) / 2
    cy <- (res$height + 1) / 2
    if (cf$n) {
      segments(cx, cy, cf$pk[, 1], cf$pk[, 2], col = cf$color)
    }
    points(cx, cy, pch = 3, col = "yellow", cex = 1.5, lwd = 2)
  }
  
  caption <- ""
  pair_modes <- c("green_app_net", "red_app_net_norm")
  
  if (mode %in% pair_modes && length(res$channels) >= 2L) {
    key <- paste(pair, collapse = "_")
    pm <- res$pairs[[key]]
    
    if (!is.null(pm)) {
      a <- res$channels[[pair[1]]]
      b <- res$channels[[pair[2]]]
      
      {
        reverse <- mode == "red_app_net_norm"
        target <- if (reverse) b else a
        source <- if (reverse) a else b
        mask <- if (reverse) pm$bdil else pm$adil
        inside <- if (reverse) pm$a_inside else pm$b_inside
        
        plot_mask(mask, target$color)
        
        if (source$n) {
          points(
            source$pk[, 1], source$pk[, 2],
            pch = ifelse(inside, 21, 4),
            bg = source$color, col = "white", cex = 1.2
          )
        }
        
        caption <- if (reverse) {
          sprintf(
            "%s -> %s: %.4f percentage points / %s object",
            a$name, b$name, pm$red_app_net_norm, b$name
          )
        } else {
          sprintf(
            "%s -> %s: %.3f percentage points",
            b$name, a$name, pm$green_app_net
          )
        }
      }
    }
  } else {
    for (ch in res$channels) {
      if (ch$n) {
        points(
          ch$pk[, 1], ch$pk[, 2], pch = 21,
          bg = ch$color, col = "white", cex = 1.1
        )
      }
    }
    caption <- switch(
      mode,
      nnd = sprintf("%s: NND CV = %.3f", cf$name, cf$ring),
      radial = sprintf("%s: radial CV = %.3f", cf$name, cf$radcv),
      paste(vapply(
        res$channels, function(ch) paste0(ch$name, ": n=", ch$n),
        character(1)
      ), collapse = " | ")
    )
  }
  
  title(names(OVERLAYS)[match(mode, OVERLAYS)], cex.main = 0.95)
  mtext(caption, side = 1, line = 0.5, cex = 0.7)
}

draw_context <- function(res, mode, pair, channel) {
  par(mar = c(5, 4, 3, 1))
  if (is.null(res)) {
    empty_plot()
    return(invisible(NULL))
  }
  
  cf <- res$channels[[min(channel, length(res$channels))]]
  
  if (mode %in% c("nnd", "radial")) {
    values <- if (mode == "nnd") cf$nn$distance else cf$radial
    if (!length(values)) {
      empty_plot("Insufficient objects for this distance plot.")
      return(invisible(NULL))
    }
    barplot(
      sort(values), col = cf$color, border = NA,
      ylim = c(0, max(values, 1) * 1.15),
      xlab = "Object rank", ylab = "Distance (pixels)",
      main = paste(cf$name, if (mode == "nnd") "nearest neighbor" else "radial")
    )
    abline(h = mean(values), lty = 2)
    return(invisible(NULL))
  }
  
  if (mode == "peaks") {
    counts <- vapply(res$channels, `[[`, numeric(1), "n")
    barplot(
      counts,
      names.arg = vapply(res$channels, `[[`, character(1), "name"),
      col = vapply(res$channels, `[[`, character(1), "color"),
      ylim = c(0, max(counts, 1) * 1.2), ylab = "Object count"
    )
    return(invisible(NULL))
  }
  
  pm <- res$pairs[[paste(pair, collapse = "_")]]
  if (is.null(pm)) {
    empty_plot("Select a valid channel pair.")
    return(invisible(NULL))
  }
  {
    reverse <- mode == "red_app_net_norm"
    gross <- if (reverse) pm$gross_ab else pm$gross_ba
    chance <- if (reverse) pm$chance_b else pm$chance_a
    net <- if (reverse) {
      pm$red_app_net_norm * res$channels[[pair[2]]]$n
    } else pm$green_app_net
    
    values <- c(Gross = gross, Chance = chance, Net = net)
    ylabel <- "Gross/chance: percent\nNet: percentage points"
  }
  
  finite <- values[is.finite(values)]
  if (!length(finite)) {
    empty_plot("Metric unavailable.")
    return(invisible(NULL))
  }
  
  limits <- range(c(0, finite))
  padding <- max(diff(limits) * 0.15, 1)
  barplot(
    values, col = c("#0072B2", "#999999", "#E69F00"),
    ylim = limits + c(-padding, padding),
    ylab = ylabel, main = "Selected channel pair"
  )
  abline(h = 0)
}

# ---- Basic self-tests -------------------------------------------------------

run_coordinate_tests <- function() {
  m <- matrix(-5, nrow = 13L, ncol = 8L)
  lab <- matrix(0L, 13L, 8L)
  lab[10:11, 2:3] <- 1L
  m[10, 3] <- -1
  
  p <- peak_coords(lab, m)
  stopifnot(
    p[1, "m.cx"] == 10,
    p[1, "m.cy"] == 3,
    inside_mask(lab > 0, p),
    !inside_mask(lab > 0, matrix(
      c(3, 7), nrow = 1,
      dimnames = list(NULL, c("m.cx", "m.cy"))
    ))
  )
  
  points <- rbind(c(1, 1), c(4, 5))
  nn <- nearest_info(points)
  stopifnot(all(nn$distance == 5))
  stopifnot(trailing_id("channel_Gn003.tif") == "3")
  stopifnot(is.na(trailing_id("channel.tif")))
  invisible(TRUE)
}

run_coordinate_tests()
# ---- Segmented-overlap display additions ------------------------------------

OVERLAYS <- c(
  OVERLAYS,
  "Segmented-object overlap" = "segcoloc"
)

draw_image_before_segcoloc <- draw_image

draw_image <- function(
    res, mode = "raw", pair = c(1L, 2L), channel = 1L
) {
  if (mode != "segcoloc") {
    return(draw_image_before_segcoloc(
      res, mode = mode, pair = pair, channel = channel
    ))
  }
  
  par(mar = c(3.5, 0.5, 2.5, 0.5), bg = "white")
  
  if (is.null(res)) {
    empty_plot()
    return(invisible(NULL))
  }
  
  pm <- res$pairs[[paste(pair, collapse = "_")]]
  
  if (is.null(pm) || is.null(pm$segcoloc)) {
    empty_plot("Select a valid channel pair and rerun analysis.")
    return(invisible(NULL))
  }
  
  overlap <- pm$segcoloc
  a <- res$channels[[pair[1]]]
  b <- res$channels[[pair[2]]]
  
  # Existing composite and segmentation outlines.
  raster <- as.matrix(make_raster(res, outlines = TRUE))
  
  # EBImage labels use [x,y]; raster matrices use [y,x].
  # Replace every qualifying shared pixel with white.
  raster[t(overlap$qualifying_mask)] <- "#FFFFFF"
  
  plot(
    NA,
    xlim = c(0.5, res$width + 0.5),
    ylim = c(res$height + 0.5, 0.5),
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    asp = 1,
    xaxs = "i",
    yaxs = "i"
  )
  
  rasterImage(
    grDevices::as.raster(raster),
    0.5, res$height + 0.5,
    res$width + 0.5, 0.5,
    interpolate = FALSE
  )
  
  title("Segmented-object overlap", cex.main = 0.95)
  
  mtext(
    sprintf(
      "White: qualifying shared pixels | minimum = %d | pairs = %d",
      overlap$min_shared_px,
      overlap$segcoloc_pairs_n
    ),
    side = 1, line = 0.5, cex = 0.65
  )
  
  mtext(
    sprintf(
      "%s -> %s: %.1f%% | %s -> %s: %.1f%%",
      b$name, a$name, overlap$segcoloc_b_to_a_pct,
      a$name, b$name, overlap$segcoloc_a_to_b_pct
    ),
    side = 1, line = 1.5, cex = 0.65
  )
  
  invisible(NULL)
}

draw_context_before_segcoloc <- draw_context

draw_context <- function(res, mode, pair, channel) {
  if (mode != "segcoloc") {
    return(draw_context_before_segcoloc(
      res, mode = mode, pair = pair, channel = channel
    ))
  }
  
  par(mar = c(5, 4, 3, 1))
  
  if (is.null(res)) {
    empty_plot()
    return(invisible(NULL))
  }
  
  pm <- res$pairs[[paste(pair, collapse = "_")]]
  
  if (is.null(pm) || is.null(pm$segcoloc)) {
    empty_plot("Select a valid channel pair and rerun analysis.")
    return(invisible(NULL))
  }
  
  overlap <- pm$segcoloc
  a <- res$channels[[pair[1]]]
  b <- res$channels[[pair[2]]]
  
  values <- c(
    overlap$segcoloc_b_to_a_pct,
    overlap$segcoloc_a_to_b_pct
  )
  
  labels <- c(
    paste(b$name, "->", a$name),
    paste(a$name, "->", b$name)
  )
  
  # Use a zero-height placeholder to display NA explicitly.
  plotted_values <- values
  plotted_values[!is.finite(plotted_values)] <- 0
  
  positions <- barplot(
    plotted_values,
    names.arg = labels,
    col = c(b$color, a$color),
    border = NA,
    ylim = c(0, 115),
    ylab = "Source objects with qualifying overlap (%)",
    main = sprintf(
      "Segmented overlap: %d qualifying pairs",
      overlap$segcoloc_pairs_n
    )
  )
  
  text(
    positions,
    plotted_values + 5,
    labels = ifelse(
      is.finite(values),
      sprintf("%.1f%%", values),
      "NA: no source objects"
    ),
    cex = 0.8
  )
  
  invisible(NULL)
}
# ---- Segmented-overlap self-tests -------------------------------------------

run_segmented_overlap_tests <- function() {
  empty <- matrix(0L, 10L, 10L)
  
  # 1. Exactly three shared pixels: qualifies.
  a <- empty
  b <- empty
  a[2:4, 2] <- 1L
  b[2:4, 2] <- 1L
  
  result <- segmented_overlap(a, b, 3L)
  
  stopifnot(
    result$segcoloc_pairs_n == 1L,
    result$segcoloc_a_to_b_pct == 100,
    result$segcoloc_b_to_a_pct == 100,
    sum(result$qualifying_mask) == 3L
  )
  
  # 2. Same overlap below a four-pixel threshold: rejected.
  result <- segmented_overlap(a, b, 4L)
  
  stopifnot(
    result$segcoloc_pairs_n == 0L,
    result$segcoloc_a_to_b_pct == 0,
    result$segcoloc_b_to_a_pct == 0,
    !any(result$qualifying_mask)
  )
  
  # 3. One A object overlaps two B objects.
  # A second A object and third B object do not overlap.
  a <- empty
  b <- empty
  
  a[2:7, 2] <- 1L
  a[2:3, 7] <- 2L
  
  b[2:3, 2] <- 1L
  b[6:7, 2] <- 2L
  b[8:9, 8] <- 3L
  
  result <- segmented_overlap(a, b, 2L)
  
  stopifnot(
    result$segcoloc_pairs_n == 2L,
    result$segcoloc_a_to_b_pct == 50,
    isTRUE(all.equal(
      result$segcoloc_b_to_a_pct,
      100 * 2 / 3
    )),
    sum(result$qualifying_mask) == 4L
  )
  
  # 4. Separate two-pixel overlaps must not be pooled.
  result <- segmented_overlap(a, b, 3L)
  
  stopifnot(
    result$segcoloc_pairs_n == 0L,
    !any(result$qualifying_mask)
  )
  
  # 5. A empty, B nonempty.
  result <- segmented_overlap(empty, b, 2L)
  
  stopifnot(
    result$segcoloc_pairs_n == 0L,
    is.na(result$segcoloc_a_to_b_pct),
    result$segcoloc_b_to_a_pct == 0
  )
  
  # 6. A nonempty, B empty.
  result <- segmented_overlap(a, empty, 2L)
  
  stopifnot(
    result$segcoloc_pairs_n == 0L,
    result$segcoloc_a_to_b_pct == 0,
    is.na(result$segcoloc_b_to_a_pct)
  )
  
  # 7. Both channels empty.
  result <- segmented_overlap(empty, empty, 3L)
  
  stopifnot(
    result$segcoloc_pairs_n == 0L,
    is.na(result$segcoloc_a_to_b_pct),
    is.na(result$segcoloc_b_to_a_pct)
  )
  
  # 8. Nonconsecutive labels must not inflate object counts.
  a <- empty
  b <- empty
  a[2:4, 2] <- 7L
  b[2:4, 2] <- 19L
  
  result <- segmented_overlap(a, b, 3L)
  
  stopifnot(
    result$total_a == 1L,
    result$total_b == 1L,
    result$segcoloc_pairs_n == 1L
  )
  
  # 9. Reject invalid thresholds and mismatched dimensions.
  stopifnot(
    inherits(
      try(segmented_overlap(a, b, 0), silent = TRUE),
      "try-error"
    ),
    inherits(
      try(segmented_overlap(a, b, 2.5), silent = TRUE),
      "try-error"
    ),
    inherits(
      try(
        segmented_overlap(a, matrix(0L, 5L, 5L), 3),
        silent = TRUE
      ),
      "try-error"
    )
  )
  
  message("Segmented-overlap self-tests passed.")
  invisible(TRUE)
}

run_segmented_overlap_tests()
# ---- Export helpers (Visualize tab) -----------------------------------------
GG_FONT <- "Helvetica"   # standard PDF font; no Cairo/XQuartz needed
FS_TXT  <- 11            # axis titles / base text (pt)
FS_TICK <- 10            # tick labels (pt)
LW      <- 0.5           # axis, box, error-bar line width
PT_SZ   <- 2.4           # point size
EXPORT_FORMATS <- c("PNG" = "png", "PDF" = "pdf")

open_export_device <- function(file, fmt, width_in, height_in, dpi) {
  if (identical(fmt, "pdf")) {
    grDevices::pdf(
      file, width = width_in, height = height_in, family = GG_FONT,
      paper = "special", onefile = FALSE, useDingbats = FALSE, bg = "white"
    )
  } else if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      file, width = width_in, height = height_in, units = "in",
      res = dpi, background = "white"
    )
  } else {
    grDevices::png(
      file, width = width_in, height = height_in, units = "in",
      res = dpi, bg = "white"
    )
  }
  grDevices::dev.cur()
}

# Returns a ggplot, or a character string explaining why there is no plot.
# pub = FALSE reproduces the on-screen styling; pub = TRUE applies figure styling.
build_distribution_plot <- function(run, feature, plot_type, stats, pub = FALSE) {
  if (is.null(run) || is.null(feature) || !feature %in% run$feature_cols) {
    return("No feature available for the current selection.")
  }
  df <- run$features
  df$value <- df[[feature]]
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) return("No finite values.")
  
  df$condition <- factor(df$condition, levels = unique(df$condition))
  palette <- setNames(
    grDevices::hcl.colors(nlevels(df$condition), "Dark 3"),
    levels(df$condition)
  )
  
  pt  <- if (pub) PT_SZ else 1.8
  lw  <- if (pub) LW else 0.5
  thm <- if (pub) {
    theme_classic(base_size = FS_TXT, base_family = GG_FONT) +
      theme(
        axis.title   = element_text(size = FS_TXT, colour = "black"),
        axis.text    = element_text(size = FS_TICK, colour = "black"),
        axis.line    = element_line(linewidth = LW),
        axis.ticks   = element_line(linewidth = LW),
        plot.caption = element_text(size = FS_TICK - 2, colour = "black"),
        legend.position = "none"
      )
  } else {
    theme_minimal(base_size = 13) + theme(legend.position = "none")
  }
  
  base <- ggplot(df, aes(condition, value, fill = condition)) +
    scale_fill_manual(values = palette) +
    labs(x = NULL, y = feature) + thm
  
  if (plot_type == "bar") {
    summary <- bind_rows(lapply(split(df, df$condition), function(part) {
      data.frame(
        condition = part$condition[1],
        mean = mean(part$value),
        se = if (nrow(part) >= 2L) sd(part$value) / sqrt(nrow(part)) else NA_real_
      )
    }))
    p <- ggplot(summary, aes(condition, mean, fill = condition)) +
      geom_col(width = 0.6, linewidth = lw) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                    width = 0.15, linewidth = lw, na.rm = TRUE) +
      scale_fill_manual(values = palette) +
      labs(x = NULL, y = paste("Mean +/- SE:", feature)) + thm
  } else if (plot_type == "violin") {
    counts <- table(df$condition)
    viable <- names(counts[counts >= 2L])
    p <- base +
      geom_violin(data = df[df$condition %in% viable, , drop = FALSE],
                  alpha = 0.35, scale = "width", linewidth = lw, na.rm = TRUE) +
      geom_point(position = position_jitter(width = 0.1, seed = 1), size = pt)
  } else {
    p <- base +
      geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = lw) +
      geom_point(position = position_jitter(width = 0.1, seed = 1), size = pt)
  }
  
  if (!is.null(stats)) {
    row <- stats[stats$feature == feature, , drop = FALSE]
    if (nrow(row) && is.finite(row$adj_p[1])) {
      p <- p + labs(caption = sprintf(
        "%s vs %s | %s | %s p = %s",
        row$condition_A[1], row$condition_B[1],
        row$test[1], row$correction[1],
        format.pval(row$adj_p[1], digits = 4, eps = 1e-5)
      ))
    }
  }
  p
}

# Writes a distribution plot to `file`.
export_distribution <- function(file, p, fmt, width_mm, height_mm, dpi) {
  dev <- open_export_device(file, fmt, width_mm / 25.4, height_mm / 25.4, dpi)
  on.exit(grDevices::dev.off(which = dev), add = TRUE)
  print(p)
}

# Clean image: no title/caption, one image pixel = one output pixel (PNG).
draw_clean_image <- function(res, outlines) {
  par(mar = c(0, 0, 0, 0), bg = "white")
  plot(NA, xlim = c(0.5, res$width + 0.5), ylim = c(res$height + 0.5, 0.5),
       type = "n", axes = FALSE, xlab = "", ylab = "",
       asp = 1, xaxs = "i", yaxs = "i")
  rasterImage(make_raster(res, outlines = outlines),
              0.5, res$height + 0.5, res$width + 0.5, 0.5,
              interpolate = FALSE)
}

# mode: "raw" / "seg" (clean = TRUE gives title-free output) or any overlay mode
# (annotated, same drawing code as the screen, with title and caption).
export_overlay <- function(file, res, mode, pair, channel, fmt,
                           width_mm, dpi, clean = FALSE) {
  if (is.null(res)) stop("No extracted image selected.")
  width_in <- width_mm / 25.4
  aspect   <- res$height / res$width
  
  if (clean && mode %in% c("raw", "seg")) {
    dev <- open_export_device(
      file, fmt, width_in, width_in * aspect,
      dpi = res$width / width_in           # 1 image px = 1 output px
    )
    on.exit(grDevices::dev.off(which = dev), add = TRUE)
    draw_clean_image(res, outlines = identical(mode, "seg"))
  } else {
    lh <- 0.2 * FS_TXT / 12                # one margin line, inches
    height_in <- (width_in - lh) * aspect + 6 * lh
    dev <- open_export_device(file, fmt, width_in, height_in, dpi)
    on.exit(grDevices::dev.off(which = dev), add = TRUE)
    par(ps = FS_TXT)
    draw_image(res, mode, pair, channel)
  }
  invisible(file)
}
# ---- UI --------------------------------------------------------------------

channel_card <- function(i) {
  default_names <- c("PGL-1", "MUT-16", "ALG-4", "PRG-1")
  default_colors <- c("Red", "Green", "Yellow", "Magenta")
  default_methods <- c("tophat", "dog", "dog", "tophat")
  
  conditionalPanel(
    condition = sprintf("Number(input.n_channels) >= %d", i),
    wellPanel(
      h4(paste("Channel", i)),
      fluidRow(
        column(3, textInput(
          paste0("ch_name_", i), "Marker name", default_names[i]
        )),
        column(3, selectInput(
          paste0("ch_src_", i), "RGB source",
          SOURCES, selected = c(1, 2, 3, 1)[i]
        )),
        column(3, selectInput(
          paste0("ch_color_", i), "Display color",
          names(COLORS), selected = default_colors[i]
        )),
        column(3, selectInput(
          paste0("ch_method_", i), "Segmentation",
          METHODS, selected = default_methods[i]
        ))
      )
    )
  )
}

parameter_card <- function(i) {
  conditionalPanel(
    condition = sprintf("Number(input.n_channels) >= %d", i),
    wellPanel(
      h4(paste("Channel", i, "parameters")),
      numericInput(paste0("p_thr_", i), "Threshold multiplier",
                   if (i %in% c(2, 3)) 2 else 1, min = 0, step = 0.1),
      numericInput(paste0("p_minarea_", i), "Minimum area (pixels)",
                   if (i %in% c(2, 3)) 20 else 10, min = 1),
      numericInput(paste0("p_maxarea_", i), "Maximum area (pixels)",
                   if (i %in% c(2, 3)) 1000 else 3000, min = 1),
      
      conditionalPanel(
        condition = sprintf("input.ch_method_%d == 'tophat'", i),
        numericInput(paste0("p_brush_", i), "Top-hat brush diameter",
                     51, min = 1, step = 2)
      ),
      conditionalPanel(
        condition = sprintf("input.ch_method_%d == 'dog'", i),
        numericInput(paste0("p_sigma1_", i), "Fine sigma", 1, min = 0.1),
        numericInput(paste0("p_sigma2_", i), "Coarse sigma", 8, min = 0.2),
        numericInput(paste0("p_globmult_", i), "Global floor multiplier",
                     1.3, min = 0, step = 0.05),
        numericInput(paste0("p_ecc_", i), "Maximum eccentricity",
                     0.90, min = 0, max = 1, step = 0.01),
        numericInput(paste0("p_zmin_", i), "Local contrast z minimum",
                     1, min = 0, step = 0.1)
      )
    )
  )
}

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body {
      font-family: 'Segoe UI', Arial, sans-serif;
      font-size: 15px;
      color: #24323d;
      background: #f5f7fa;
    }

    body > .container-fluid {
      padding: 20px 28px 32px;
    }

    /* Logo header and subtitle */
    body > .container-fluid > h2 {
      margin-top: 0;
      margin-bottom: 8px;
    }

    body > .container-fluid > p {
      color: #64748b;
      font-size: 15px;
      margin-bottom: 22px;
    }

    /* Navigation */
    #main_tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      padding: 7px;
      background: #eaf0f5;
      border: none;
      border-radius: 14px;
    }

    #main_tabs > li {
      float: none;
      margin: 0;
    }

    #main_tabs > li > a {
      margin: 0;
      padding: 10px 16px;
      border: none;
      border-radius: 10px;
      color: #526174;
      font-weight: 600;
      transition: background 0.18s ease, color 0.18s ease;
    }

    #main_tabs > li > a:hover {
      background: #ffffff;
      color: #0f766e;
      border: none;
    }

    #main_tabs > li.active > a,
    #main_tabs > li.active > a:hover,
    #main_tabs > li.active > a:focus {
      background: #0f766e;
      color: #ffffff;
      border: none;
      box-shadow: 0 3px 8px rgba(15, 118, 110, 0.18);
    }

    .tab-content {
      padding-top: 24px;
    }

    /* Cards */
    .well {
      background: #ffffff;
      border: 1px solid #e1e7ef;
      border-top: 3px solid #8b5cf6;
      border-radius: 14px;
      padding: 18px 20px;
      margin-bottom: 18px;
      box-shadow: 0 4px 14px rgba(31, 41, 55, 0.05);
    }

    .well h4 {
      margin-top: 0;
      margin-bottom: 16px;
      color: #24323d;
      font-size: 17px;
      font-weight: 700;
    }

    .well .form-group {
      margin-bottom: 8px;
    }

    /* Input labels and controls */
    .control-label {
      color: #475569;
      font-size: 13px;
      font-weight: 600;
      margin-bottom: 7px;
    }

    .form-control,
    .selectize-input {
      border: 1px solid #d7e0e9;
      border-radius: 8px;
      background: #ffffff;
      color: #24323d;
      font-size: 14px;
      box-shadow: none;
    }

    .form-control {
      min-height: 38px;
    }

    .selectize-input {
      min-height: 38px;
      padding: 8px 12px;
    }

    .selectize-input.full {
      background: #ffffff;
    }

    .form-control:focus,
    .selectize-input.focus {
      border-color: #14b8a6;
      box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.12);
    }

    .selectize-dropdown {
      border: 1px solid #d7e0e9;
      border-radius: 8px;
      overflow: hidden;
      box-shadow: 0 6px 18px rgba(31, 41, 55, 0.10);
    }

    .selectize-dropdown .active {
      background: #e6f5f2;
      color: #0f766e;
    }

    /* Buttons */
    .btn {
      padding: 9px 15px;
      margin-bottom: 6px;
      border-radius: 8px;
      font-size: 14px;
      font-weight: 600;
      transition: background 0.18s ease;
    }

    .btn-primary {
      background: #0f766e;
      border-color: #0f766e;
      color: #ffffff;
    }

    .btn-primary:hover,
    .btn-primary:focus,
    .btn-primary:active {
      background: #115e59;
      border-color: #115e59;
      color: #ffffff;
    }

    .btn-default {
      background: #ffffff;
      border-color: #d7e0e9;
      color: #475569;
    }

    .btn-default:hover {
      background: #f0f5f8;
      border-color: #b9c8d7;
    }

    /* Guidance and status boxes */
    .help-block {
      color: #718096;
      font-size: 13px;
      line-height: 1.55;
      margin-top: 8px;
      margin-bottom: 14px;
    }

    pre.shiny-text-output {
      background: #edf4f5;
      border: 1px solid #d9e7e8;
      border-left: 3px solid #14b8a6;
      border-radius: 8px;
      padding: 12px;
      color: #475569;
      font-size: 12px;
      white-space: pre-wrap;
      overflow-wrap: anywhere;
    }

    .shiny-output-error-validation {
      color: #8b4100;
    }

    /* Smaller screens */
    @media (max-width: 767px) {
      body > .container-fluid {
        padding: 12px;
      }

      #main_tabs > li > a {
        padding: 8px 11px;
        font-size: 13px;
      }

      .well {
        padding: 16px;
      }
    }
  "))),
  titlePanel(
    tagList(
      div(
        style = paste(
          "display: flex;",
          "align-items: center;",
          "justify-content: space-between;",
          "width: 100%;",
          "gap: 24px;"
        ),
        
        # Main GGProfiler logo on the left
        tags$img(
          src = "image.png",
          alt = "GGProfiler",
          style = paste(
            "width: 540px;",
            "max-width: 55%;",
            "height: auto;",
            "display: block;"
          )
        ),
        
        # Silva Garcia Lab logo (right), enlarged
        tags$img(
          src = "Lab Logo full.png",
          alt = "Silva Garcia Lab",
          style = paste(
            "height: 130px;",
            "width: auto;",
            "max-width: 40vw;",
            "display: block;",
            "margin-left: auto;",
            "flex-shrink: 0;"
          )
        )
      ),
      
      div(
        style = paste(
          "font-size: 23px;",
          "color: #6c7a89;",
          "margin-top: 12px;"
        ),
        "Multi-channel germ granule analysis and spatial profiling"
      )
    ),
    windowTitle = "GGProfiler"
  ),
  
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel(
      "Setup", value = "setup",
      fluidRow(
        column(
          3,
          selectInput(
            "n_channels", "Number of channels",
            1:4, selected = 2
          ),
          helpText(
            "Configure your markers, then click Apply setup. ",
            "Expand Setup instructions below for guidance."
          ),
          actionButton(
            "apply_setup", "Apply setup",
            class = "btn-primary"
          ),
          verbatimTextOutput("setup_status"),
          
          tags$details(
            style = "margin-top: 12px;",
            
            tags$summary(
              style = paste(
                "cursor: pointer;",
                "color: #0f766e;",
                "font-weight: 600;"
              ),
              "Setup instructions"
            ),
            
            div(
              style = "padding-top: 12px; line-height: 1.7;",
              
              h3("Setup instructions"),
              p(
                "Use Setup to define the markers you want to analyze ",
                "and tell GGProfiler where each marker's signal is stored. ",
                "Complete this step before adding conditions or uploading images."
              ),
              
              h4("1. Choose the number of channels"),
              p(
                "Select the number of markers you want to analyze, from 1 to 4. ",
                "This is the number of analysis channels, not the number of ",
                "nuclei, images or experimental conditions."
              ),
              p(
                "For example, to analyze PGL-1 and MUT-16, select 2 channels. ",
                "A separate settings card appears for each channel."
              ),
              
              h4("2. Enter a marker name"),
              p(
                "Enter the marker name in each channel card, for example ",
                "PGL-1 or MUT-16. Replace the default names if necessary. ",
                "Every channel must have a non-empty, distinct name."
              ),
              p(
                "These names identify your upload fields and exported measurements. ",
                "Avoid names that differ only by spaces or punctuation, because ",
                "R may convert them to the same name in the results."
              ),
              
              h4("3. Select the RGB source"),
              p(
                "RGB source determines which intensity signal is analyzed. ",
                "Choose it according to how your image files store the marker."
              ),
              
              tags$dl(
                tags$dt("Separate grayscale images"),
                tags$dd(
                  p(
                    "Recommended: use one 2D grayscale TIFF per marker and nucleus. ",
                    "The app uses the file's single intensity plane automatically. ",
                    "The RGB source selection is ignored for these files."
                  )
                ),
                
                tags$dt("Separate RGB color images"),
                tags$dd(
                  p(
                    "Choose R, G or B to match the plane containing the marker's ",
                    "signal. For example, choose R for a marker stored in the ",
                    "red plane and G for a marker stored in the green plane."
                  )
                ),
                
                tags$dt("Merged RGB/composite images"),
                tags$dd(
                  p(
                    "Configure one analysis channel for each marker. ",
                    "Choose the RGB plane containing that marker's signal: ",
                    "for example, R for PGL-1 and G for MUT-16 if those markers ",
                    "are stored in the red and green planes."
                  ),
                  p(
                    "Later, in Conditions, upload the same complete batch of ",
                    "composite images under each corresponding marker's upload ",
                    "field. You can select multiple images together; you do not ",
                    "need to upload each nucleus separately. Keep the filenames ",
                    "unchanged across the upload fields."
                  ),
                  p(
                    "This works only when the markers occupy distinct RGB planes. ",
                    "Do not assume an arbitrary pseudocolored composite preserves ",
                    "each marker in a separate plane. Use separate grayscale ",
                    "marker images when the signals are mixed across RGB planes."
                  )
                )
              ),
              h4("4. Choose a display color"),
              p(
                "Select the color used to display each marker in the app's ",
                "composite images and overlays. Choose distinct colors so ",
                "the markers are easy to distinguish."
              ),
              p(
                tags$strong("Display color does not select the input signal. "),
                "For example, you can analyze the R source and display it in cyan. ",
                "Changing the display color does not change the measured intensities."
              ),
              
              h4("5. Choose a segmentation method"),
              p(
                "Segmentation identifies candidate granules or objects in each ",
                "marker channel. Choose a starting method for each marker:"
              ),
              tags$ul(
                tags$li(
                  tags$strong("Top-hat: "),
                  "A starting option for larger objects against a varying ",
                  "background. The brush diameter and detection thresholds ",
                  "are adjusted in QC."
                ),
                tags$li(
                  tags$strong("DoG (Difference of Gaussians): "),
                  "A starting option for small, punctate objects. ",
                  "The fine and coarse scales, thresholds and object filters ",
                  "are adjusted in QC."
                )
              ),
              p(
                "Neither method is guaranteed to suit every image. ",
                "Check segmentation in QC and adjust the parameters before ",
                "running extraction. The two methods may be used for different ",
                "markers in the same analysis."
              ),
              
              h4("6. Apply the setup"),
              tags$ol(
                tags$li(
                  "Check the number of channels, marker names, RGB sources, ",
                  "display colors and segmentation methods."
                ),
                tags$li(
                  "Click Apply setup. The app saves these selections ",
                  "and switches to the Conditions tab."
                ),
                tags$li(
                  "Add your conditions and upload the images under the ",
                  "corresponding marker fields."
                )
              ),
              p(
                tags$strong("Changing the setup later: "),
                "Remove all existing conditions before applying changes to ",
                "the channel setup. Removing conditions removes their current ",
                "uploads from the analysis, so you will need to add the ",
                "conditions and upload their images again."
              ),
              
              h4("Image checklist before uploading"),
              tags$ul(
                tags$li(
                  "Use single 2D grayscale or RGB images, not Z-stacks ",
                  "or time-series stacks."
                ),
                tags$li(
                  "Images for different markers from the same nucleus must ",
                  "have matching dimensions, alignment and crop."
                ),
                tags$li(
                  "For separate marker images, use the same ending numeric ",
                  "ID for the same nucleus, for example PGL1_1.tif and MUT16_1.tif."
                ),
                tags$li(
                  "For composite images, upload the same batch under each ",
                  "relevant marker, for example Composite_1.tif, Composite_2.tif ",
                  "and Composite_3.tif."
                ),
                tags$li(
                  "Use a unique nucleus number within each condition, ",
                  "including across animals or imaging fields."
                )
              )
            )
          )
        ),
        
        column(
          9,
          lapply(seq_len(MAX_CHANNELS), channel_card)
        )
      )
    ),
    tabPanel(
      "Crop nuclei", value = "crop",
      tags$details(
        style = "margin-top: 12px; margin-bottom: 16px;",
        tags$summary(
          style = paste(
            "cursor: pointer;",
            "color: #0f766e;",
            "font-weight: 600;"
          ),
          "Crop nuclei instructions"
        ),
        div(
          style = "padding-top: 12px; line-height: 1.6;",
          
          h4("1. Purpose"),
          p(
            "Use this optional tab to cut individual nuclei out of a ",
            "larger field and create the matched per-marker images that ",
            "Conditions expects. Cropping uses raw pixel values; the ",
            "on-screen contrast stretch is for viewing only."
          ),
          
          h4("2. Upload the field images"),
          p(
            "Select every image of ONE field together: for example one ",
            "grayscale TIFF per marker, or one merged RGB image. All ",
            "selected images must have identical pixel dimensions and be ",
            "aligned. Images whose dimensions differ from the reference ",
            "image are skipped and reported in the crop log. This tab ",
            "reads TIFF, PNG and JPEG files; convert other formats to ",
            "TIFF first."
          ),
          
          h4("3. Draw and save a crop"),
          tags$ol(
            tags$li("Choose the reference image shown for drawing."),
            tags$li(
              "Optional zoom: use the Zoom and Pan sliders, or drag a box ",
              "around the area of interest and click Zoom to dragged box ",
              "(this sets the sliders for you). Reset zoom returns to the ",
              "full image. For multi-plane TIFFs, use the Z plane slider to ",
              "choose the plane shown. Zoom and plane change the view only, ",
              "not the pixel data."
            ),
            tags$li(
              "Drag a rectangle around one nucleus (in the zoomed view if ",
              "you zoomed), then click Apply rectangle."
            ),
            tags$li(
              "Enter a nucleus number and click Save crop. The same ",
              "rectangle is applied to every uploaded image. Repeat for ",
              "each nucleus; saved rectangles are outlined in cyan."
            )
          ),
          
          h4("4. Output file names"),
          p(
            "Each output is named from its source file plus the nucleus ",
            "number, for example ",
            tags$code("PGL1_field_1.tif"), " and ",
            tags$code("MUT16_field_1.tif"),
            ". The ending number is the nucleus ID used for pairing in ",
            "Conditions. Nucleus numbers must be unique within a ",
            "condition, so number crops consistently across fields of the ",
            "same condition."
          ),
          
          h4("5. Saving the crops"),
          p(
            "Two options: (a) Save to folder writes the cropped TIFFs ",
            "(16-bit) directly into the output folder you choose, plus a ",
            "timestamped crop_log_*.csv. Type or paste the folder path, or ",
            "click Browse... to pick it. This works only when the app runs ",
            "on your own computer. (b) Download cropped images (ZIP) ",
            "bundles the TIFFs and crop_log.csv. By default all planes of a ",
            "multi-plane TIFF are saved; tick 'Save only the displayed ",
            "plane (2D output)' to save just the plane shown. Existing files ",
            "are skipped unless 'Overwrite existing files' is ticked."
          ),
          
          h4("6. Reproducibility"),
          p(
            "The log records the source file, rectangle coordinates ",
            "(EBImage convention: x = column, y = row, 1-based, y = 1 is ",
            "the top row), planes saved and status of every crop. Keep it ",
            "with your data and cite the crop criterion in your Methods. ",
            "Cropping does not align images or detect nuclei automatically."
          ),
          
          h4("7. Continue"),
          p(
            "Upload the cropped files under the matching markers in ",
            "Conditions (unzip first if you used the ZIP)."
          )
        )
      ),
      fluidRow(
        column(
          3,
          fileInput(
            "crop_files",
            "Field images (select all channels of one field)",
            multiple = TRUE,
            accept = c(".tif", ".tiff", ".png", ".jpg", ".jpeg")
          ),
          selectInput("crop_ref", "Reference image for drawing", choices = character()),
          uiOutput("crop_plane_ui"),
          tags$div(
            style = "margin-bottom: 6px; font-weight: 600; color: #475569;",
            "Zoom"
          ),
          actionButton("crop_zoom", "Zoom to dragged box"),
          actionButton("crop_zoom_out", "Zoom out"),
          actionButton("crop_zoom_reset", "Reset zoom"),
          tags$hr(),
          actionButton("crop_apply", "Apply rectangle", class = "btn-primary"),
          actionButton("crop_clear", "Clear rectangle"),
          numericInput("crop_id", "Nucleus number", value = 1, min = 1, step = 1),
          actionButton("crop_save", "Save crop", class = "btn-primary"),
          actionButton("crop_remove_last", "Remove last"),
          actionButton("crop_clear_all", "Clear all"),
          tags$hr(),
          checkboxInput(
            "crop_single_plane",
            "Save only the displayed plane (2D output)", FALSE
          ),
          textInput("crop_out_dir", "Output folder", placeholder = "e.g. ~/Desktop/crops"),
          actionButton("crop_browse_dir", "Browse..."),
          checkboxInput("crop_overwrite", "Overwrite existing files", FALSE),
          actionButton("crop_save_dir", "Save to folder", class = "btn-primary"),
          helpText("Folder saving works when the app runs on your own computer. Otherwise use the ZIP."),
          tags$hr(),
          downloadButton("crop_dl_zip", "Download cropped images (ZIP)"),
          downloadButton("crop_dl_log", "Crop coordinates CSV")
        ),
        column(
          9,
          verbatimTextOutput("crop_status"),
          plotOutput(
            "crop_img",
            height = "520px",
            brush = brushOpts(
              id = "crop_brush",
              delay = 200,
              delayType = "debounce",
              resetOnNew = FALSE
            )
          ),
          h4("Preview of the current rectangle"),
          div(
            style = "width: 280px; max-width: 100%;",
            plotOutput("crop_preview", height = "260px")
          ),
          h4("Saved crops"),
          DTOutput("crop_table")
        )
      )
    ),
    
    tabPanel(
      "Conditions", value = "conditions",
      fluidRow(
        column(
          3,
          textInput("new_cond_name", "Condition name"),
          actionButton("add_cond", "Add condition", class = "btn-primary"),
          
          helpText(
            "Add a condition, then upload its images under each marker. ",
            "Expand Conditions instructions below for guidance."
          ),
          
          tags$details(
            style = "margin-top: 12px; margin-bottom: 16px;",
            
            tags$summary(
              style = paste(
                "cursor: pointer;",
                "color: #0f766e;",
                "font-weight: 600;"
              ),
              "Conditions instructions"
            ),
            
            div(
              style = "padding-top: 12px; line-height: 1.7;",
              
              h3("Conditions instructions"),
              p(
                "Use Conditions to organize your experimental groups and ",
                "upload the images belonging to each group. Complete Setup ",
                "and click Apply setup before adding conditions."
              ),
              
              h4("1. Add an experimental condition"),
              tags$ol(
                tags$li(
                  "Enter a clear condition name, such as Control or Treatment."
                ),
                tags$li(
                  "Click Add condition. A condition card appears on the right ",
                  "with a separate upload field for each marker defined in Setup."
                ),
                tags$li(
                  "Repeat for each experimental group. Each condition must ",
                  "have a different name."
                )
              ),
              p(
                "Upload only images belonging to that condition into its card. ",
                "You can include multiple nuclei in the same condition."
              ),
              
              h4("2. Prepare your images"),
              tags$ul(
                tags$li(
                  "Use single 2D grayscale or RGB images, not Z-stacks ",
                  "or time-series stacks. Separate grayscale TIFFs ",
                  "are recommended."
                ),
                tags$li(
                  "For this nucleus-based workflow, prepare a matching image ",
                  "crop for each marker from each nucleus."
                ),
                tags$li(
                  "Images from the same nucleus must have matching dimensions, ",
                  "alignment and crop. The app does not align images or ",
                  "identify nuclear boundaries automatically."
                )
              ),
              
              h4("3. Name files so the channels match"),
              p(
                tags$strong(
                  "Use a consistent naming pattern and the same ending ",
                  "numeric ID for all marker images from the same nucleus."
                ),
                " Place the nucleus number immediately before the file extension."
              ),
              tags$ul(
                tags$li(
                  "Nucleus 1: ",
                  tags$code("PGL1_1.tif"), " and ",
                  tags$code("MUT16_1.tif"), "."
                ),
                tags$li(
                  "Nucleus 2: ",
                  tags$code("PGL1_2.tif"), " and ",
                  tags$code("MUT16_2.tif"), "."
                )
              ),
              p(
                "The app matches files using this ending number, not the ",
                "order in which you select or upload them."
              ),
              tags$ul(
                tags$li(
                  "Use a unique nucleus number within each condition, ",
                  "including across different animals or imaging fields. ",
                  "Do not restart numbering for each animal in the same condition."
                ),
                tags$li(
                  "Provide exactly one file per marker for each nucleus number."
                ),
                tags$li(
                  "Numbers can be reused in different conditions because ",
                  "pairing is performed separately within each condition."
                ),
                tags$li(
                  "Leading zeros do not create different IDs: ",
                  tags$code("001"), " and ", tags$code("1"),
                  " are treated as the same nucleus number."
                )
              ),
              
              h4("4. Upload separate marker images"),
              tags$ol(
                tags$li(
                  "Under the first marker, click Browse and select all images ",
                  "for that marker in this condition. You can select multiple ",
                  "files together."
                ),
                tags$li(
                  "Under the next marker, select the corresponding images ",
                  "for the same nuclei."
                ),
                tags$li(
                  "Repeat for every marker configured in Setup."
                )
              ),
              p(
                "For example, upload ",
                tags$code("PGL1_1.tif"), " and ",
                tags$code("PGL1_2.tif"),
                " under PGL-1, then upload ",
                tags$code("MUT16_1.tif"), " and ",
                tags$code("MUT16_2.tif"),
                " under MUT-16."
              ),
              p(
                "For grayscale files, the single intensity plane is used ",
                "automatically. For separate RGB color files, select the ",
                "correct R, G or B source for each marker in Setup."
              ),
              
              h4("5. Upload merged RGB/composite images"),
              p(
                "If each composite contains your markers in distinct RGB ",
                "planes, upload the same complete batch of composite files ",
                "under every corresponding marker."
              ),
              tags$ol(
                tags$li(
                  "In Setup, assign the correct source plane to each marker. ",
                  "For example, choose R for PGL-1 and G for MUT-16 only if ",
                  "those markers are stored in the red and green planes."
                ),
                tags$li(
                  "In the condition card, click Browse under PGL-1 and select ",
                  "the complete batch of composite images."
                ),
                tags$li(
                  "Click Browse under MUT-16 and select the same complete ",
                  "batch again. Keep the filenames unchanged."
                )
              ),
              p(
                "For example, select ",
                tags$code("Composite_1.tif"), ", ",
                tags$code("Composite_2.tif"), " and ",
                tags$code("Composite_3.tif"),
                " together under PGL-1, then select those same three files ",
                "under MUT-16. You do not need to upload each nucleus separately."
              ),
              p(
                tags$strong("Important: "),
                "This works only when each marker occupies a distinct RGB ",
                "plane. A pseudocolored composite may mix marker signals ",
                "across planes. If the signals are mixed, use separate ",
                "grayscale marker images instead."
              ),
              
              h4("6. Check the pairing manifest"),
              p(
                "The pairing manifest lists the uploaded filenames, their ",
                "condition, nucleus ID, marker and pairing status."
              ),
              tags$ul(
                tags$li(
                  tags$strong("ready: "),
                  "Exactly one file is available for every configured marker ",
                  "with that nucleus ID."
                ),
                tags$li(
                  tags$strong("excluded: "),
                  "A marker file is missing, or more than one file has the ",
                  "same nucleus ID within a marker. Read the detail column ",
                  "and correct the upload."
                ),
                tags$li(
                  tags$strong("invalid: "),
                  "The filename does not end with a numeric nucleus ID ",
                  "before the extension. Rename the file and upload again."
                )
              ),
              p(
                tags$strong("A ready status confirms filename pairing only. "),
                "It does not verify alignment, image quality or the selected ",
                "RGB plane. Image reading and dimension checks occur during ",
                "analysis; inspect the images and segmentation in QC."
              ),
              p(
                "Click Download pairing manifest to save the current ",
                "pairing table as a CSV file."
              ),
              
              h4("7. Replace uploads or remove a condition"),
              p(
                tags$strong(
                  "Uploading again replaces the previous selection for ",
                  "that marker and condition; it does not append files."
                ),
                " To add more images, select the full set you want to keep, ",
                "including previously uploaded images and the new ones."
              ),
              p(
                "If you rename files to correct pairing, upload the corrected ",
                "full set under each affected marker and check the manifest again."
              ),
              p(
                "Click Remove condition to remove that condition and its ",
                "uploads from the analysis. This does not delete the original ",
                "files on your computer."
              ),
              
              h4("8. Continue to QC"),
              p(
                "Once your intended nuclei are listed as ready, open QC ",
                "to inspect the images and preview segmentation before extraction."
              ),
              p(
                "Changing uploads or conditions makes previous calibration ",
                "and extraction results stale. Recalibrate if required by ",
                "your QC settings, then repeat preview and extraction ",
                "for the updated dataset."
              )
            )
          ),
          
          downloadButton("dl_manifest", "Download pairing manifest")
        ),
        column(
          9,
          div(id = "conditions_container"),
          h4("Pairing manifest"),
          DTOutput("manifest_table")
        )
      )
    ),
    tabPanel(
      "QC", value = "qc",
      tags$details(
        style = "margin-top: 12px; margin-bottom: 16px;",
        tags$summary(
          style = paste(
            "cursor: pointer;",
            "color: #0f766e;",
            "font-weight: 600;"
          ),
          "QC instructions"
        ),
        div(
          style = "padding-top: 12px; line-height: 1.6;",
          
          h4("1. Sample"),
          p(
            "Select the uploaded sample to inspect. QC previews one sample ",
            "at a time; it does not replace extraction of the full dataset. ",
            "After choosing a sample, click Preview segmentation."
          ),
          
          h4("2. Use legacy batch-noise scaling"),
          p(
            "This checkbox enables the existing noise-based adjustment ",
            "for channels using DoG segmentation. It adjusts the local ",
            "contrast cutoff using the sample's estimated background noise ",
            "relative to the dataset reference. It does not normalize ",
            "raw image intensities and does not affect top-hat channels."
          ),
          p(
            "When enabled and a DoG channel is used, current calibration ",
            "is required for preview and extraction. When disabled, the ",
            "entered local contrast cutoff is used without batch scaling."
          ),
          
          h4("3. Calibrate current dataset"),
          p(
            "Estimates background-noise references for DoG channels using ",
            "the complete paired samples in the current dataset. Inspect ",
            "any warning or error notifications before continuing."
          ),
          p(
            "Changing channel settings, segmentation parameters or uploads ",
            "can make calibration stale. Recalibrate when required, then ",
            "repeat the preview. Calibration is not an intensity-standardization ",
            "procedure and does not guarantee comparable acquisition settings."
          ),
          
          h4("4. Preview segmentation"),
          p(
            "Runs segmentation and measurements for the selected sample ",
            "using the current settings. The images, context plot and ",
            "Preview features table show this preview result."
          ),
          p(
            "After changing a sample or an analysis setting, click this ",
            "button again. Inspect representative samples from each ",
            "condition before running full extraction."
          ),
          
          h4("5. Overlay"),
          p(
            "Selects the annotation shown in the third image panel and ",
            "the corresponding context plot. Changing the overlay changes ",
            "the display, not the segmentation."
          ),
          tags$ul(
            tags$li(
              strong("Peak positions: "),
              "shows the brightest-pixel position selected within each ",
              "segmented object. These positions are not necessarily ",
              "the geometric centers of the objects."
            ),
            tags$li(
              strong("Nearest-neighbor links: "),
              "connects object peak positions to their nearest neighbors ",
              "within the selected channel. These are peak-to-peak ",
              "distances, not edge-to-edge gaps between masks. ",
              "At least two objects are needed for a neighbor distance."
            ),
            tags$li(
              strong("Radial spokes: image center: "),
              "connects object peaks to the image center. This is not a ",
              "distance from a segmented nuclear boundary or nuclear center."
            ),
            tags$li(
              strong("B peaks near A objects: "),
              "tests B peak positions against the expanded A-object masks ",
              "used by the existing proximity method. Circles indicate ",
              "proximity; crosses indicate no proximity. The reported net ",
              "value includes the existing chance correction."
            ),
            tags$li(
              strong("A peaks near B objects: "),
              "tests A peaks against expanded B-object masks. Its retained ",
              "legacy reported value is chance-adjusted and additionally ",
              "normalized by the number of B objects; it is not directly ",
              "equivalent to the reverse-direction percentage."
            ),
            tags$li(
              strong("Segmented-object overlap: "),
              "uses the original segmented masks without expanding them. ",
              "White pixels identify overlaps from object pairs meeting ",
              "the minimum-shared-pixels threshold below. This is distinct ",
              "from peak-based proximity."
            )
          ),
          
          h4("6. Distance-overlay channel"),
          p(
            "Chooses the channel used for nearest-neighbor links and ",
            "radial spokes. It does not change which channels are segmented. ",
            "For a pair-based overlay, use Channel pair instead."
          ),
          
          h4("7. Channel pair: A / B"),
          p(
            "Chooses the two channels used for proximity and segmented ",
            "overlap displays. The first channel is A and the second is B. ",
            "Read directional results using this order, rather than assuming ",
            "A and B always represent particular colors or markers."
          ),
          
          h4("8. Minimum shared pixels per object pair"),
          numericInput(
            "segcoloc_min_px",
            "Minimum shared pixels per object pair",
            value = 3,
            min = 1,
            step = 1
          ),
          p(
            "An A/B object pair qualifies for segmented overlap only when ",
            "its two masks share at least this many pixels. Overlaps with ",
            "different partners are not pooled to reach the threshold."
          ),
          p(
            "A higher threshold makes overlap qualification more stringent. ",
            "The default of 3 pixels is a starting setting, not a validated ",
            "cutoff for every image resolution or experiment. Check whether ",
            "conclusions depend strongly on this choice."
          ),
          p(
            "Qualifying pair count counts every qualifying A/B pair. ",
            "Directional percentages count each source object only once, ",
            "even if it overlaps several partners. Read percentages alongside ",
            "the total object counts in both channels and qualifying pair count. ",
            "A directional percentage with no source objects is undefined (NA)."
          ),
          p(
            "Changing this threshold requires a new preview and extraction, ",
            "but does not require repeating noise calibration."
          ),
          
          h4("9. Channel parameter cards"),
          p(
            "Each card controls segmentation for its corresponding channel. ",
            "The segmentation method is selected in Setup. Method-specific ",
            "controls appear only for the relevant method."
          ),
          tags$ul(
            tags$li(
              strong("Threshold multiplier: "),
              "sets the candidate threshold on the filtered image response: ",
              "the mean positive response plus this multiplier times its ",
              "standard deviation. Increasing it is more stringent, but can ",
              "also change object boundaries, merging and splitting."
            ),
            tags$li(
              strong("Minimum area (pixels): "),
              "rejects candidate objects smaller than this pixel count. ",
              "Increasing it can remove small noise objects but also real ",
              "small granules."
            ),
            tags$li(
              strong("Maximum area (pixels): "),
              "rejects candidate objects larger than this pixel count. ",
              "It does not split large or merged objects into smaller ones. ",
              "Keep it at least as large as the minimum area."
            ),
            tags$li(
              strong("Top-hat brush diameter: "),
              "controls the disc-shaped background-estimation brush for ",
              "top-hat segmentation. It is a diameter in pixels, not an ",
              "object-area cutoff. The code caps the brush to fit the image ",
              "and uses an odd diameter."
            ),
            tags$li(
              strong("Fine sigma: "),
              "sets the smaller Gaussian smoothing scale for DoG, in pixels."
            ),
            tags$li(
              strong("Coarse sigma: "),
              "sets the larger Gaussian smoothing scale for DoG, in pixels. ",
              "DoG uses the difference between fine and coarse smoothed ",
              "images. The effective coarse scale must exceed the fine ",
              "scale; the code also limits it according to image size."
            ),
            tags$li(
              strong("Global floor multiplier: "),
              "DoG only. Requires an object's mean raw intensity to be at ",
              "least the image mean plus this multiplier times the image ",
              "intensity standard deviation. Higher values are stricter."
            ),
            tags$li(
              strong("Maximum eccentricity: "),
              "DoG only. Rejects objects more elongated than this limit. ",
              "Values lie between 0 and 1; lower limits favor rounder ",
              "objects and higher limits allow more elongated objects."
            ),
            tags$li(
              strong("Local contrast z minimum: "),
              "DoG only. Requires sufficient object intensity above the ",
              "local background, scaled by the local background's robust ",
              "variation estimate. Higher values are stricter. With batch ",
              "noise scaling enabled, the effective cutoff is adjusted ",
              "from this entered value."
            )
          ),
          p(
            "Local background is estimated from a ring around each candidate. ",
            "The retained method does not exclude neighboring objects from ",
            "that ring, so crowded regions can affect local contrast."
          ),
          
          h4("10. Raw image panel"),
          p(
            "Shows the input channels for visual comparison with the ",
            "segmentation. Display contrast is for viewing only; it does ",
            "not change the measured raw intensities. Do not infer ",
            "quantitative intensity differences from display brightness alone."
          ),
          
          h4("11. Segmentation image panel"),
          p(
            "Shows the retained segmented objects. Compare this with the ",
            "raw image to identify missed granules, background incorrectly ",
            "included as objects, merged neighbors or fragmented objects. ",
            "Overlap results depend on these masks being appropriate."
          ),
          
          h4("12. Annotated image panel"),
          p(
            "Shows the selected overlay. Use Distance-overlay channel for ",
            "within-channel distance annotations and Channel pair for ",
            "between-channel annotations. Read the caption together with ",
            "the selected overlay; different overlays report different metrics."
          ),
          
          h4("13. Context plot"),
          p(
            "Provides the supporting summary for the selected overlay below ",
            "the image panels. Its content and units depend on the overlay: ",
            "read its title, axis labels and legend rather than treating ",
            "every plotted value as a colocalization percentage. Missing or ",
            "insufficient objects can make some summaries unavailable."
          ),
          
          h4("14. Preview features table"),
          p(
            "Lists feature names and numeric values for the selected preview ",
            "sample, not a condition average or a completed dataset extraction. ",
            "Use the table's search, sorting and page controls to inspect ",
            "individual measurements."
          ),
          p(
            "NA means that a measurement is unavailable or undefined, ",
            "not zero. Consult Measurement guide for feature definitions. ",
            "The legacy _ring feature describes nearest-neighbor distance ",
            "variability, not ring detection; _radcv describes radial-distance ",
            "variability about the image center."
          ),
          
          h4("15. Recommended QC sequence"),
          tags$ol(
            tags$li("Confirm the channel setup and uploaded sample pairing."),
            tags$li("Choose a sample and review the channel parameters."),
            tags$li(
              "Calibrate if batch-noise scaling is enabled for DoG channels."
            ),
            tags$li("Click Preview segmentation."),
            tags$li(
              "Compare raw and segmented images, then inspect the overlays ",
              "and numerical features."
            ),
            tags$li(
              "If settings change, recalibrate when required and repeat ",
              "the preview. Check representative samples across conditions."
            ),
            tags$li(
              "Once settings are satisfactory, use Extract to analyze ",
              "the full dataset. QC preview alone does not do this."
            )
          )
        )
      ),
      fluidRow(
        column(
          3,
          selectInput("qc_sample", "Sample", choices = character()),
          checkboxInput(
            "calibration_enabled",
            "Use legacy batch-noise scaling", TRUE
          ),
          actionButton("calibrate", "Calibrate current dataset"),
          actionButton("qc_run", "Preview segmentation", class = "btn-primary"),
          selectInput("ann_mode", "Overlay", OVERLAYS, selected = "peaks"),
          selectInput("overlay_channel", "Distance-overlay channel", character()),
          selectInput("overlay_pair", "Channel pair: A / B", character()),
          lapply(seq_len(MAX_CHANNELS), parameter_card)
        ),
        column(
          9,
          fluidRow(
            column(4, plotOutput("qc_raw", height = "290px")),
            column(4, plotOutput("qc_seg", height = "290px")),
            column(4, plotOutput("qc_ann", height = "290px"))
          ),
          plotOutput("qc_context", height = "360px"),
          h4("Preview features"),
          DTOutput("qc_features")
        )
      )
    ),
    
    tabPanel(
      "Extract", value = "extract",
      actionButton("run_extract", "Run extraction", class = "btn-primary"),
      verbatimTextOutput("extraction_summary"),
      fluidRow(
        column(3, downloadButton("dl_features", "Per-image CSV")),
        column(3, downloadButton("dl_agg", "Condition summary CSV")),
        column(3, downloadButton("dl_objects", "Object measurements CSV")),
        column(3, downloadButton("dl_log", "Processing log CSV"))
      ),
      downloadButton("dl_session", "Analysis snapshot RDS"),
      helpText(
        "The RDS snapshot records results, labels, settings, calibration,",
        "manifest and session information. It is not a project importer",
        "and does not preserve original uploaded image files."
      ),
      h4("Per-image features"),
      DTOutput("features_table"),
      h4("Condition summaries: counts = mean; other features = median"),
      DTOutput("aggregate_table"),
      h4("Processing log"),
      DTOutput("log_table")
    ),
    
    tabPanel(
      "Compare", value = "compare",
      tags$details(
        style = "margin-top: 12px; margin-bottom: 16px;",
        tags$summary(
          style = paste(
            "cursor: pointer;",
            "color: #0f766e;",
            "font-weight: 600;"
          ),
          "Compare instructions"
        ),
        div(
          style = "padding-top: 12px; line-height: 1.7;",
          
          h4("1. What this section compares"),
          p(
            "Compares per-image feature measurements between two conditions. ",
            "Each observation is an image-level measurement, such as an ",
            "object count, median granule area or segmented-overlap percentage. ",
            "This section does not test individual granule rows or the ",
            "condition-summary table."
          ),
          p(
            "Run extraction first. To analyze individual granules separately, ",
            "download Object measurements CSV from the Extract tab."
          ),
          
          h4("2. Condition pair"),
          p(
            "Choose the two conditions to compare. The results identify them ",
            "as condition_A and condition_B. Differences are calculated as ",
            "B minus A: positive values indicate a higher value in B, and ",
            "negative values indicate a lower value in B."
          ),
          
          h4("3. Test"),
          tags$ul(
            tags$li(
              strong("Wilcoxon rank-sum: "),
              "a two-sided, rank-based comparison of two independent groups. ",
              "It does not require normally distributed measurements. ",
              "It is not generally a test of medians alone; differences in ",
              "distribution shape or spread can affect its interpretation."
            ),
            tags$li(
              strong("Welch t-test: "),
              "a two-sided comparison of group means that allows unequal ",
              "variances. Outliers and strongly non-normal measurements ",
              "can affect its reliability, particularly with small samples."
            )
          ),
          p(
            "Both tests assume independent observations. Multiple images ",
            "from the same animal or sample may be related; this section ",
            "does not account for pairing or hierarchical grouping. ",
            "This limitation does not prevent measurement extraction."
          ),
          
          h4("4. BH correction"),
          p(
            "When enabled, Benjamini-Hochberg correction adjusts p-values ",
            "across the features in the selected condition-pair comparison ",
            "to control the false discovery rate under its assumptions. ",
            "It does not correct across other condition pairs."
          ),
          p(
            "When disabled, adj_p is the unadjusted p-value and the ",
            "correction column reports none."
          ),
          
          h4("5. Alpha"),
          p(
            "Sets the significance threshold. A feature is marked ",
            "significant when its finite adj_p value is strictly below ",
            "alpha. With BH correction disabled, this threshold is applied ",
            "to the unadjusted p-value. Statistical significance does not ",
            "measure the size or biological importance of a difference."
          ),
          
          h4("6. Run comparison"),
          p(
            "Click Run comparison after choosing the condition pair, test, ",
            "correction and alpha. Rerun the comparison after changing ",
            "these settings or generating a new extraction."
          ),
          p(
            "For each feature, non-finite measurements are excluded. ",
            "The app requires at least three finite image-level observations ",
            "in each condition before attempting a test. Meeting this ",
            "minimum does not guarantee adequate statistical power."
          ),
          
          h4("7. Results table"),
          tags$ul(
            tags$li(
              "feature: the exported image-level measurement being compared."
            ),
            tags$li(
              "condition_A and condition_B: the selected comparison groups."
            ),
            tags$li(
              "test: wilcox for Wilcoxon rank-sum or ttest for Welch t-test."
            ),
            tags$li(
              "n_A and n_B: finite image-level observations used for that ",
              "feature, not total granule counts."
            ),
            tags$li(
              "mean_A, mean_B, median_A and median_B: summaries of the ",
              "available image-level measurements in each condition."
            ),
            tags$li(
              "delta_mean and delta_median: the corresponding difference ",
              "between B and A, in the feature's original units."
            ),
            tags$li(
              "p_value: the unadjusted test p-value. adj_p: the value used ",
              "for significance after the selected correction."
            ),
            tags$li(
              "correction, alpha and significant: the correction method, ",
              "threshold and resulting significance flag."
            ),
            tags$li(
              "error and warnings: diagnostic messages. An unavailable ",
              "p-value is not evidence that the conditions are equivalent."
            )
          ),
          p(
            "NA indicates an unavailable or undefined measurement, not zero. ",
            "Use the Measurement guide for feature definitions."
          ),
          
          h4("8. Difference plot"),
          p(
            "Shows median differences, B minus A, for features marked ",
            "significant that have a finite median difference. ",
            "The plot uses median differences even when Welch t-test ",
            "is selected; Welch t-test itself compares means."
          ),
          p(
            "Differences are not standardized. Features have different ",
            "units, so bar lengths should not be compared as a ranking ",
            "of biological importance. No significant features does not ",
            "demonstrate that the conditions are identical."
          ),
          
          h4("9. Download results"),
          p(
            "Use Numeric statistics CSV to export the comparison table. ",
            "Keep the per-image measurements alongside it. For segmented ",
            "colocalization, also examine total object counts and qualifying ",
            "pair counts when interpreting percentage differences."
          )
        )
      ),
      fluidRow(
        column(
          3,
          selectInput("compare_pair", "Condition pair", character()),
          radioButtons(
            "stat_test", "Test",
            c("Wilcoxon rank-sum" = "wilcox", "Welch t-test" = "ttest")
          ),
          checkboxInput("fdr_correct", "BH correction across this pair's features", TRUE),
          numericInput("alpha_thresh", "Alpha", 0.05, min = 0.001, max = 1),
          actionButton("run_stats", "Run comparison", class = "btn-primary"),
          downloadButton("dl_stats", "Numeric statistics CSV")
        ),
        column(
          9,
          DTOutput("stats_table"),
          plotOutput("effect_plot", height = "340px")
        )
      )
    ),
    
    tabPanel(
      "Visualize", value = "visualise",
      fluidRow(
        column(
          3,
          checkboxInput("show_sig_only", "Significant features only", FALSE),
          selectInput("viz_feature", "Feature", character()),
          radioButtons(
            "plot_type", "Distribution plot",
            c("Box + points" = "box", "Violin + points" = "violin",
              "Mean + SE" = "bar")
          ),
          selectInput("browse_sample", "Extracted image", character()),
          selectInput("ann_mode_browse", "Browser overlay", OVERLAYS),
          helpText("The browser uses stored extraction results, not a new segmentation."),
          tags$hr(),
          h4("Export"),
          selectInput("exp_format", "Format", EXPORT_FORMATS, selected = "pdf"),
          numericInput("exp_width_mm", "Width (mm)", 90, min = 20, max = 500, step = 1),
          numericInput("exp_height_mm", "Distribution plot height (mm)", 70,
                       min = 20, max = 500, step = 1),
          numericInput("exp_dpi", "PNG resolution (dpi)", 600,
                       min = 72, max = 1200, step = 50),
          checkboxInput("exp_clean",
                        "Raw / segmentation: clean image (no title)", TRUE),
          downloadButton("dl_viz_dist", "Distribution plot"),
          downloadButton("dl_viz_raw",  "Raw composite"),
          downloadButton("dl_viz_seg",  "Segmentation"),
          downloadButton("dl_viz_ann",  "Annotated overlay")
        ),
        column(
          9,
          plotOutput("distribution_plot", height = "340px"),
          fluidRow(
            column(4, plotOutput("browse_raw", height = "270px")),
            column(4, plotOutput("browse_seg", height = "270px")),
            column(4, plotOutput("browse_ann", height = "270px"))
          )
        )
      )
    ),
    
    tabPanel(
      "Measurement guide",
      h4("Interpret exported feature names carefully"),
      tags$ul(
        tags$li("_n: detected object count."),
        tags$li("_size: median segmented object area in pixels; NA if no objects."),
        tags$li("_elong: median EBImage eccentricity."),
        tags$li("_ring: nearest-neighbor distance CV; requires at least four objects."),
        tags$li("_radcv: distance CV about the image center; requires at least three objects."),
        tags$li("<B>_near_<A>_net: B-peaks-near-A-objects gross proximity minus legacy chance percentage, in percentage points. A and B are the marker names of the selected channel pair (formerly _green_app_net)."),
        tags$li("<A>_near_<B>_net_norm: A-peaks-near-B-objects net proximity divided by B object count (formerly _red_app_net_norm)."),
        tags$li("_sizeratio: median A object area divided by median B object area."),
        tags$li(
          "_segcoloc_pairs_n: number of A-B object pairs sharing at least ",
          "the configured minimum number of pixels; masks are not dilated."
        ),
        tags$li(
          "_segcoloc_b_to_a_pct: percentage of B objects with qualifying ",
          "overlap with at least one A object."
        ),
        tags$li(
          "_segcoloc_a_to_b_pct: percentage of A objects with qualifying ",
          "overlap with at least one B object."
        )
      ),
      tags$details(
        style = "margin-top: 12px; margin-bottom: 16px;",
        tags$summary(
          style = paste(
            "cursor: pointer;",
            "color: #0f766e;",
            "font-weight: 600;"
          ),
          "How segmented-object colocalization is calculated"
        ),
        div(
          style = "padding-top: 12px; line-height: 1.7;",
          
          p(
            "This MatCol-inspired metric counts shared pixels between ",
            "individual labeled objects in two channels. It uses the ",
            "existing segmentation masks without dilation or additional ",
            "intensity thresholding."
          ),
          
          p(
            "The minimum shared-pixel criterion applies to each A-B pair ",
            "separately. Subthreshold overlaps with different objects are ",
            "not added together to create a qualifying pair."
          ),
          
          p(
            "If one A object overlaps two B objects and both pairs meet ",
            "the threshold, the pair count is two. For directional ",
            "percentages, that A object is counted once and the two B ",
            "objects are each counted once."
          ),
          
          p(
            "If either channel has no objects, the pair count is zero. ",
            "A directional percentage is NA when its source channel has ",
            "no objects; it is zero when source objects exist but none ",
            "has qualifying overlap."
          ),
          
          p(
            "This is two-dimensional image overlap, not evidence of direct ",
            "molecular interaction. It is separate from the existing ",
            "proximity and chance-adjusted proximity measurements."
          ),
          
          p(
            "This implementation does not reproduce MatCol background ",
            "thresholding, Monte Carlo randomization or within-image ",
            "significance testing."
          ),
          
          p(
            "Reference: Khushi et al. (2017), MatCol: a tool to measure ",
            "fluorescence signal colocalisation in biological systems. ",
            "Scientific Reports 7, 8879."
          )
        )
      ),
      p(
        "Proximity uses a disc brush of diameter 7 pixels, capped to fit",
        "the image. The legacy chance mask is the union of channel pixels",
        "above mean + 0.5 SD. Adding channels can therefore change",
        "chance-adjusted measurements."
      ),
      p(
        "All spatial measurements are in pixels. Images require matching",
        "registration and appropriate acquisition settings. The app does",
        "not perform chromatic registration, deconvolution or nuclear segmentation."
      ),
      p(
        "Legacy batch-noise scaling adjusts a contrast threshold that",
        "already uses local MAD. Its scientific appropriateness must be",
        "validated for your data."
      )
    )
  )
)

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {
  rv <- reactiveValues(
    setup = NULL,
    conditions = list(),
    next_id = 0L,
    calibration = NULL,
    preview = NULL,
    run = NULL,
    stats = NULL
  )
  
  upload_root <- tempfile("ggprofiler_")
  dir.create(upload_root)
  session$onSessionEnded(function() {
    unlink(upload_root, recursive = TRUE)
  })
  
  notify_error <- function(e) {
    showNotification(conditionMessage(e), type = "error", duration = 12)
  }
  
  setup_draft <- reactive({
    n <- as.integer(input$n_channels)
    lapply(seq_len(n), function(i) {
      list(
        name = trimws(input[[paste0("ch_name_", i)]] %||% ""),
        src = as.integer(input[[paste0("ch_src_", i)]] %||% 1),
        color = unname(COLORS[input[[paste0("ch_color_", i)]] %||% "Cyan"]),
        method = input[[paste0("ch_method_", i)]] %||% "tophat"
      )
    })
  })
  
  observeEvent(input$apply_setup, {
    tryCatch({
      cfg <- setup_draft()
      names <- vapply(cfg, `[[`, character(1), "name")
      
      if (any(!nzchar(names))) stop("Enter every channel name.")
      if (anyDuplicated(names) || anyDuplicated(make.names(names))) {
        stop("Channel names must also be unique after R name sanitization.")
      }
      if (length(rv$conditions) && !identical(cfg, rv$setup)) {
        stop("Remove all conditions before changing channel setup.")
      }
      
      rv$setup <- cfg
      updateTabsetPanel(session, "main_tabs", selected = "crop")
    }, error = notify_error)
  })
  
  output$setup_status <- renderText({
    if (is.null(rv$setup)) return("Setup not yet applied.")
    if (!identical(setup_draft(), rv$setup)) {
      return("Unapplied setup changes. Restore or apply them before analysis.")
    }
    paste(
      "Applied:",
      paste(vapply(rv$setup, `[[`, character(1), "name"), collapse = ", ")
    )
  })
  
  configs <- reactive({
    req(rv$setup)
    lapply(seq_along(rv$setup), function(i) {
      cfg <- rv$setup[[i]]
      getnum <- function(prefix, default) {
        input[[paste0(prefix, i)]] %||% default
      }
      
      cfg$params <- list(
        thr_mult = getnum("p_thr_", 1),
        min_area = getnum("p_minarea_", 10),
        max_area = getnum("p_maxarea_", 3000),
        brush_size = getnum("p_brush_", 51),
        sigma1 = getnum("p_sigma1_", 1),
        sigma2 = getnum("p_sigma2_", 8),
        global_mult = getnum("p_globmult_", 1.3),
        max_eccentricity = getnum("p_ecc_", 0.90),
        contrast_z_min = getnum("p_zmin_", 1)
      )
      cfg
    })
  })
  
  validate_settings <- function(cfg) {
    if (!identical(setup_draft(), rv$setup)) {
      stop("Setup has unapplied changes.")
    }
    
    for (ch in cfg) {
      p <- ch$params
      if (any(!is.finite(unlist(p)))) stop("Parameters must be finite numbers.")
      if (p$min_area < 1 || p$max_area < p$min_area) {
        stop(ch$name, ": invalid area bounds.")
      }
      if (p$thr_mult < 0 || p$global_mult < 0 || p$contrast_z_min < 0) {
        stop(ch$name, ": thresholds cannot be negative.")
      }
      if (p$brush_size < 1) stop(ch$name, ": brush size must be positive.")
      if (p$max_eccentricity < 0 || p$max_eccentricity > 1) {
        stop(ch$name, ": eccentricity must be between 0 and 1.")
      }
      if (ch$method == "dog" &&
          (p$sigma1 <= 0 || p$sigma2 <= p$sigma1)) {
        stop(ch$name, ": coarse sigma must exceed positive fine sigma.")
      }
    }
    invisible(TRUE)
  }
  
  # Persistent inserted upload controls: no renderUI rebuild on upload.
  observeEvent(input$add_cond, {
    tryCatch({
      req(rv$setup)
      if (!identical(setup_draft(), rv$setup)) stop("Apply or restore setup first.")
      name <- trimws(input$new_cond_name %||% "")
      if (!nzchar(name)) stop("Enter a condition name.")
      if (length(rv$conditions) >= MAX_CONDS) stop("Maximum condition count reached.")
      existing <- vapply(rv$conditions, `[[`, character(1), "name")
      if (name %in% existing) stop("Condition name already exists.")
      
      rv$next_id <- rv$next_id + 1L
      id <- paste0("condition_", rv$next_id)
      cfg <- rv$setup
      
      cond <- list(
        id = id, name = name,
        uploads = vector("list", length(cfg))
      )
      conditions <- rv$conditions
      conditions[[id]] <- cond
      rv$conditions <- conditions
      
      insertUI(
        selector = "#conditions_container",
        where = "beforeEnd",
        ui = div(
          id = id, class = "well",
          h4(name),
          fluidRow(lapply(seq_along(cfg), function(i) {
            column(
              max(3L, 12L %/% length(cfg)),
              fileInput(
                paste0(id, "_upload_", i),
                cfg[[i]]$name, multiple = TRUE,
                accept = c(".tif", ".tiff", ".png", ".jpg", ".jpeg")
              )
            )
          })),
          actionButton(paste0(id, "_remove"), "Remove condition")
        )
      )
      
      local({
        cid <- id
        observeEvent(input[[paste0(cid, "_remove")]], {
          conditions <- rv$conditions
          conditions[[cid]] <- NULL
          rv$conditions <- conditions
          removeUI(selector = paste0("#", cid))
        }, ignoreInit = TRUE)
        
        for (i in seq_along(cfg)) {
          local({
            channel_index <- i
            upload_id <- paste0(cid, "_upload_", channel_index)
            
            observeEvent(input[[upload_id]], {
              tryCatch({
                if (is.null(rv$conditions[[cid]])) return()
                files <- input[[upload_id]]
                req(files, nrow(files) > 0)
                
                # Copy uploads so their lifetime is owned by this session.
                destination <- vapply(files$name, function(filename) {
                  extension <- tools::file_ext(filename)
                  tempfile(
                    pattern = "image_", tmpdir = upload_root,
                    fileext = if (nzchar(extension)) paste0(".", extension) else ""
                  )
                }, character(1))
                
                copied <- file.copy(files$datapath, destination, overwrite = FALSE)
                if (!all(copied)) {
                  unlink(destination)
                  stop("Could not preserve uploaded files.")
                }
                
                conditions <- rv$conditions
                previous <- conditions[[cid]]$uploads[[channel_index]]
                conditions[[cid]]$uploads[[channel_index]] <- data.frame(
                  name = files$name, path = destination,
                  stringsAsFactors = FALSE
                )
                rv$conditions <- conditions
                
                if (!is.null(previous)) unlink(previous$path)
              }, error = notify_error)
            }, ignoreNULL = TRUE)
          })
        }
      })
      
      updateTextInput(session, "new_cond_name", value = "")
    }, error = notify_error)
  })
  
  manifest <- reactive({
    req(rv$setup)
    build_manifest(rv$conditions, rv$setup)
  })
  
  output$manifest_table <- renderDT({
    DT::datatable(
      manifest()$table, rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  observe({
    req(rv$setup)
    samples <- manifest()$samples
    labels <- vapply(samples, function(s) {
      paste(s$condition, s$stem, sep = " / ")
    }, character(1))
    
    choices <- setNames(as.character(seq_along(samples)), labels)
    updateSelectInput(session, "qc_sample", choices = choices)
    
    cfg <- rv$setup
    updateSelectInput(
      session, "overlay_channel",
      choices = setNames(
        seq_along(cfg), vapply(cfg, `[[`, character(1), "name")
      )
    )
    
    if (length(cfg) >= 2L) {
      pairs <- combn(seq_along(cfg), 2, simplify = FALSE)
      updateSelectInput(
        session, "overlay_pair",
        choices = setNames(
          vapply(pairs, paste, character(1), collapse = "_"),
          vapply(pairs, function(ij) {
            paste(cfg[[ij[1]]]$name, cfg[[ij[2]]]$name, sep = " / ")
          }, character(1))
        )
      )
    } else {
      updateSelectInput(session, "overlay_pair", choices = character())
    }
  })
  
  # Exact settings/data snapshots prevent stale results being presented
  # as current, even if a reactive invalidation event is delayed.
  calibration_key <- reactive({
    list(configs = configs(), conditions = rv$conditions)
  })
  
  overlap_minimum <- reactive({
    validate_overlap_min(input$segcoloc_min_px)
  })
  analysis_key <- reactive({
    list(
      calibration_key = calibration_key(),
      calibration_enabled = isTRUE(input$calibration_enabled),
      draft = setup_draft(),
      segcoloc_min_px = input$segcoloc_min_px
    )
  })
  
  calibration_current <- reactive({
    cal <- rv$calibration
    if (is.null(cal) || !identical(cal$key, calibration_key())) return(NULL)
    cal
  })
  
  run_current <- reactive({
    run <- rv$run
    if (is.null(run) || !identical(run$key, analysis_key())) return(NULL)
    run
  })
  
  preview_current <- reactive({
    preview <- rv$preview
    if (is.null(preview) || !identical(preview$key, analysis_key())) return(NULL)
    if (!identical(preview$selection, input$qc_sample)) return(NULL)
    preview$result
  })
  
  selected_pair <- reactive({
    value <- input$overlay_pair %||% "1_2"
    result <- suppressWarnings(as.integer(strsplit(value, "_", fixed = TRUE)[[1]]))
    if (length(result) != 2L || anyNA(result)) c(1L, 2L) else result
  })
  
  selected_channel <- reactive({
    x <- suppressWarnings(as.integer(input$overlay_channel %||% 1))
    if (length(x) != 1L || is.na(x)) 1L else x
  })
  
  observeEvent(input$calibrate, {
    tryCatch({
      cfg <- configs()
      validate_settings(cfg)
      samples <- manifest()$samples
      if (!length(samples)) stop("No complete paired samples.")
      dog <- which(vapply(cfg, function(ch) ch$method == "dog", logical(1)))
      reference <- rep(NA_real_, length(cfg))
      logs <- list()
      
      withProgress(message = "Calibrating noise", value = 0, {
        total <- max(length(dog) * length(samples), 1L)
        done <- 0L
        
        for (i in dog) {
          values <- numeric(0)
          for (s in samples) {
            captured <- capture_analysis({
              m <- read_channel(s$paths[i], cfg[[i]]$src)
              noise_scale(m, cfg[[i]]$params)
            })
            if (is.null(captured$error)) {
              values <- c(values, captured$result)
            }
            
            logs[[length(logs) + 1L]] <- data.frame(
              condition = s$condition, sample = s$stem,
              channel = cfg[[i]]$name,
              error = captured$error %||% "",
              warnings = captured$warnings
            )
            
            done <- done + 1L
            setProgress(done / total)
          }
          reference[i] <- finite_median(values)
        }
      })
      
      log <- bind_rows(logs)
      rv$calibration <- list(
        key = calibration_key(), reference = reference,
        log = log, timestamp = Sys.time()
      )
      
      if (length(dog) && any(!is.finite(reference[dog]))) {
        stop("Calibration failed for at least one DoG channel; see status.")
      }
      
      if (nrow(log) && any(nzchar(log$error) | nzchar(log$warnings))) {
        showNotification(
          "Calibration completed with errors/warnings. Inspect calibration status.",
          type = "warning", duration = 12
        )
      } else {
        showNotification("Calibration complete.", type = "message")
      }
    }, error = notify_error)
  })
  
  output$calibration_status <- renderPrint({
    if (!isTRUE(input$calibration_enabled)) {
      cat("Legacy batch scaling disabled for both preview and extraction.\n")
      return()
    }
    
    cal <- calibration_current()
    if (is.null(cal)) {
      cat("Calibration absent or stale. Click Calibrate current dataset.\n")
      return()
    }
    
    print(data.frame(
      channel = vapply(configs(), `[[`, character(1), "name"),
      reference_MAD = cal$reference
    ))
    
    if (nrow(cal$log)) {
      issues <- cal$log[nzchar(cal$log$error) | nzchar(cal$log$warnings), ]
      if (nrow(issues)) print(issues, row.names = FALSE)
    }
  })
  
  get_reference <- function(cfg) {
    if (!isTRUE(input$calibration_enabled)) return(NULL)
    dog <- which(vapply(cfg, function(ch) ch$method == "dog", logical(1)))
    if (!length(dog)) return(NULL)
    
    cal <- calibration_current()
    if (is.null(cal)) stop("Calibrate the current dataset/settings first.")
    if (any(!is.finite(cal$reference[dog]) | cal$reference[dog] <= 0)) {
      stop("Calibration is invalid for a DoG channel.")
    }
    cal$reference
  }
  
  observeEvent(input$qc_run, {
    rv$preview <- NULL
    tryCatch({
      cfg <- configs()
      validate_settings(cfg)
      reference <- get_reference(cfg)
      samples <- manifest()$samples
      index <- suppressWarnings(as.integer(input$qc_sample))
      if (length(index) != 1L || is.na(index) ||
          index < 1L || index > length(samples)) stop("Select a valid sample.")
      
      captured <- NULL
      withProgress(message = "Previewing segmentation", value = 0.5, {
        captured <- capture_analysis(
          extract_sample(
            samples[[index]],
            cfg,
            reference,
            min_shared_px = overlap_minimum()
          )
        )
      })
      
      if (!is.null(captured$error)) stop(captured$error)
      rv$preview <- list(
        key = analysis_key(), selection = input$qc_sample,
        result = captured$result
      )
      
      if (nzchar(captured$warnings)) {
        showNotification(captured$warnings, type = "warning", duration = 12)
      }
    }, error = notify_error)
  })
  
  output$qc_raw <- renderPlot(draw_image(preview_current(), "raw"))
  output$qc_seg <- renderPlot(draw_image(preview_current(), "seg"))
  output$qc_ann <- renderPlot({
    draw_image(
      preview_current(), input$ann_mode,
      selected_pair(), selected_channel()
    )
  })
  output$qc_context <- renderPlot({
    draw_context(
      preview_current(), input$ann_mode,
      selected_pair(), selected_channel()
    )
  })
  
  output$qc_features <- renderDT({
    res <- preview_current()
    req(res)
    df <- data.frame(
      feature = names(res$feats),
      value = as.numeric(res$feats[1, ])
    )
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 10)) |>
      DT::formatSignif("value", digits = 5)
  })
  
  output$effective_settings <- renderPrint({
    req(preview_current())
    print(lapply(preview_current()$effective, function(ch) {
      list(channel = ch$name, method = ch$method, parameters = ch$params)
    }))
  })
  
  observeEvent(input$run_extract, {
    tryCatch({
      cfg <- configs()
      validate_settings(cfg)
      samples <- manifest()$samples
      if (!length(samples)) stop("No complete paired samples.")
      reference <- get_reference(cfg)
      key <- analysis_key()
      
      rows <- list()
      objects <- list()
      logs <- list()
      results <- list()
      result_labels <- character(0)
      
      withProgress(message = "Extracting features", value = 0, {
        for (i in seq_along(samples)) {
          s <- samples[[i]]
          setProgress(
            i / length(samples),
            detail = paste(s$condition, s$stem, sep = " / ")
          )
          
          captured <- capture_analysis(
            extract_sample(
              s,
              cfg,
              reference,
              min_shared_px = overlap_minimum()
            )
          )
          res <- captured$result
          
          status <- if (!is.null(captured$error)) {
            "failed"
          } else if (all(vapply(res$channels, `[[`, numeric(1), "n") == 0)) {
            "ok_zero_objects"
          } else {
            "ok"
          }
          
          logs[[i]] <- data.frame(
            condition = s$condition, image_name = s$stem,
            status = status, error = captured$error %||% "",
            warnings = captured$warnings,
            stringsAsFactors = FALSE
          )
          
          if (!is.null(captured$error)) next
          
          rows[[length(rows) + 1L]] <- cbind(
            data.frame(
              condition = s$condition, image_name = s$stem,
              stringsAsFactors = FALSE
            ),
            res$feats
          )
          
          if (nrow(res$objects)) {
            objects[[length(objects) + 1L]] <- cbind(
              data.frame(
                condition = s$condition, image_name = s$stem,
                stringsAsFactors = FALSE
              ),
              res$objects
            )
          }
          
          results[[length(results) + 1L]] <- res
          result_labels <- c(
            result_labels, paste(s$condition, s$stem, sep = " / ")
          )
        }
      })
      
      df <- bind_rows(rows)
      feature_cols <- setdiff(names(df), c("condition", "image_name"))
      
      aggregated <- if (nrow(df)) {
        bind_rows(lapply(split(df, df$condition), function(part) {
          values <- lapply(feature_cols, function(feature) {
            if (grepl("_n$", feature)) {
              finite_mean(part[[feature]])
            } else {
              finite_median(part[[feature]])
            }
          })
          names(values) <- feature_cols
          cbind(
            data.frame(condition = part$condition[1], n_images = nrow(part)),
            as.data.frame(values, check.names = FALSE)
          )
        }))
      } else data.frame()
      
      rv$run <- list(
        key = key,
        id = paste(format(Sys.time(), "%Y%m%d%H%M%OS6"), length(samples)),
        features = df,
        feature_cols = feature_cols,
        aggregated = aggregated,
        objects = bind_rows(objects),
        log = bind_rows(logs),
        results = results,
        result_labels = result_labels,
        manifest = manifest()$table,
        configs = cfg,
        calibration = calibration_current(),
        timestamp = Sys.time()
      )
      rv$stats <- NULL
      
      failed <- sum(rv$run$log$status == "failed")
      showNotification(
        sprintf("%d successful; %d failed.", nrow(df), failed),
        type = if (failed) "warning" else "message",
        duration = 10
      )
    }, error = notify_error)
  })
  
  output$extraction_summary <- renderText({
    run <- run_current()
    if (is.null(run)) {
      return("No current extraction. Previous results, if any, are stale.")
    }
    paste0(
      "Processed: ", nrow(run$log),
      "\nSuccessful: ", nrow(run$features),
      "\nFailed: ", sum(run$log$status == "failed"),
      "\nZero-object images: ", sum(run$log$status == "ok_zero_objects")
    )
  })
  
  output$features_table <- renderDT({
    req(run_current())
    DT::datatable(
      run_current()$features, rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  output$aggregate_table <- renderDT({
    req(run_current())
    DT::datatable(
      run_current()$aggregated, rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  output$log_table <- renderDT({
    req(run_current())
    DT::datatable(
      run_current()$log, rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  observe({
    run <- run_current()
    if (is.null(run)) {
      updateSelectInput(session, "compare_pair", choices = character())
      updateSelectInput(session, "browse_sample", choices = character())
      return()
    }
    
    updateSelectInput(
      session, "browse_sample",
      choices = setNames(seq_along(run$results), run$result_labels)
    )
    
    conditions <- unique(run$features$condition)
    if (length(conditions) >= 2L) {
      pairs <- combn(conditions, 2, simplify = FALSE)
      updateSelectInput(
        session, "compare_pair",
        choices = setNames(
          seq_along(pairs),
          vapply(pairs, paste, character(1), collapse = " vs ")
        )
      )
    } else {
      updateSelectInput(session, "compare_pair", choices = character())
    }
  })
  
  stats_key <- reactive({
    run <- run_current()
    if (is.null(run)) return(NULL)
    list(
      run_id = run$id,
      pair = input$compare_pair,
      test = input$stat_test,
      fdr = input$fdr_correct,
      alpha = input$alpha_thresh
    )
  })
  
  stats_current <- reactive({
    stats <- rv$stats
    if (is.null(stats) || !identical(stats$key, stats_key())) return(NULL)
    stats$table
  })
  
  observeEvent(input$run_stats, {
    tryCatch({
      run <- run_current()
      if (is.null(run) || !nrow(run$features)) stop("Run extraction first.")
      
      alpha <- input$alpha_thresh
      if (!is.finite(alpha) || alpha <= 0 || alpha > 1) stop("Invalid alpha.")
      
      conditions <- unique(run$features$condition)
      if (length(conditions) < 2L) stop("Two conditions with data are required.")
      pairs <- combn(conditions, 2, simplify = FALSE)
      index <- suppressWarnings(as.integer(input$compare_pair))
      
      if (length(index) != 1L || is.na(index) ||
          index < 1L || index > length(pairs)) stop("Select a valid comparison.")
      
      pair <- pairs[[index]]
      df <- run$features
      
      rows <- lapply(run$feature_cols, function(feature) {
        a <- df[df$condition == pair[1], feature]
        b <- df[df$condition == pair[2], feature]
        a <- a[is.finite(a)]
        b <- b[is.finite(b)]
        
        error <- ""
        warnings <- ""
        pvalue <- NA_real_
        
        if (length(a) < 3L || length(b) < 3L) {
          error <- "At least three finite observations per condition required."
        } else {
          captured <- capture_analysis({
            if (input$stat_test == "wilcox") {
              wilcox.test(a, b, exact = FALSE)$p.value
            } else {
              t.test(a, b, var.equal = FALSE)$p.value
            }
          })
          pvalue <- captured$result %||% NA_real_
          error <- captured$error %||% ""
          warnings <- captured$warnings
          if (!is.finite(pvalue) && !nzchar(error)) {
            error <- "Test returned an undefined p-value, possibly constant data."
          }
        }
        
        data.frame(
          feature = feature,
          condition_A = pair[1], condition_B = pair[2],
          test = input$stat_test,
          n_A = length(a), n_B = length(b),
          mean_A = finite_mean(a), mean_B = finite_mean(b),
          median_A = finite_median(a), median_B = finite_median(b),
          delta_mean = finite_mean(b) - finite_mean(a),
          delta_median = finite_median(b) - finite_median(a),
          p_value = pvalue,
          error = error, warnings = warnings,
          stringsAsFactors = FALSE
        )
      })
      
      table <- bind_rows(rows)
      table$adj_p <- if (isTRUE(input$fdr_correct)) {
        p.adjust(table$p_value, method = "BH")
      } else table$p_value
      table$correction <- if (isTRUE(input$fdr_correct)) "BH" else "none"
      table$alpha <- alpha
      table$significant <- is.finite(table$adj_p) & table$adj_p < alpha
      
      rv$stats <- list(key = stats_key(), table = table)
    }, error = notify_error)
  })
  
  output$stats_table <- renderDT({
    table <- stats_current()
    req(table)
    DT::datatable(
      table, rownames = FALSE,
      options = list(pageLength = 15, scrollX = TRUE)
    ) |>
      DT::formatSignif(c("p_value", "adj_p"), digits = 5)
  })
  
  output$effect_plot <- renderPlot({
    table <- stats_current()
    if (is.null(table)) {
      empty_plot("Run a current comparison.")
      return()
    }
    table <- table[table$significant & is.finite(table$delta_median), ]
    if (!nrow(table)) {
      empty_plot("No significant features.")
      return()
    }
    
    par(mar = c(4, 15, 3, 1))
    barplot(
      table$delta_median, names.arg = table$feature,
      horiz = TRUE, las = 1, cex.names = 0.7,
      col = ifelse(table$delta_median >= 0, "#009E73", "#D55E00"),
      xlab = "Median difference: B - A (feature-specific units)",
      main = "Significant features: unstandardized differences"
    )
    abline(v = 0)
  })
  
  observe({
    run <- run_current()
    table <- stats_current()
    features <- if (is.null(run)) character() else run$feature_cols
    
    if (isTRUE(input$show_sig_only)) {
      features <- if (is.null(table)) {
        character()
      } else table$feature[table$significant]
    }
    
    selected <- isolate(input$viz_feature)
    if (is.null(selected) || !selected %in% features) {
      selected <- if (length(features)) features[1] else character()
    }
    
    updateSelectInput(
      session, "viz_feature", choices = features, selected = selected
    )
  })
  
  # One plot object feeds both the screen and the export.
  distribution_obj <- function(pub = FALSE) {
    build_distribution_plot(
      run_current(), input$viz_feature, input$plot_type,
      stats_current(), pub = pub
    )
  }
  
  output$distribution_plot <- renderPlot({
    p <- distribution_obj(pub = FALSE)
    if (is.character(p)) {
      empty_plot(p)
      return(invisible(NULL))
    }
    print(p)
  })
  
  # ---- Visualize exports ----------------------------------------------------
  export_dims <- function() {
    w <- input$exp_width_mm
    h <- input$exp_height_mm
    d <- input$exp_dpi
    if (length(w) != 1L || !is.finite(w) || w < 20 || w > 500)
      stop("Width must be 20-500 mm.")
    if (length(h) != 1L || !is.finite(h) || h < 20 || h > 500)
      stop("Height must be 20-500 mm.")
    if (length(d) != 1L || !is.finite(d) || d < 72 || d > 1200)
      stop("Resolution must be 72-1200 dpi.")
    list(w = w, h = h, dpi = d)
  }
  
  export_stem <- function(suffix) {
    run <- run_current()
    i <- suppressWarnings(as.integer(input$browse_sample))
    lab <- if (!is.null(run) && length(i) == 1L && !is.na(i) &&
               i >= 1L && i <= length(run$result_labels)) {
      run$result_labels[i]
    } else "image"
    paste0("GGProfiler_", gsub("[^A-Za-z0-9_-]+", "_", lab), "_", suffix)
  }
  
  output$dl_viz_dist <- downloadHandler(
    filename = function() {
      paste0("GGProfiler_",
             gsub("[^A-Za-z0-9_-]+", "_", input$viz_feature %||% "feature"),
             "_", input$plot_type, ".", input$exp_format)
    },
    content = function(file) {
      p <- distribution_obj(pub = TRUE)
      if (is.character(p)) stop(p)
      d <- export_dims()
      export_distribution(file, p, input$exp_format, d$w, d$h, d$dpi)
    }
  )
  
  make_overlay_handler <- function(suffix, mode_fun) {
    downloadHandler(
      filename = function() paste0(export_stem(suffix()), ".", input$exp_format),
      content = function(file) {
        res <- browse_result()
        if (is.null(res)) stop("Run extraction and select an extracted image first.")
        d <- export_dims()
        export_overlay(
          file, res, mode_fun(), selected_pair(), selected_channel(),
          input$exp_format, d$w, d$dpi, clean = isTRUE(input$exp_clean)
        )
      }
    )
  }
  
  output$dl_viz_raw <- make_overlay_handler(
    function() if (isTRUE(input$exp_clean)) "raw_clean" else "raw",
    function() "raw"
  )
  output$dl_viz_seg <- make_overlay_handler(
    function() if (isTRUE(input$exp_clean)) "segmentation_clean" else "segmentation",
    function() "seg"
  )
  output$dl_viz_ann <- make_overlay_handler(
    function() paste0("overlay_", input$ann_mode_browse),
    function() input$ann_mode_browse
  )
  browse_result <- reactive({
    run <- run_current()
    if (is.null(run)) return(NULL)
    i <- suppressWarnings(as.integer(input$browse_sample))
    if (length(i) != 1L || is.na(i) ||
        i < 1L || i > length(run$results)) return(NULL)
    run$results[[i]]
  })
  
  output$browse_raw <- renderPlot(draw_image(browse_result(), "raw"))
  output$browse_seg <- renderPlot(draw_image(browse_result(), "seg"))
  output$browse_ann <- renderPlot({
    draw_image(
      browse_result(), input$ann_mode_browse,
      selected_pair(), selected_channel()
    )
  })
  # ---- Crop nuclei module ----------------------------------------------------
  # Display-only drawing; cropping uses raw pixel values. Coordinates follow
  # the EBImage convention: x = column (1..width), y = row (1..height, top = 1).
  
  crop_state <- reactiveValues(
    images = list(),
    rect = NULL,
    saved = list(),
    # Zoom window in plot coordinates: list(x0, x1, y0, y1), NULL = full image.
    view = NULL,
    message = "Upload the images of one field to begin."
  )
  
  crop_dims <- function(im) {
    d <- dim(EBImage::imageData(im))
    c(w = as.integer(d[1]), h = as.integer(d[2]))
  }
  
  crop_stretch <- function(m) {
    lims <- as.numeric(stats::quantile(m, c(0.001, 0.999), names = FALSE))
    if (!all(is.finite(lims)) || lims[2] <= lims[1]) lims <- range(m)
    if (!all(is.finite(lims)) || lims[2] <= lims[1]) {
      return(matrix(0.5, nrow(m), ncol(m)))
    }
    pmin(pmax((m - lims[1]) / (lims[2] - lims[1]), 0), 1)
  }
  
  # Returns a [row = y, column = x] display array for rasters.
  crop_display_array <- function(im, plane = 1L) {
    a <- EBImage::imageData(im)
    d <- dim(a)
    
    if (length(d) == 2L) {
      return(crop_stretch(t(matrix(as.numeric(a), d[1], d[2]))))
    }
    
    if (EBImage::colorMode(im) == EBImage::Color && d[3] >= 3L) {
      out <- array(0, dim = c(d[2], d[1], 3L))
      for (k in 1:3) {
        out[, , k] <- crop_stretch(
          t(matrix(as.numeric(a[, , k]), d[1], d[2]))
        )
      }
      return(out)
    }
    
    k <- min(max(as.integer(plane), 1L), d[3])
    crop_stretch(t(matrix(as.numeric(a[, , k]), d[1], d[2])))
  }
  
  crop_boxes <- function(lst, H) {
    if (!length(lst)) return(NULL)
    do.call(rbind, lapply(lst, function(r) {
      data.frame(
        xmin = r$x0 - 0.5,
        xmax = r$x1 + 0.5,
        ymin = H - r$y1 + 0.5,
        ymax = H - r$y0 + 1.5,
        label = as.character(r$id %||% ""),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  crop_plot <- function(im, plane = 1L, rect = NULL, saved = list(),
                        view = NULL) {
    d <- crop_dims(im)
    W <- d[["w"]]
    H <- d[["h"]]
    
    p <- ggplot2::ggplot() +
      ggplot2::annotation_raster(
        grDevices::as.raster(crop_display_array(im, plane)),
        xmin = 0.5, xmax = W + 0.5,
        ymin = 0.5, ymax = H + 0.5,
        interpolate = FALSE
      )
    
    same_size <- Filter(function(r) r$w == W && r$h == H, saved)
    saved_boxes <- crop_boxes(same_size, H)
    if (!is.null(saved_boxes)) {
      p <- p +
        ggplot2::geom_rect(
          data = saved_boxes,
          ggplot2::aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
          ),
          fill = NA, colour = "#00BCD4", linewidth = 0.6,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = saved_boxes,
          ggplot2::aes(x = xmin, y = ymax, label = label),
          hjust = -0.2, vjust = 1.3, colour = "#00BCD4", size = 4,
          inherit.aes = FALSE
        )
    }
    
    if (!is.null(rect) && rect$w == W && rect$h == H) {
      current <- crop_boxes(list(rect), H)
      p <- p +
        ggplot2::geom_rect(
          data = current,
          ggplot2::aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
          ),
          fill = NA, colour = "#D55E00", linetype = "dashed",
          linewidth = 0.8, inherit.aes = FALSE
        )
    }
    
    xlim <- c(0.5, W + 0.5)
    ylim <- c(0.5, H + 0.5)
    if (!is.null(view)) {
      xlim <- c(max(0.5, view$x0), min(W + 0.5, view$x1))
      ylim <- c(max(0.5, view$y0), min(H + 0.5, view$y1))
    }
    
    # Zooming only changes the visible window; scale limits stay at the full
    # image, so brush coordinates remain in full-image plot coordinates.
    p +
      ggplot2::scale_x_continuous(limits = c(0.5, W + 0.5), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.5, H + 0.5), expand = c(0, 0)) +
      ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::theme_void()
  }
  
  # Crop every plane/channel of one image to the rectangle.
  # If `plane` is given, a multi-plane (non-RGB) image is reduced to that
  # single plane (clamped to the plane count). RGB and 2D images are unchanged.
  # Pixel values are never rescaled.
  crop_subimage <- function(im, r, plane = NULL) {
    a <- EBImage::imageData(im)
    d <- dim(a)
    if (length(d) == 2L) {
      sub <- a[r$x0:r$x1, r$y0:r$y1, drop = FALSE]
      return(EBImage::Image(sub, colormode = EBImage::colorMode(im)))
    }
    
    is_rgb <- EBImage::colorMode(im) == EBImage::Color && d[3] >= 3L
    if (!is.null(plane) && !is_rgb) {
      k <- min(max(as.integer(plane), 1L), d[3])
      sub <- matrix(
        as.numeric(a[r$x0:r$x1, r$y0:r$y1, k]),
        nrow = r$x1 - r$x0 + 1L, ncol = r$y1 - r$y0 + 1L
      )
      return(EBImage::Image(sub, colormode = EBImage::Grayscale))
    }
    
    sub <- a[r$x0:r$x1, r$y0:r$y1, , drop = FALSE]
    EBImage::Image(sub, colormode = EBImage::colorMode(im))
  }
  
  # Plane chosen in the "Plane to display" selector (1 if none shown).
  crop_current_plane <- function() {
    p <- suppressWarnings(as.integer(input$crop_plane %||% 1L))
    if (length(p) != 1L || is.na(p) || p < 1L) 1L else p
  }
  
  crop_ref_name <- reactive({
    imgs <- crop_state$images
    req(length(imgs) > 0L)
    nm <- input$crop_ref
    if (is.null(nm) || !nzchar(nm) || !nm %in% names(imgs)) {
      nm <- names(imgs)[1]
    }
    nm
  })
  
  crop_ref_image <- reactive({
    crop_state$images[[crop_ref_name()]]
  })
  
  observeEvent(input$crop_files, {
    files <- input$crop_files
    req(files, nrow(files) > 0)
    
    names_unique <- make.unique(files$name)
    imgs <- list()
    problems <- character(0)
    
    for (i in seq_len(nrow(files))) {
      im <- tryCatch(
        EBImage::readImage(files$datapath[i]),
        error = function(e) e
      )
      if (inherits(im, "error")) {
        problems <- c(
          problems,
          paste0(files$name[i], ": ", conditionMessage(im))
        )
        next
      }
      
      d <- dim(EBImage::imageData(im))
      if (!length(d) %in% c(2L, 3L) || min(d[1:2]) < 7L) {
        problems <- c(
          problems,
          paste0(files$name[i], ": unsupported size ", paste(d, collapse = " x "))
        )
        next
      }
      imgs[[names_unique[i]]] <- im
    }
    
    crop_state$images <- imgs
    crop_state$rect <- NULL
    crop_state$saved <- list()
    updateNumericInput(session, "crop_id", value = 1)
    updateSelectInput(
      session, "crop_ref",
      choices = names(imgs),
      selected = if (length(imgs)) names(imgs)[1] else character()
    )
    
    dims <- if (length(imgs)) {
      vapply(imgs, function(im) {
        d <- crop_dims(im)
        paste(d[["w"]], "x", d[["h"]])
      }, character(1))
    } else character(0)
    
    crop_state$message <- paste(
      c(
        sprintf("Loaded %d image(s).", length(imgs)),
        if (length(imgs)) paste0("  ", names(imgs), ": ", dims),
        if (length(unique(dims)) > 1L) {
          "WARNING: images have different dimensions; only images matching the reference are cropped."
        },
        if (length(problems)) c("Could not load:", paste0("  ", problems))
      ),
      collapse = "\n"
    )
  })
  
  output$crop_plane_ui <- renderUI({
    im <- crop_ref_image()
    d <- dim(EBImage::imageData(im))
    if (length(d) == 3L && EBImage::colorMode(im) != EBImage::Color) {
      sliderInput(
        "crop_plane",
        "Z plane to display (crop applies to all planes)",
        min = 1, max = d[3], value = 1, step = 1, ticks = TRUE
      )
    } else {
      NULL
    }
  })
  
  output$crop_status <- renderText(crop_state$message)
  
  output$crop_img <- renderPlot({
    if (!length(crop_state$images)) {
      empty_plot("Upload field images to begin cropping.")
      return(invisible(NULL))
    }
    plane <- suppressWarnings(as.integer(input$crop_plane %||% 1L))
    if (length(plane) != 1L || is.na(plane)) plane <- 1L
    im <- crop_ref_image()
    d <- crop_dims(im)
    W <- d[["w"]]
    H <- d[["h"]]
    
    # Base graphics: Shiny maps brush coordinates to the plot's user
    # coordinates (here: image x = column, y flipped so that plot y = H - row + 1).
    xlim <- c(0.5, W + 0.5)
    ylim <- c(0.5, H + 0.5)
    v <- crop_state$view
    if (!is.null(v)) {
      xlim <- c(max(0.5, v$x0), min(W + 0.5, v$x1))
      ylim <- c(max(0.5, v$y0), min(H + 0.5, v$y1))
    }
    
    par(mar = c(0, 0, 0, 0), bg = "black")
    plot(
      NA, xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
      xlab = "", ylab = "", asp = 1, xaxs = "i", yaxs = "i"
    )
    rasterImage(
      grDevices::as.raster(crop_display_array(im, plane)),
      0.5, 0.5, W + 0.5, H + 0.5, interpolate = FALSE
    )
    
    same_size <- Filter(function(r) r$w == W && r$h == H, crop_state$saved)
    sb <- crop_boxes(same_size, H)
    if (!is.null(sb)) {
      rect(sb$xmin, sb$ymin, sb$xmax, sb$ymax, border = "#00BCD4", lwd = 2)
      text(sb$xmin, sb$ymax, labels = sb$label, col = "#00BCD4",
           adj = c(-0.3, 1.3), cex = 1.2)
    }
    
    r <- crop_state$rect
    if (!is.null(r) && r$w == W && r$h == H) {
      cb <- crop_boxes(list(r), H)
      rect(cb$xmin, cb$ymin, cb$xmax, cb$ymax,
           border = "#D55E00", lty = 2, lwd = 2)
    }
  })
  
  # ---- Zoom controls ----------------------------------------------------
  # Drag a box, click "Zoom to dragged box", then draw the crop rectangle
  # inside the zoomed view for precise selection.
  observeEvent(input$crop_zoom, {
    tryCatch({
      if (!length(crop_state$images)) stop("Upload field images first.")
      br <- input$crop_brush
      if (is.null(br)) stop("Drag a box over the area to zoom into first.")
      
      d <- crop_dims(crop_ref_image())
      x0 <- max(0.5, br$xmin)
      x1 <- min(d[["w"]] + 0.5, br$xmax)
      y0 <- max(0.5, br$ymin)
      y1 <- min(d[["h"]] + 0.5, br$ymax)
      
      if (x1 - x0 < 7 || y1 - y0 < 7) {
        stop("Zoom area is too small; drag a larger box (at least 7 x 7 pixels).")
      }
      
      crop_state$view <- list(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
      session$resetBrush("crop_brush")
      crop_state$message <- sprintf(
        paste0(
          "Zoomed to x %d-%d, y %d-%d. Drag the crop rectangle inside this ",
          "view, then click Apply rectangle."
        ),
        as.integer(round(x0)), as.integer(round(x1)),
        as.integer(round(d[["h"]] - y1 + 1)), as.integer(round(d[["h"]] - y0 + 1))
      )
    }, error = notify_error)
  })
  
  observeEvent(input$crop_zoom_out, {
    v <- crop_state$view
    if (is.null(v) || !length(crop_state$images)) return()
    d <- crop_dims(crop_ref_image())
    cx <- (v$x0 + v$x1) / 2
    cy <- (v$y0 + v$y1) / 2
    hw <- (v$x1 - v$x0)  # doubles the width (half-width * 2)
    hh <- (v$y1 - v$y0)
    nx0 <- cx - hw
    nx1 <- cx + hw
    ny0 <- cy - hh
    ny1 <- cy + hh
    if (nx0 <= 0.5 && nx1 >= d[["w"]] + 0.5 &&
        ny0 <= 0.5 && ny1 >= d[["h"]] + 0.5) {
      crop_state$view <- NULL
    } else {
      crop_state$view <- list(
        x0 = max(0.5, nx0), x1 = min(d[["w"]] + 0.5, nx1),
        y0 = max(0.5, ny0), y1 = min(d[["h"]] + 0.5, ny1)
      )
    }
    session$resetBrush("crop_brush")
  })
  
  observeEvent(input$crop_zoom_reset, {
    crop_state$view <- NULL
    session$resetBrush("crop_brush")
    crop_state$message <- "Zoom reset: showing the full image."
  })
  
  observeEvent(list(input$crop_ref, input$crop_files), {
    crop_state$view <- NULL
  }, ignoreInit = TRUE)
  
  observeEvent(input$crop_apply, {
    tryCatch({
      if (!length(crop_state$images)) stop("Upload field images first.")
      br <- input$crop_brush
      if (is.null(br)) stop("Drag a rectangle on the image first.")
      
      d <- crop_dims(crop_ref_image())
      W <- d[["w"]]
      H <- d[["h"]]
      
      x0 <- max(1L, as.integer(round(br$xmin)))
      x1 <- min(W, as.integer(round(br$xmax)))
      y0 <- max(1L, as.integer(round(H - br$ymax + 1)))
      y1 <- min(H, as.integer(round(H - br$ymin + 1)))
      
      if (x1 - x0 + 1L < 7L || y1 - y0 + 1L < 7L) {
        stop(sprintf(
          paste0(
            "Rectangle must be at least 7 x 7 pixels. Your selection is ",
            "%d x %d pixels (image is %d x %d). Brush x: %.1f-%.1f, y: %.1f-%.1f."
          ),
          x1 - x0 + 1L, y1 - y0 + 1L, W, H,
          br$xmin, br$xmax, br$ymin, br$ymax
        ))
      }
      
      crop_state$rect <- list(
        x0 = x0, x1 = x1, y0 = y0, y1 = y1,
        w = W, h = H, ref = crop_ref_name(), id = NULL
      )
      crop_state$message <- sprintf(
        "Rectangle: x %d-%d, y %d-%d (%d x %d pixels). Enter a nucleus number and click Save crop.",
        x0, x1, y0, y1, x1 - x0 + 1L, y1 - y0 + 1L
      )
    }, error = notify_error)
  })
  
  # Live readout of the dragged box so the size in image pixels is visible
  # before clicking Apply rectangle.
  observeEvent(input$crop_brush, {
    br <- input$crop_brush
    if (is.null(br) || !length(crop_state$images)) return()
    d <- crop_dims(crop_ref_image())
    crop_state$message <- sprintf(
      "Dragged box: about %d x %d pixels (image is %d x %d). Click Apply rectangle.",
      as.integer(round(br$xmax - br$xmin)),
      as.integer(round(br$ymax - br$ymin)),
      d[["w"]], d[["h"]]
    )
  }, ignoreNULL = TRUE)
  
  observeEvent(input$crop_clear, {
    crop_state$rect <- NULL
    crop_state$message <- "Rectangle cleared."
  })
  
  observeEvent(input$crop_save, {
    tryCatch({
      r <- crop_state$rect
      if (is.null(r)) stop("Apply a rectangle first.")
      
      id <- input$crop_id
      if (length(id) != 1L || !is.finite(id) || id < 1 || id != floor(id)) {
        stop("Nucleus number must be a positive whole number.")
      }
      id <- as.integer(id)
      
      used <- vapply(crop_state$saved, function(s) s$id, integer(1))
      if (id %in% used) {
        stop("Nucleus number ", id, " is already used; choose another.")
      }
      
      r$id <- id
      crop_state$saved[[length(crop_state$saved) + 1L]] <- r
      crop_state$rect <- NULL
      updateNumericInput(session, "crop_id", value = id + 1L)
      crop_state$message <- sprintf(
        "Saved nucleus %d. %d crop(s) saved.", id, length(crop_state$saved)
      )
    }, error = notify_error)
  })
  
  observeEvent(input$crop_remove_last, {
    n <- length(crop_state$saved)
    if (n) {
      crop_state$saved[[n]] <- NULL
      crop_state$message <- sprintf("Removed last crop. %d saved.", n - 1L)
    }
  })
  
  observeEvent(input$crop_clear_all, {
    crop_state$saved <- list()
    crop_state$rect <- NULL
    updateNumericInput(session, "crop_id", value = 1)
    crop_state$message <- "All saved crops cleared."
  })
  
  output$crop_preview <- renderPlot({
    r <- crop_state$rect
    if (is.null(r) || !length(crop_state$images)) {
      empty_plot("Apply a rectangle to preview it.")
      return(invisible(NULL))
    }
    plane <- suppressWarnings(as.integer(input$crop_plane %||% 1L))
    if (length(plane) != 1L || is.na(plane)) plane <- 1L
    print(crop_plot(crop_subimage(crop_ref_image(), r), plane))
  })
  
  crop_saved_table <- reactive({
    saved <- crop_state$saved
    if (!length(saved)) {
      return(data.frame(
        nucleus_id = integer(), x_min = integer(), x_max = integer(),
        y_min = integer(), y_max = integer(),
        width_px = integer(), height_px = integer(),
        reference_image = character(), stringsAsFactors = FALSE
      ))
    }
    do.call(rbind, lapply(saved, function(s) {
      data.frame(
        nucleus_id = s$id,
        x_min = s$x0, x_max = s$x1, y_min = s$y0, y_max = s$y1,
        width_px = s$x1 - s$x0 + 1L, height_px = s$y1 - s$y0 + 1L,
        reference_image = s$ref,
        stringsAsFactors = FALSE
      )
    }))
  })
  
  output$crop_table <- renderDT({
    DT::datatable(
      crop_saved_table(), rownames = FALSE, selection = "none",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  # Writes every saved crop for every loaded image into `dir`.
  # Returns the per-file log.
  # `plane` = NULL saves all planes; an integer saves only that plane of
  # multi-plane images. `overwrite` replaces existing files of the same name.
  crop_write_all <- function(dir, plane = NULL, overwrite = FALSE) {
    imgs <- crop_state$images
    saved <- crop_state$saved
    rows <- list()
    
    for (s in saved) {
      for (nm in names(imgs)) {
        im <- imgs[[nm]]
        d <- crop_dims(im)
        out_file <- ""
        status <- "written"
        
        if (d[["w"]] != s$w || d[["h"]] != s$h) {
          status <- sprintf(
            "skipped: dimensions %d x %d differ from reference %d x %d",
            d[["w"]], d[["h"]], s$w, s$h
          )
        } else {
          stem <- gsub(
            "[^A-Za-z0-9_.-]+", "_",
            tools::file_path_sans_ext(basename(nm))
          )
          out_file <- sprintf("%s_%d.tif", stem, s$id)
          
          if (!overwrite && file.exists(file.path(dir, out_file))) {
            status <- "skipped: output file name already exists"
            out_file <- ""
          } else {
            result <- tryCatch({
              EBImage::writeImage(
                crop_subimage(im, s, plane),
                file.path(dir, out_file),
                type = "tiff",
                bits.per.sample = 16L
              )
              NULL
            }, error = function(e) conditionMessage(e))
            if (!is.null(result)) {
              status <- paste("failed:", result)
              out_file <- ""
            }
          }
        }
        
        rows[[length(rows) + 1L]] <- data.frame(
          nucleus_id = s$id,
          source_file = nm,
          output_file = out_file,
          x_min = s$x0, x_max = s$x1, y_min = s$y0, y_max = s$y1,
          width_px = s$x1 - s$x0 + 1L,
          height_px = s$y1 - s$y0 + 1L,
          reference_image = s$ref,
          planes_saved = {
            dd <- dim(EBImage::imageData(im))
            multi <- length(dd) == 3L &&
              !(EBImage::colorMode(im) == EBImage::Color && dd[3] >= 3L)
            if (!multi) "all (2D or RGB image)"
            else if (is.null(plane)) paste0("all ", dd[3], " planes")
            else paste0("plane ", min(max(plane, 1L), dd[3]), " of ", dd[3])
          },
          status = status,
          app_version = APP_VERSION,
          coordinate_convention = paste(
            "EBImage: x = column, y = row, 1-based, y = 1 is top row"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
    
    bind_rows(rows)
  }
  
  output$crop_dl_zip <- downloadHandler(
    filename = function() paste0("GGProfiler_crops_", Sys.Date(), ".zip"),
    content = function(file) {
      if (!length(crop_state$saved)) stop("Save at least one crop first.")
      
      work <- tempfile("ggprofiler_crops_")
      dir.create(work)
      on.exit(unlink(work, recursive = TRUE), add = TRUE)
      
      plane <- if (isTRUE(input$crop_single_plane)) crop_current_plane() else NULL
      log <- crop_write_all(work, plane = plane)
      utils::write.csv(
        log, file.path(work, "crop_log.csv"),
        row.names = FALSE, na = ""
      )
      
      members <- list.files(work)
      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zipr(file, files = members, root = work)
      } else {
        old <- setwd(work)
        on.exit(setwd(old), add = TRUE)
        status <- utils::zip(file, files = members, flags = "-q")
        if (!identical(as.integer(status), 0L)) {
          stop("ZIP creation failed. Install the 'zip' package: install.packages('zip')")
        }
      }
    },
    contentType = "application/zip"
  )
  
  # ---- Save directly to a chosen folder ----------------------------------
  # Works when the app runs on your own computer (local R / RStudio). If the
  # app is hosted on a server, a folder here would be on the server, so use
  # the ZIP download instead.
  crop_pick_folder <- function() {
    sysname <- Sys.info()[["sysname"]]
    path <- tryCatch({
      if (identical(sysname, "Darwin")) {
        out <- suppressWarnings(system2(
          "osascript",
          c("-e", shQuote(
            'POSIX path of (choose folder with prompt "Choose output folder")'
          )),
          stdout = TRUE, stderr = FALSE
        ))
        if (length(out) >= 1L) trimws(out[1]) else ""
      } else if (identical(sysname, "Windows")) {
        p <- utils::choose.dir(caption = "Choose output folder")
        if (is.na(p)) "" else p
      } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
                 rstudioapi::isAvailable()) {
        rstudioapi::selectDirectory("Choose output folder") %||% ""
      } else if (requireNamespace("tcltk", quietly = TRUE)) {
        p <- as.character(tcltk::tk_choose.dir())
        if (length(p) != 1L || is.na(p)) "" else p
      } else ""
    }, error = function(e) "")
    path
  }
  
  observeEvent(input$crop_browse_dir, {
    p <- crop_pick_folder()
    if (nzchar(p)) {
      updateTextInput(session, "crop_out_dir", value = p)
    } else {
      showNotification(
        paste(
          "No folder chosen (or a folder dialog is unavailable here).",
          "Type or paste the folder path instead."
        ),
        type = "warning", duration = 10
      )
    }
  })
  
  observeEvent(input$crop_save_dir, {
    tryCatch({
      if (!length(crop_state$saved)) stop("Save at least one crop first.")
      
      dir <- trimws(input$crop_out_dir %||% "")
      if (!nzchar(dir)) stop("Choose or enter an output folder first.")
      dir <- path.expand(dir)
      
      if (!dir.exists(dir)) {
        ok <- suppressWarnings(dir.create(dir, recursive = TRUE))
        if (!ok) stop("Could not create folder: ", dir)
      }
      if (file.access(dir, 2) != 0) stop("Folder is not writable: ", dir)
      
      plane <- if (isTRUE(input$crop_single_plane)) crop_current_plane() else NULL
      log <- crop_write_all(
        dir, plane = plane, overwrite = isTRUE(input$crop_overwrite)
      )
      
      log_file <- sprintf("crop_log_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
      utils::write.csv(log, file.path(dir, log_file), row.names = FALSE, na = "")
      
      n_ok <- sum(log$status == "written")
      n_bad <- nrow(log) - n_ok
      crop_state$message <- sprintf(
        "Wrote %d file(s) to:\n%s\n%s%s",
        n_ok, normalizePath(dir),
        paste0("Log: ", log_file),
        if (n_bad) sprintf("\n%d file(s) skipped or failed; see the log.", n_bad) else ""
      )
      showNotification(
        sprintf("%d cropped file(s) written to %s", n_ok, dir),
        type = if (n_bad) "warning" else "message", duration = 10
      )
    }, error = notify_error)
  })
  
  output$crop_dl_log <- downloadHandler(
    filename = function() paste0("GGProfiler_crop_coordinates_", Sys.Date(), ".csv"),
    content = function(file) {
      table <- crop_saved_table()
      if (!nrow(table)) stop("No saved crops to export.")
      table$coordinate_convention <- paste(
        "EBImage: x = column, y = row, 1-based, y = 1 is top row"
      )
      write.csv(table, file, row.names = FALSE, na = "")
    },
    contentType = "text/csv"
  )
  
  # ---- End Crop nuclei module ----------------------------------------------
  
  # ---- Figure export: display only -----------------------------------------
  # This section reads raw channels independently.
  # It never modifies analysis arrays, segmentation or measurements.
  
  session$onFlushed(function() {
    appendTab(
      inputId = "main_tabs",
      tab = tabPanel(
        "Figure export",
        value = "figure_export",
        
        tags$details(
          style = "margin-bottom: 18px;",
          tags$summary(
            style = paste(
              "cursor: pointer;",
              "color: #0f766e;",
              "font-weight: 600;"
            ),
            "Figure export instructions"
          ),
          tags$div(
            style = "padding-top: 12px; line-height: 1.6;",
            tags$p(
              "Select a paired sample uploaded in Conditions. ",
              "Segmentation and calibration are not required."
            ),
            tags$p(
              "Display minimum maps to black; display maximum maps ",
              "to full channel color. Values outside this range are ",
              "clipped for display only. Moving both limits changes ",
              "brightness; reducing their separation increases contrast."
            ),
            tags$p(
              "Auto uses the selected image's 0.5th and 99.5th ",
              "intensity percentiles. Reset uses its full intensity ",
              "range. Both buttons affect every configured channel."
            ),
            tags$p(
              "Limits and colors remain unchanged when switching ",
              "samples. For comparable figures, use the same limits ",
              "for each marker across images acquired comparably. ",
              "Applying Auto separately to every sample can hide ",
              "real intensity differences."
            ),
            tags$p(
              "The histogram shows the selected channel before ",
              "display adjustments. Dashed lines mark the display limits. ",
              "Clipping percentages count values strictly below or ",
              "above those limits."
            ),
            tags$p(
              "PDF width controls physical page width; height follows ",
              "the image aspect ratio. The PDF embeds a display RGB ",
              "raster at the original pixel dimensions, without ",
              "interpolation. It does not add resolution or preserve ",
              "the original quantitative intensity precision."
            ),
            tags$p(
              "All controls here affect figure display only. Raw files, ",
              "segmentation, QC and exported measurements are unchanged."
            )
          )
        ),
        
        fluidRow(
          column(
            3,
            selectInput(
              "fig_sample",
              "Sample",
              choices = character()
            ),
            
            actionButton(
              "fig_auto",
              "Auto contrast: all channels"
            ),
            actionButton(
              "fig_reset",
              "Reset to full range"
            ),
            
            helpText(
              "Initial display limits are 0–1. Use Reset or enter ",
              "appropriate limits for your loaded image."
            ),
            
            div(
              style = paste(
                "height: 420px; max-height: 55vh; overflow-y: auto;",
                "overflow-x: hidden; overscroll-behavior-y: contain;",
                "padding: 4px 12px 4px 4px;"
              ),
              uiOutput("fig_channel_controls")
            ),
            
            numericInput(
              "fig_width_mm",
              "PDF width (mm)",
              value = 90,
              min = 10,
              max = 500,
              step = 1
            ),
            
            downloadButton(
              "fig_download_pdf",
              "Download figure PDF"
            )
          ),
          
          column(
            9,
            div(
              style = "width: 240px; max-width: 100%; margin: 0 auto;",
              plotOutput(
                "fig_preview",
                width = "100%",
                height = "220px"
              )
            ),
            verbatimTextOutput("fig_dimensions"),
            
            selectInput(
              "fig_hist_channel",
              "Histogram channel",
              choices = character()
            ),
            plotOutput(
              "fig_histogram",
              height = "220px"
            ),
            
            h4("Display clipping"),
            tableOutput("fig_clipping")
          )
        )
      ),
      select = FALSE,
      session = session
    )
  }, once = TRUE)
  
  
  # Populate the sample selector from the existing pairing manifest.
  observe({
    samples <- manifest()$samples
    
    choices <- if (length(samples)) {
      labels <- vapply(
        samples,
        function(s) paste(s$condition, s$stem, sep = " / "),
        character(1)
      )
      
      stats::setNames(
        as.character(seq_along(samples)),
        labels
      )
    } else {
      character()
    }
    
    previous <- isolate(input$fig_sample)
    
    selected <- if (
      length(previous) == 1L &&
      previous %in% unname(choices)
    ) {
      previous
    } else if (length(choices)) {
      unname(choices[1])
    } else {
      character()
    }
    
    updateSelectInput(
      session,
      "fig_sample",
      choices = choices,
      selected = selected
    )
  })
  
  
  # Controls depend on channel setup, not the selected sample.
  # Therefore changing sample does not automatically rescale it.
  output$fig_channel_controls <- renderUI({
    cfg <- rv$setup
    req(length(cfg) > 0L)
    
    tagList(lapply(seq_along(cfg), function(i) {
      wellPanel(
        h4(cfg[[i]]$name),
        
        checkboxInput(
          paste0("fig_show_", i),
          "Show channel",
          value = TRUE
        ),
        
        selectInput(
          paste0("fig_colour_", i),
          "Figure color",
          choices = COLORS,
          selected = cfg[[i]]$color
        ),
        
        numericInput(
          paste0("fig_min_", i),
          "Base display minimum (black at neutral sliders)",
          value = 0,
          step = 0.001
        ),
        
        numericInput(
          paste0("fig_max_", i),
          "Base display maximum (full colour at neutral sliders)",
          value = 1,
          step = 0.001
        ),
        sliderInput(
          paste0("fig_brightness_", i),
          "Brightness (%)",
          min = -100, max = 100, value = 0, step = 1,
          post = "%"
        ),
        sliderInput(
          paste0("fig_contrast_", i),
          "Contrast (multiplier)",
          min = 0.1, max = 5, value = 1, step = 0.1,
          post = "x"
        )
      )
    }))
  })
  
  
  observe({
    cfg <- rv$setup
    req(length(cfg) > 0L)
    
    updateSelectInput(
      session,
      "fig_hist_channel",
      choices = stats::setNames(
        as.character(seq_along(cfg)),
        vapply(cfg, `[[`, character(1), "name")
      ),
      selected = "1"
    )
  })
  
  
  # Read original channels directly, without running segmentation.
  fig_raw <- reactive({
    cfg <- rv$setup
    samples <- manifest()$samples
    
    req(length(cfg) > 0L)
    req(length(samples) > 0L)
    req(input$fig_sample)
    
    index <- suppressWarnings(as.integer(input$fig_sample))
    
    validate(need(
      length(index) == 1L &&
        is.finite(index) &&
        index >= 1L &&
        index <= length(samples),
      "Select a valid sample."
    ))
    
    sample <- samples[[index]]
    
    result <- tryCatch({
      mats <- lapply(seq_along(cfg), function(i) {
        read_channel(sample$paths[i], cfg[[i]]$src)
      })
      
      dimensions <- lapply(mats, dim)
      
      if (!all(vapply(
        dimensions,
        identical,
        logical(1),
        dimensions[[1]]
      ))) {
        stop(
          "Channel dimensions differ. ",
          "Use matching, aligned image crops."
        )
      }
      
      list(
        mats = mats,
        cfg = cfg,
        sample = sample,
        width = nrow(mats[[1]]),
        height = ncol(mats[[1]])
      )
    }, error = function(e) {
      list(error = conditionMessage(e))
    })
    
    validate(need(
      is.null(result$error),
      paste("Cannot load figure:", result$error)
    ))
    
    result
  })
  
  
  # Ensure a strictly positive display window, including constant images.
  fig_valid_range <- function(m, auto = FALSE) {
    limits <- if (auto) {
      as.numeric(stats::quantile(
        m,
        probs = c(0.005, 0.995),
        names = FALSE
      ))
    } else {
      range(m)
    }
    
    if (limits[2] <= limits[1]) {
      limits <- range(m)
    }
    
    if (limits[2] <= limits[1]) {
      limits[2] <- limits[1] +
        max(abs(limits[1]), 1) * 1e-6
    }
    
    limits
  }
  
  
  fig_update_ranges <- function(auto = FALSE) {
    raw <- fig_raw()
    
    for (i in seq_along(raw$mats)) {
      limits <- fig_valid_range(raw$mats[[i]], auto)
      
      updateNumericInput(
        session,
        paste0("fig_min_", i),
        value = limits[1]
      )
      
      updateNumericInput(
        session,
        paste0("fig_max_", i),
        value = limits[2]
      )
      updateSliderInput(
        session, paste0("fig_brightness_", i), value = 0
      )
      updateSliderInput(
        session, paste0("fig_contrast_", i), value = 1
      )
    }
  }
  
  observeEvent(input$fig_auto, {
    fig_update_ranges(auto = TRUE)
  })
  
  observeEvent(input$fig_reset, {
    fig_update_ranges(auto = FALSE)
  })
  
  
  # Read and validate current display controls.
  fig_settings <- reactive({
    raw <- fig_raw()
    
    settings <- lapply(seq_along(raw$mats), function(i) {
      low <- input[[paste0("fig_min_", i)]]
      high <- input[[paste0("fig_max_", i)]]
      colour <- input[[paste0("fig_colour_", i)]]
      visible <- input[[paste0("fig_show_", i)]]
      
      req(
        !is.null(low),
        !is.null(high),
        !is.null(colour),
        !is.null(visible)
      )
      
      validate(need(
        length(low) == 1L &&
          length(high) == 1L &&
          is.finite(low) &&
          is.finite(high) &&
          high > low,
        paste(
          raw$cfg[[i]]$name,
          ": display maximum must exceed display minimum."
        )
      ))
      
      validate(need(
        colour %in% unname(COLORS),
        "Select a valid figure color."
      ))
      
      brightness <- input[[paste0("fig_brightness_", i)]]
      contrast <- input[[paste0("fig_contrast_", i)]]
      req(!is.null(brightness), !is.null(contrast))
      validate(need(
        length(brightness) == 1L && is.finite(brightness) &&
          brightness >= -100 && brightness <= 100 &&
          length(contrast) == 1L && is.finite(contrast) &&
          contrast >= 0.1 && contrast <= 5,
        "Brightness must be -100 to 100%; contrast must be 0.1 to 5."
      ))
      
      # Shift the base-window center for brightness and change its width
      # for contrast. Never modify the original channel matrix.
      span <- high - low
      center <- low + span / 2 - (brightness / 100) * span
      half_width <- span / (2 * contrast)
      adjusted_low <- center - half_width
      adjusted_high <- center + half_width
      validate(need(
        is.finite(adjusted_low) && is.finite(adjusted_high) &&
          adjusted_high > adjusted_low,
        "Adjusted display limits must be finite and distinct."
      ))
      
      list(
        low = adjusted_low,
        high = adjusted_high,
        colour = colour,
        visible = isTRUE(visible)
      )
    })
    
    settings
  })
  
  
  # Construct a NEW display raster.
  # EBImage is [x,y]; raster rows/columns are [y,x].
  fig_raster <- reactive({
    raw <- fig_raw()
    settings <- fig_settings()
    
    validate(need(
      any(vapply(settings, `[[`, logical(1), "visible")),
      "Enable at least one channel."
    ))
    
    rgb_array <- array(
      0,
      dim = c(raw$height, raw$width, 3L)
    )
    
    for (i in seq_along(raw$mats)) {
      s <- settings[[i]]
      if (!s$visible) next
      
      intensity <- t(
        (raw$mats[[i]] - s$low) / (s$high - s$low)
      )
      
      intensity[intensity < 0] <- 0
      intensity[intensity > 1] <- 1
      
      colour_rgb <- as.numeric(
        grDevices::col2rgb(s$colour)
      ) / 255
      
      for (k in seq_len(3L)) {
        rgb_array[, , k] <- pmin(
          1,
          rgb_array[, , k] + intensity * colour_rgb[k]
        )
      }
    }
    
    grDevices::as.raster(rgb_array)
  })
  
  
  # Use physical device dimensions to preserve the image aspect ratio.
  # Raster row 1 remains the top image row.
  fig_draw <- function(raster, width_px, height_px) {
    device <- grDevices::dev.size("in")
    image_ratio <- width_px / height_px
    device_ratio <- device[1] / device[2]
    
    if (image_ratio >= device_ratio) {
      viewport_width <- 1
      viewport_height <- device_ratio / image_ratio
    } else {
      viewport_width <- image_ratio / device_ratio
      viewport_height <- 1
    }
    
    grid::grid.newpage()
    
    grid::grid.raster(
      raster,
      x = 0.5,
      y = 0.5,
      width = grid::unit(viewport_width, "npc"),
      height = grid::unit(viewport_height, "npc"),
      interpolate = FALSE
    )
  }
  
  
  output$fig_preview <- renderPlot({
    raw <- fig_raw()
    
    fig_draw(
      fig_raster(),
      raw$width,
      raw$height
    )
  }, res = 96)
  
  
  output$fig_dimensions <- renderText({
    raw <- fig_raw()
    width_mm <- input$fig_width_mm
    
    req(
      length(width_mm) == 1L,
      is.finite(width_mm),
      width_mm > 0
    )
    
    height_mm <- width_mm * raw$height / raw$width
    effective_ppi <- raw$width / (width_mm / 25.4)
    
    sprintf(
      paste0(
        "Native image: %d x %d pixels\n",
        "PDF size: %.1f x %.1f mm\n",
        "Effective raster resolution at this size: %.0f ppi"
      ),
      raw$width,
      raw$height,
      width_mm,
      height_mm,
      effective_ppi
    )
  })
  
  
  output$fig_histogram <- renderPlot({
    raw <- fig_raw()
    settings <- fig_settings()
    
    req(input$fig_hist_channel)
    i <- suppressWarnings(as.integer(input$fig_hist_channel))
    
    req(
      length(i) == 1L,
      is.finite(i),
      i >= 1L,
      i <= length(raw$mats)
    )
    
    graphics::par(mar = c(4, 4, 2, 1))
    
    graphics::hist(
      as.numeric(raw$mats[[i]]),
      breaks = 100,
      col = settings[[i]]$colour,
      border = NA,
      main = raw$cfg[[i]]$name,
      xlab = "Loaded raw intensity",
      ylab = "Pixel count"
    )
    
    graphics::abline(
      v = c(settings[[i]]$low, settings[[i]]$high),
      col = "#333333",
      lty = 2,
      lwd = 2
    )
  })
  
  
  output$fig_clipping <- renderTable({
    raw <- fig_raw()
    settings <- fig_settings()
    
    do.call(rbind, lapply(seq_along(raw$mats), function(i) {
      m <- raw$mats[[i]]
      s <- settings[[i]]
      
      data.frame(
        Channel = raw$cfg[[i]]$name,
        Visible = s$visible,
        Minimum = s$low,
        Maximum = s$high,
        Below_min_pct = 100 * mean(m < s$low),
        Above_max_pct = 100 * mean(m > s$high),
        check.names = FALSE
      )
    }))
  }, digits = 4, rownames = FALSE)
  
  
  output$fig_download_pdf <- downloadHandler(
    filename = function() {
      raw <- fig_raw()
      
      stem <- paste(
        raw$sample$condition,
        raw$sample$stem,
        sep = "_"
      )
      
      stem <- gsub("[^A-Za-z0-9_-]+", "_", stem)
      
      paste0(stem, "_figure.pdf")
    },
    
    content = function(file) {
      raw <- fig_raw()
      raster <- fig_raster()
      width_mm <- input$fig_width_mm
      
      if (
        length(width_mm) != 1L ||
        !is.finite(width_mm) ||
        width_mm < 10 ||
        width_mm > 500
      ) {
        stop("PDF width must be between 10 and 500 mm.")
      }
      
      width_in <- width_mm / 25.4
      height_in <- width_in * raw$height / raw$width
      
      grDevices::pdf(
        file = file,
        width = width_in,
        height = height_in,
        paper = "special",
        onefile = FALSE,
        bg = "white",
        useDingbats = FALSE,
        compress = TRUE
      )
      
      pdf_device <- grDevices::dev.cur()
      
      on.exit(
        grDevices::dev.off(which = pdf_device),
        add = TRUE
      )
      
      fig_draw(
        raster,
        raw$width,
        raw$height
      )
    },
    
    contentType = "application/pdf"
  )
  
  # ---- End Figure export --------------------------------------------------
  
  # ---- Downloads -----------------------------------------------------------
  
  output$dl_manifest <- downloadHandler(
    filename = function() paste0("GGProfiler_manifest_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(manifest()$table, file, row.names = FALSE, na = "")
    }
  )
  
  for (field in c(
    features = "features",
    agg = "aggregated",
    objects = "objects",
    log = "log"
  )) {
    local({
      run_field <- unname(field)
      output_id <- switch(
        run_field,
        features = "dl_features",
        aggregated = "dl_agg",
        objects = "dl_objects",
        log = "dl_log"
      )
      output[[output_id]] <- downloadHandler(
        filename = function() {
          paste0("GGProfiler_", run_field, "_", Sys.Date(), ".csv")
        },
        content = function(file) {
          run <- run_current()
          req(run)
          write.csv(run[[run_field]], file, row.names = FALSE, na = "")
        }
      )
    })
  }
  
  output$dl_stats <- downloadHandler(
    filename = function() paste0("GGProfiler_statistics_", Sys.Date(), ".csv"),
    content = function(file) {
      table <- stats_current()
      req(table)
      write.csv(table, file, row.names = FALSE, na = "")
    }
  )
  
  output$dl_session <- downloadHandler(
    filename = function() paste0("GGProfiler_snapshot_", Sys.Date(), ".rds"),
    content = function(file) {
      run <- run_current()
      req(run)
      
      saveRDS(
        list(
          app_version = APP_VERSION,
          created = Sys.time(),
          run = run,
          statistics = stats_current(),
          session_info = sessionInfo(),
          coordinate_convention = "EBImage [x,y]; integer pixel centers",
          notes = c(
            "Original uploaded image files are not embedded.",
            "Temporary input paths are provenance only and may expire.",
            "Spatial units are pixels.",
            "Legacy feature names and chance model retained."
          )
        ),
        file = file
      )
    }
  )
}

shinyApp(ui, server)
