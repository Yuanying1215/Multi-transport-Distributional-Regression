# ==============================================================================
# Script: reference_sensitivity.R
# Description: Reference Sensitivity Analysis for the Real Data Application.
#              Uses mortality data from the Human Mortality Database (HMD) to
#              predict male age-at-death distributions in 2010 from male and
#              female distributions in 2005. Independently refits MTDR under
#              three reference choices using leave-one-out cross-validation.
#              Supports Section S5 and Tables S7-S8 of the Supplementary Material.
#
# Method: Multi-transport Distributional Regression (MTDR)
#         Calls MOT2() in Functions.R; OT and GOT models are not fitted here.
#
# Reference Distributions (on the common age interval [0, 100]):
#   1. FM: Wasserstein Frechet mean of the training response distributions.
#   2. U:  Uniform distribution on [0, 100].
#   3. TN: N(50, 25^2) truncated to [0, 100].
#
# Data:
#   - MortMale.RData: Male mortality data.
#   - MortFemale.RData: Female mortality data.
#
# Dependencies:
#   - Functions.R (the original fitting functions are not modified on disk).
#   - Packages: pracma, fdadensity.
#
# Usage:
#   Set the working directory to the folder containing this script, Functions.R,
#   and the two data files. Alternatively, specify --data-dir=/path/to/Real_data.
#
#   Terminal (full analysis):
#     Rscript --vanilla reference_sensitivity.R
#
#   R console / RStudio (source alone defines functions but does not run fits):
#     source("reference_sensitivity.R")
#     opt <- parse_options(character())
#     run_analysis(opt)
#
#   Quick check (Bulgaria only; NOT the full analysis):
#     Rscript --vanilla reference_sensitivity.R --folds=5 --full-fit=false
#
#   All command-line options use --name=value. Defaults:
#     --data-dir=.          --folds=all       --full-fit=true
#     --max-iter=500        --tol=1e-7        --quantiles=50
#     --map-grid=101        --resume=false
#   Output defaults to a new reference_results_YYYYMMDD_HHMMSS folder in the
#   working directory; override with --output-dir=/path/to/results.
#   To resume, supply the same output directory and analysis options together
#   with --resume=true. Saved fits, including flagged fits, are not rerun.
#
# Output:
#   1. loo_by_country.csv: Held-out prediction errors, reference sensitivity,
#      estimated weights, and fitting diagnostics for each country/reference.
#   2. loo_summary.csv: Average Wasserstein Distance (AWD) and descriptive
#      summaries across all leave-one-out folds for each reference.
#   3. full_sample_weights.csv: Separate full-data fits (if full_fit = TRUE).
#   4. fits/ and logs/: Numeric checkpoints and complete optimization logs.
#   5. prepared_data.rds, config.rds/txt, sessionInfo.txt, Functions_used.R:
#      Processed data, settings, input checksums, and reproducibility records.
#
# Interpretation:
#   - Weight order: alpha0 = reference, alpha1 = male 2005, alpha2 = female 2005.
#   - All Wasserstein distances are in years. Cross-reference differences use
#     the FM fit in the SAME fold as the baseline.
#   - Fold-level weights are common-model estimates from the training countries,
#     not country-specific coefficients. Across-fold SDs are not standard errors.
#   - max_iter marks fits that did not meet the original stopping criterion.
#     Such fits are retained; failed fits produce NA rather than being omitted.
#   - ready_for_review checks completeness and diagnostics, not invariance or
#     global optimality. Review flagged fits before reporting numerical results.
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 1: Analysis Options
# ------------------------------------------------------------------------------
# Defaults follow the original real-data implementation, apart from the three
# reference choices. File paths are relative to the working directory by default.

parse_options <- function(args) {
  opt <- list(
    data_dir = ".",
    output_dir = file.path(getwd(), paste0("reference_results_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
    folds = "all", full_fit = TRUE, resume = FALSE,
    max_iter = 500L, tol = 1e-7, quantiles = 50L, map_grid = 101L
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) stop("Use --name=value: ", arg)
    key <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", arg))
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(opt)) stop("Unknown option: ", key)
    if (key %in% c("full_fit", "resume")) {
      if (!tolower(value) %in% c("true", "false")) stop(key, " must be true or false")
      value <- tolower(value) == "true"
    }
    if (key %in% c("max_iter", "quantiles", "map_grid")) {
      value <- suppressWarnings(as.numeric(value))
      if (!is.finite(value) || value < 3 || value != floor(value)) stop("Invalid integer: ", key)
      value <- as.integer(value)
    }
    if (key == "tol") {
      value <- suppressWarnings(as.numeric(value))
      if (!is.finite(value) || value <= 0) stop("Invalid tolerance")
    }
    opt[[key]] <- value
  }
  opt
}

# ------------------------------------------------------------------------------
# Step 2: Wasserstein Distance and Quantile Validation
# ------------------------------------------------------------------------------
# Evaluate one-dimensional W2 distances by trapezoidal integration on the shared
# quantile grid, and reject nonfinite or nonmonotone numerical representations.

trap_integral <- function(x, y) sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
w2_quantile <- function(a, b, q) {
  stopifnot(
    length(a) == length(q), length(b) == length(q),
    all(is.finite(a)), all(is.finite(b))
  )
  sqrt(max(0, trap_integral(q, (a - b)^2)))
}
check_quantile <- function(q, label, bounds = c(0, 100)) {
  if (any(!is.finite(q)) || any(diff(q) < -1e-8) ||
    min(q) < bounds[1] - 1e-8 || max(q) > bounds[2] + 1e-8) {
    stop("Invalid quantile/map: ", label)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Step 3: Data Pre-processing (Male and Female)
# ------------------------------------------------------------------------------
# Same age cutoff, density normalization and dens2quantile as Results_figure.R.
# Only the required years are converted; no pooled response-dependent preprocessing.
load_mortality <- function(path, q, years) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!all(c("age", "country", "mort", "year") %in% ls(e))) stop("Incomplete RData: ", path)
  countries <- as.character(e$country)
  if (anyNA(countries) || anyDuplicated(countries) || length(countries) != length(e$mort)) {
    stop("Invalid country identifiers")
  }
  keep <- which(e$age <= 100)
  age <- e$age[keep]
  if (min(age) != 0 || max(age) != 100 || any(diff(age) <= 0)) stop("Expected age grid [0,100]")
  result <- lapply(years, function(yr) {
    mat <- vapply(seq_along(countries), function(i) {
      col <- which(e$year[[i]] == yr)
      if (length(col) != 1L) stop("Missing/duplicate year for ", countries[i], ": ", yr)
      density <- e$mort[[i]][keep, col]
      if (any(!is.finite(density)) || any(density < 0)) stop("Invalid mortality density")
      mass <- pracma::trapz(age, density)
      if (mass <= 0) stop("Zero density mass")
      quant <- as.numeric(fdadensity::dens2quantile(density / mass, dSup = age, qSup = q))
      check_quantile(quant, paste(countries[i], yr))
      quant
    }, numeric(length(q)))
    colnames(mat) <- countries
    mat
  })
  names(result) <- as.character(years)
  list(country = countries, quantiles = result)
}

prepare_data <- function(data_dir, q) {
  male <- load_mortality(file.path(data_dir, "MortMale.RData"), q, c(2005, 2010))
  female <- load_mortality(file.path(data_dir, "MortFemale.RData"), q, 2005)
  if (!setequal(male$country, female$country)) stop("Male/female country sets differ")
  # Join explicitly by country, never assume the two files have identical ordering.
  ix <- match(male$country, female$country)
  list(
    country = male$country, Qy = male$quantiles[["2010"]],
    Qx = list(
      male_2005 = male$quantiles[["2005"]],
      female_2005 = female$quantiles[["2005"]][, ix, drop = FALSE]
    )
  )
}

# ------------------------------------------------------------------------------
# Step 4: Reference Distributions
# ------------------------------------------------------------------------------
# The FM reference uses training responses only. Uniform and truncated-normal
# references are specified in advance and remain fixed across folds.

reference_quantiles <- function(Qy_train, q) {
  lo <- pnorm(0, mean = 50, sd = 25)
  hi <- pnorm(100, mean = 50, sd = 25)
  normal <- qnorm(lo + q * (hi - lo), mean = 50, sd = 25)
  normal[c(1L, length(q))] <- c(0, 100)
  refs <- list(
    response_mean = rowMeans(Qy_train), uniform = 100 * q,
    truncated_normal = normal
  )
  for (name in names(refs)) check_quantile(refs[[name]], name)
  refs
}

# ------------------------------------------------------------------------------
# Step 5: MTDR Fitting and Optimization Diagnostics
# ------------------------------------------------------------------------------
# Source Functions.R in a private environment and use the original MOT2 routine.
# Fit each reference independently with identity-map initialization. Preserve
# warnings, errors, iteration counts, and stopping status for subsequent review.

load_solver <- function(path, max_iter) {
  solver <- new.env(parent = globalenv())
  sys.source(path, envir = solver)
  if (!all(c("MOT2", "Main2") %in% ls(solver))) stop("Missing MOT2/Main2 in Functions.R")
  # Change only the iteration-budget default in this private environment.
  # For max_iter=500, the original algorithm and numerical settings are unchanged.
  formals(solver$Main2)$max_iter <- max_iter
  solver
}

fit_one <- function(solver, Qy, Qx, bar, q, x, tol, log_path) {
  warnings <- character()
  error <- NULL
  fit <- NULL
  elapsed <- system.time({
    trace <- capture.output({
      fit <- tryCatch(withCallingHandlers(
        solver$MOT2(Qy, Qx, bar, q, x,
          tol = tol, inti_type = 1,
          al = 0, ar = 1, zero = 1e-6
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ), error = function(e) {
        error <<- conditionMessage(e)
        NULL
      })
    })
  })[["elapsed"]]
  iter_lines <- grep("^Iteration [0-9]+:", trace, value = TRUE)
  iterations <- if (length(iter_lines)) as.integer(sub("^Iteration ([0-9]+):.*", "\\1", tail(iter_lines, 1))) else 0L
  converged <- any(grepl("^Converged\\.", trace))
  if (!is.null(fit)) {
    error <- tryCatch(
      {
        alpha <- as.numeric(fit$a_res)
        if (length(alpha) != 3L || any(!is.finite(alpha)) || any(alpha < 0) || abs(sum(alpha) - 1) > 1e-7) {
          stop("Invalid simplex weights")
        }
        if (!is.finite(fit$loss_res)) stop("Nonfinite training loss")
        maps <- vapply(fit$T_res, function(f) f(x), numeric(length(x)))
        for (j in 1:3) check_quantile(maps[, j], paste("T", j - 1))
        mu0 <- fit$T_res[[1]](bar)
        check_quantile(mu0, "transported reference")
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }
  writeLines(c(trace, paste("WARNING:", warnings), if (!is.null(error)) paste("ERROR:", error)), log_path)
  if (!is.null(error)) {
    return(list(
      status = "failed", error = error, warnings = warnings,
      iterations = iterations, elapsed = elapsed
    ))
  }
  list(
    status = if (converged) "converged" else "max_iter", error = "",
    warnings = warnings, iterations = iterations, elapsed = elapsed,
    alpha = as.numeric(fit$a_res), maps = maps, mu0 = as.numeric(mu0),
    reference_quantile = bar, training_loss = fit$loss_res,
    # Functions are used only during this run; checkpoints store numeric grids.
    functions = fit$T_res
  )
}

# ------------------------------------------------------------------------------
# Step 6: Result Tables and Descriptive Summaries
# ------------------------------------------------------------------------------
# Keep all folds in the summaries. Componentwise comparisons describe numerical
# stability and do not establish identifiability or exact parameter invariance.

make_row <- function(fit, fold, country, reference) {
  ok <- fit$status != "failed"
  a <- if (ok) fit$alpha else rep(NA_real_, 3)
  data.frame(
    fold = fold, country = country, reference = reference, status = fit$status,
    iterations = fit$iterations, elapsed_seconds = fit$elapsed,
    alpha0_reference = a[1], alpha1_male2005 = a[2], alpha2_female2005 = a[3],
    training_loss = if (ok) fit$training_loss else NA_real_,
    error = fit$error, warnings = paste(unique(fit$warnings), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

summarize_folds <- function(rows, expected_n, selected_n) {
  # Never silently omit failed countries. NA propagates into all affected metrics.
  do.call(rbind, lapply(unique(rows$reference), function(ref) {
    z <- rows[rows$reference == ref, ]
    metrics <- c(
      "prediction_W2", "prediction_difference_W2", "reference_difference_W2",
      "alpha0_reference", "alpha1_male2005", "alpha2_female2005", "alpha_difference_L2"
    )
    ans <- list(
      reference = ref, n_countries = nrow(z), n_expected = expected_n,
      scope = if (selected_n == expected_n) "all_LOOCV_folds" else "partial_test_NOT_final",
      n_converged = sum(z$status == "converged"), n_max_iter = sum(z$status == "max_iter"),
      n_failed = sum(z$status == "failed"),
      n_warnings = sum(nzchar(z$warnings)),
      ready_for_review = selected_n == expected_n && nrow(z) == expected_n &&
        all(z$status == "converged") && !any(nzchar(z$warnings)) &&
        all(vapply(z[metrics], function(v) all(is.finite(v)), logical(1)))
    )
    for (metric in metrics) {
      v <- z[[metric]]
      ans[[paste0(metric, "_mean")]] <- mean(v)
      ans[[paste0(metric, "_sd")]] <- if (length(v) > 1) sd(v) else NA_real_
    }
    ans$prediction_difference_W2_max <- max(z$prediction_difference_W2)
    as.data.frame(ans, stringsAsFactors = FALSE)
  }))
}

# ------------------------------------------------------------------------------
# Step 7: Leave-One-Out Cross-Validation and Optional Full-Data Fits
# ------------------------------------------------------------------------------
# Align countries by name and keep the male/female predictor order fixed.
# Recompute the FM reference within each training fold. Checkpoint every fit,
# then compare predictions and transported reference measures with the FM fit.
# Full-data weights are saved separately from leave-one-out summaries.

run_analysis <- function(opt) {
  for (pkg in c("pracma", "fdadensity")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
  }
  data_dir <- normalizePath(opt$data_dir, mustWork = TRUE)
  source_paths <- file.path(data_dir, c("Functions.R", "MortMale.RData", "MortFemale.RData"))
  if (!all(file.exists(source_paths))) stop("Missing source functions/data")
  q <- seq(0, 1, length.out = opt$quantiles)
  x <- seq(0, 100, length.out = opt$map_grid)
  data <- prepare_data(data_dir, q)
  n <- length(data$country)
  folds <- if (opt$folds == "all") seq_len(n) else suppressWarnings(as.numeric(strsplit(opt$folds, ",", fixed = TRUE)[[1]]))
  if (!length(folds) || anyNA(folds) || any(folds != floor(folds)) || any(!folds %in% seq_len(n)) || anyDuplicated(folds)) {
    stop("folds must be all or comma-separated distinct country indices")
  }
  folds <- as.integer(folds)
  out <- path.expand(opt$output_dir)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)) && !opt$resume) {
    stop("Output directory is nonempty. Use a new directory or --resume=true.")
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, mustWork = TRUE)
  dir.create(file.path(out, "logs"), showWarnings = FALSE)
  dir.create(file.path(out, "fits"), showWarnings = FALSE)
  config <- opt[setdiff(names(opt), c("output_dir", "resume", "data_dir"))]
  config$data_dir <- data_dir
  config$fold_indices <- folds
  config$input_md5 <- tools::md5sum(source_paths)
  config$protocol_version <- "1.0"
  config$r_version <- R.version.string
  config$package_versions <- vapply(c("pracma", "fdadensity"), function(p) as.character(utils::packageVersion(p)), character(1))
  config_path <- file.path(out, "config.rds")
  if (file.exists(config_path)) {
    if (!identical(readRDS(config_path), config)) stop("Resume configuration/data mismatch; use a new output directory")
  } else {
    saveRDS(config, config_path)
    file.copy(source_paths[1], file.path(out, "Functions_used.R"), overwrite = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(capture.output(dput(config)), file.path(out, "config.txt"))
  saveRDS(list(q = q, x = x, data = data), file.path(out, "prepared_data.rds"))
  solver <- load_solver(source_paths[1], opt$max_iter)
  refs_order <- names(reference_quantiles(data$Qy[, -1, drop = FALSE], q))
  rows <- list()
  tasks <- c(as.list(folds), if (opt$full_fit) list(0L) else list())
  for (i in tasks) {
    full <- i == 0L
    train <- if (full) seq_len(n) else setdiff(seq_len(n), i)
    label <- if (full) "full_sample" else sprintf("loo_%02d", i)
    country <- if (full) "ALL_COUNTRIES" else data$country[i]
    refs <- reference_quantiles(data$Qy[, train, drop = FALSE], q)
    fits <- list()
    for (ref in refs_order) {
      message(sprintf("[%s / %s] %s", label, country, ref))
      checkpoint <- file.path(out, "fits", paste0(label, "_", ref, ".rds"))
      if (opt$resume && file.exists(checkpoint)) {
        fit <- readRDS(checkpoint)
      } else {
        fit <- fit_one(
          solver, data$Qy[, train, drop = FALSE],
          lapply(data$Qx, function(z) z[, train, drop = FALSE]),
          refs[[ref]], q, x, opt$tol, file.path(out, "logs", paste0(label, "_", ref, ".log"))
        )
        if (fit$status != "failed" && !full) {
          inputs <- list(refs[[ref]], data$Qx[[1]][, i], data$Qx[[2]][, i])
          fit$prediction <- Reduce(`+`, lapply(1:3, function(j) fit$alpha[j] * fit$functions[[j]](inputs[[j]])))
          check_quantile(fit$prediction, paste(label, ref, "prediction"))
          fit$prediction_W2 <- w2_quantile(fit$prediction, data$Qy[, i], q)
        }
        fit$functions <- NULL
        saveRDS(fit, checkpoint)
      }
      fits[[ref]] <- fit
    }
    base <- fits$response_mean
    group_rows <- lapply(refs_order, function(ref) {
      fit <- fits[[ref]]
      row <- make_row(fit, i, country, ref)
      pair_ok <- fit$status != "failed" && base$status != "failed"
      row$reference_difference_W2 <- if (pair_ok) w2_quantile(fit$mu0, base$mu0, q) else NA_real_
      row$alpha_difference_L2 <- if (pair_ok) sqrt(sum((fit$alpha - base$alpha)^2)) else NA_real_
      if (!full) {
        row$prediction_W2 <- if (fit$status != "failed") fit$prediction_W2 else NA_real_
        row$prediction_difference_W2 <- if (pair_ok) w2_quantile(fit$prediction, base$prediction, q) else NA_real_
      }
      row
    })
    if (full) {
      write.csv(do.call(rbind, group_rows), file.path(out, "full_sample_weights.csv"), row.names = FALSE)
    } else {
      rows <- c(rows, group_rows)
      tab <- do.call(rbind, rows)
      write.csv(tab, file.path(out, "loo_by_country.csv"), row.names = FALSE)
      write.csv(summarize_folds(tab, n, length(folds)), file.path(out, "loo_summary.csv"), row.names = FALSE)
    }
  }
  summary <- summarize_folds(do.call(rbind, rows), n, length(folds))
  print(summary[, c("reference", "n_countries", "n_converged", "n_max_iter", "n_failed", "prediction_W2_mean")], row.names = FALSE)
  message("Saved to: ", out)
  if (any(!summary$ready_for_review)) message("Some results are partial/nonconverged/flagged: inspect diagnostics before reporting.")
  invisible(list(output_dir = out, summary = summary))
}

# ------------------------------------------------------------------------------
# Step 8: Script Entry Point
# ------------------------------------------------------------------------------
# Rscript runs the analysis directly. When sourced in an interactive session,
# call run_analysis(parse_options(character())) explicitly after configuration.

if (sys.nframe() == 0L) run_analysis(parse_options(commandArgs(trailingOnly = TRUE)))
