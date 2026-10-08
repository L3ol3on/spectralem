library(testthat)

# Six well-separated Voigt peaks on a flat offset, with a little noise.
six_peaks <- function(noise = 0.001, seed = 1) {
  set.seed(seed)
  params <- list(
    amp = c(10, 6, 8, 5, 9, 7),
    pos = c(150, 300, 450, 600, 750, 900),
    gwidth = rep(6, 6),
    lwidth = rep(2, 6),
    a = 0,
    b = 0.5
  )
  x <- seq(0, 1000, 1)
  y <- stats::rnorm(length(x), voigt.model(x, params), noise)
  list(x = x, y = y, params = params)
}


test_that("estimate_peak_count finds well-separated peaks", {
  data <- six_peaks()
  expect_equal(estimate_peak_count(data$x, data$y), 6)
})


test_that("estimate_peak_count ignores x order and returns at least one", {
  data <- six_peaks()
  o <- rev(seq_along(data$x))
  expect_equal(estimate_peak_count(data$x[o], data$y[o]), 6)
  expect_equal(estimate_peak_count(seq(0, 1, length.out = 50), rep(1, 50)), 1)
})


test_that("failed_conditions names every enabled condition a model fails", {
  ref <- list(criterion = -1000, sse = 0.05)
  fc <- function(cand, abs_tol = NULL, rel_tol = 1e-3, sse_tol = 0.01, mode = "absolute") {
    failed_conditions(cand, ref, abs_tol, rel_tol, sse_tol, mode)
  }
  # criterion 2 % (20 units) better, error 0.02 lower -> passes everything
  expect_equal(fc(list(criterion = -1020, sse = 0.03), abs_tol = 10), character(0))
  # error drop 0.005 < 0.01
  expect_equal(fc(list(criterion = -1020, sse = 0.045)), "sse")
  # ... but 10 % of the remaining error > 5 % in relative mode
  expect_equal(fc(list(criterion = -1020, sse = 0.045), sse_tol = 0.05, mode = "relative"), character(0))
  # criterion worse -> fails both criterion conditions, not the error one
  expect_equal(fc(list(criterion = -990, sse = 0.01), abs_tol = 10), c("criterion_absolute", "criterion_relative"))
  # 5 units better: passes relative 0.1 % (0.5 %) but not absolute 10
  expect_equal(fc(list(criterion = -1005, sse = 0.03), abs_tol = 10), "criterion_absolute")
  # a NULL tolerance switches the condition off
  expect_equal(fc(list(criterion = -1020, sse = 0.045), sse_tol = NULL), character(0))
})


test_that("spectralem_select refuses to run with every condition switched off", {
  data <- six_peaks()
  expect_error(
    spectralem_select(data$x, data$y, criterion_relative_tolerance = NULL, sse_tolerance = NULL),
    "enable at least one condition"
  )
})


test_that("spectralem_select recovers the true number of peaks", {
  data <- six_peaks()
  res <- spectralem_select(data$x, data$y, criterion = "bic",
                           sse_tolerance = NULL, print_progress = FALSE)
  expect_equal(res$n_peaks, 6)
  expect_equal(sort(res$fit_params$pos), data$params$pos, tol = 1e-3)
  expect_equal(sum(res$selection$selected), 1)
  expect_true(all(c("n_peaks", "phase", "sse", "sse_raw", "bic", "aic", "compared_to", "failed",
                    "accepted", "selected") %in% names(res$selection)))
  # the conditions, not the cap, stopped the search, and the record says which
  expect_equal(res$decision$reason, "conditions")
  rejected_above <- res$selection[res$selection$n_peaks > res$n_peaks, ]
  expect_true(all(nzchar(rejected_above$failed)))
  expect_equal(sum(unlist(res$decision$limiting_conditions)) > 0, TRUE)
})


test_that("the strictest enabled condition decides", {
  data <- six_peaks()
  # An absolute criterion threshold no extra peak can reach stops the search
  # at the bottom of the window (or below), and is named as the reason.
  res <- spectralem_select(data$x, data$y, criterion_absolute_tolerance = 1e9,
                           criterion_relative_tolerance = NULL, sse_tolerance = NULL,
                           print_progress = FALSE, max_peaks = 10)
  expect_lt(res$n_peaks, 6)
  expect_equal(names(res$decision$limiting_conditions), "criterion_absolute")
})


test_that("spectralem_select searches downward from a coarse guess that is too high", {
  data <- six_peaks()
  # A curvature threshold this low counts noise ripples, so K0 is far too high.
  res <- spectralem_select(data$x, data$y, curvature_threshold = 1e-4, window = 1,
                           sse_tolerance = NULL, print_progress = FALSE, max_peaks = 12)
  expect_gt(res$coarse_estimate, 8)
  expect_true("down" %in% res$selection$phase || res$n_peaks < min(res$selection$n_peaks) + 2)
  expect_lte(res$n_peaks, 7)
})


test_that("spectralem_select respects max_peaks and refuses K", {
  data <- six_peaks()
  res <- spectralem_select(data$x, data$y, max_peaks = 4, print_progress = FALSE)
  expect_lte(max(res$selection$n_peaks), 4)
  expect_error(spectralem_select(data$x, data$y, K = 3), "chooses K itself")
})


test_that("patience accepts a later peak that beats the last accepted model", {
  # Five peaks with a sixth so close to the fifth that the pair needs two
  # Voigt profiles: the search must look past one weak addition.
  data <- six_peaks()
  res0 <- spectralem_select(data$x, data$y, patience = 0, sse_tolerance = NULL,
                            print_progress = FALSE)
  res3 <- spectralem_select(data$x, data$y, patience = 3, sse_tolerance = NULL,
                            print_progress = FALSE)
  expect_gte(res3$n_peaks, res0$n_peaks)
  # with patience, the search tries up to `patience` counts past the selected one
  expect_lte(max(res3$selection$n_peaks), res3$n_peaks + 4)
})


test_that("fitting peak counts in parallel batches changes nothing but the time", {
  data <- six_peaks()
  one <- spectralem_select(data$x, data$y, n_cores = 1, print_progress = FALSE)
  four <- spectralem_select(data$x, data$y, n_cores = 4, print_progress = FALSE)
  expect_identical(one$selection, four$selection)
  expect_identical(one$decision, four$decision)
  expect_equal(one$fit, four$fit)
})


test_that("parallel_map's socket-cluster path (used on Windows) matches lapply", {
  # Socket-cluster workers are fresh R processes. With the sources loaded into
  # a plain environment (as pyaxact does), that environment is serialized to
  # them whole; an installed package is instead loaded there by name. A
  # load_all() session is neither, so load the sources the pyaxact way.
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package sources not available")
  sem <- new.env(parent = globalenv())
  for (f in list.files(r_dir, pattern = "[.][Rr]$", full.names = TRUE)) sys.source(f, envir = sem)

  data <- six_peaks()
  fit_k <- sem$make_fitter(data$x, data$y, list(print_progress = FALSE))
  sequential <- lapply(5:6, fit_k)
  cluster <- sem$parallel_map(5:6, fit_k, 2, fork = FALSE)
  expect_equal(lapply(cluster, `[[`, "fit"), lapply(sequential, `[[`, "fit"))
})
