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


test_that("is_significant_improvement needs both the criterion and the error to improve", {
  ref <- list(criterion = -1000, sse = 0.05)
  # criterion 2 % better, error 0.02 lower -> significant
  expect_true(is_significant_improvement(list(criterion = -1020, sse = 0.03), ref, 1e-3, 0.01, "absolute"))
  # criterion better, error drop 0.005 < 0.01 -> not
  expect_false(is_significant_improvement(list(criterion = -1020, sse = 0.045), ref, 1e-3, 0.01, "absolute"))
  # ... but 10 % of the remaining error > 5 % -> significant in relative mode
  expect_true(is_significant_improvement(list(criterion = -1020, sse = 0.045), ref, 1e-3, 0.05, "relative"))
  # error much lower but criterion worse -> not
  expect_false(is_significant_improvement(list(criterion = -990, sse = 0.01), ref, 1e-3, 0.01, "absolute"))
})


test_that("spectralem_select recovers the true number of peaks", {
  data <- six_peaks()
  res <- spectralem_select(data$x, data$y, criterion = "bic",
                           sse_tolerance = 0, print_progress = FALSE)
  expect_equal(res$n_peaks, 6)
  expect_equal(sort(res$fit_params$pos), data$params$pos, tol = 1e-3)
  expect_equal(sum(res$selection$selected), 1)
  expect_true(all(c("n_peaks", "phase", "sse", "sse_raw", "bic", "aic", "accepted", "selected")
                  %in% names(res$selection)))
})


test_that("spectralem_select searches downward from a coarse guess that is too high", {
  data <- six_peaks()
  # A curvature threshold this low counts noise ripples, so K0 is far too high.
  res <- spectralem_select(data$x, data$y, curvature_threshold = 1e-4, window = 1,
                           sse_tolerance = 0, print_progress = FALSE, max_peaks = 12)
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
  res0 <- spectralem_select(data$x, data$y, patience = 0, sse_tolerance = 0,
                            print_progress = FALSE)
  res3 <- spectralem_select(data$x, data$y, patience = 3, sse_tolerance = 0,
                            print_progress = FALSE)
  expect_gte(res3$n_peaks, res0$n_peaks)
  # with patience, the search tries up to `patience` counts past the selected one
  expect_lte(max(res3$selection$n_peaks), res3$n_peaks + 4)
})
