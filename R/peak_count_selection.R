# Automatic selection of the number of peaks.
#
# Implements the model-selection procedure of
#   M. Kasterke, L. Kaufmann, M. Kateri, T. Brands, "An expectation-maximization
#   algorithm for spectral reconstruction under the spectral hard model",
#   Chemometr. Intell. Lab. Syst. 267 (2025) 105518, Section 5,
# extended by (a) AIC as an alternative to BIC, (b) a downward search when the
# coarse guess is too high, and (c) a look-ahead ("patience") over successive
# insignificant peak additions. Choices the paper leaves open are marked
# ASSUMPTION below.

#' Coarse estimate of the number of peaks in a spectrum
#'
#' @description
#' Fits a smoothing spline to the max-normalized spectrum and counts the
#' local minima of its second derivative that lie below zero (negative
#' curvature, i.e. a peak or a shoulder), as in Kasterke et al. (2025),
#' Section 5. Shoulders count as candidates too, since an asymmetric band
#' needs several Voigt profiles.
#'
#' ASSUMPTION: the paper does not say how noise ripples are told apart from
#' peaks. Here a minimum counts only if its depth is at least
#' \code{curvature_threshold} times that of the deepest one.
#'
#' @param x the \code{x} coordinates of the signal
#' @param y the function values of the signal at \code{x}
#' @param spar smoothing parameter of \code{stats::smooth.spline}; \code{NULL}
#'   chooses it by generalized cross-validation
#' @param curvature_threshold minimum depth of a second-derivative minimum,
#'   as a fraction of the deepest one, for it to count as a peak
#' @return the estimated number of peaks, at least 1
#' @export
estimate_peak_count <- function(x, y, spar = NULL, curvature_threshold = 0.05) {
  o <- order(x)
  x <- x[o]
  y <- y[o] / max(abs(y))

  fit <- stats::smooth.spline(x, y, spar = spar)
  d2 <- stats::predict(fit, x, deriv = 2)$y

  n <- length(d2)
  if (n < 3) {
    return(1L)
  }
  inner <- 2:(n - 1)
  is_min <- c(FALSE, d2[inner] < d2[inner - 1] & d2[inner] <= d2[inner + 1], FALSE)
  depth <- -d2[is_min & d2 < 0]
  # Curvature made dimensionless by the x range: on a featureless signal the
  # deepest "minimum" is rounding noise, and a relative threshold alone would
  # count every ripple of it.
  if (length(depth) == 0 || max(depth) * diff(range(x))^2 < 1e-6) {
    return(1L)
  }
  max(1L, sum(depth >= curvature_threshold * max(depth)))
}

# Information criterion of a fit with `p` free parameters, `sse` measured on
# the max-normalized spectrum over `h` points. BIC as in the paper; AIC
# replaces the log(h) penalty per parameter by 2.
information_criterion <- function(sse, h, p, criterion) {
  penalty <- if (criterion == "bic") log(h) else 2
  h * log(max(sse, .Machine$double.xmin) / h) + penalty * p
}

# Which of the enabled conditions does model `cand` (more peaks) fail against
# `ref`? A condition whose tolerance is NULL is disabled. An empty result means
# `cand` is a significant improvement; since every enabled condition must
# hold, the strictest one decides.
#   criterion_absolute:  C_ref - C_cand                > criterion_absolute_tolerance
#   criterion_relative: (C_ref - C_cand) / |C_ref|     > criterion_relative_tolerance
#   sse (absolute):      S_ref - S_cand                > sse_tolerance
#   sse (relative):     (S_ref - S_cand) / S_ref       > sse_tolerance
failed_conditions <- function(cand, ref, criterion_absolute_tolerance,
                              criterion_relative_tolerance, sse_tolerance,
                              sse_tolerance_mode) {
  crit_gain <- ref$criterion - cand$criterion
  failed <- character(0)
  if (!is.null(criterion_absolute_tolerance) && !(crit_gain > criterion_absolute_tolerance)) {
    failed <- c(failed, "criterion_absolute")
  }
  rel_gain <- crit_gain / max(abs(ref$criterion), .Machine$double.eps)
  if (!is.null(criterion_relative_tolerance) && !(rel_gain > criterion_relative_tolerance)) {
    failed <- c(failed, "criterion_relative")
  }
  sse_gain <- ref$sse - cand$sse
  if (sse_tolerance_mode == "relative") {
    sse_gain <- sse_gain / max(ref$sse, .Machine$double.xmin)
  }
  if (!is.null(sse_tolerance) && !(sse_gain > sse_tolerance)) {
    failed <- c(failed, "sse")
  }
  failed
}

#' Fit Voigt peaks with an automatically selected number of peaks
#'
#' @description
#' Runs \code{spectralem()} for several peak counts and picks one by an
#' information criterion (Kasterke et al. 2025, Section 5):
#' \enumerate{
#'   \item Normalize the spectrum by its maximum.
#'   \item Estimate a coarse peak count \eqn{K_0} with \code{estimate_peak_count()}.
#'   \item Fit every \eqn{K} in \eqn{K_0 \pm} \code{window} (shifted to stay
#'         within [\code{min_peaks}, \code{max_peaks}]), in parallel with
#'         \code{n_cores} -- the paper's "minimum of five spectral analyses".
#'   \item Starting from the smallest of them, step up one peak at a time
#'         (reusing the window's fits, then fitting further counts). A count
#'         is accepted only if it is a significant improvement on the last
#'         accepted model (see below). After one that is not, up to
#'         \code{patience} further counts are tried; the first significant
#'         one is accepted and the search goes on, otherwise the last
#'         accepted model is selected. (The paper instead takes the window's
#'         lowest criterion and extends only from its upper edge; walking the
#'         window with the conditions means the enabled conditions always
#'         decide.)
#'   \item If nothing above the smallest window count was accepted, remove
#'         one peak at a time as long as the larger model is NOT a
#'         significant improvement on the smaller one -- the mirror image of
#'         the upward rule.
#' }
#' Model \eqn{B} is a significant improvement on model \eqn{A} (fewer peaks)
#' when it passes every ENABLED condition (a \code{NULL} tolerance disables
#' one), so the strictest enabled condition decides:
#' \itemize{
#'   \item absolute criterion gain: \eqn{C_A - C_B >} \code{criterion_absolute_tolerance}
#'   \item relative criterion gain: \eqn{(C_A - C_B) / |C_A| >} \code{criterion_relative_tolerance}
#'   \item error drop: \eqn{S_A - S_B} (absolute) or \eqn{(S_A - S_B) / S_A}
#'         (relative) \eqn{>} \code{sse_tolerance}
#' }
#' Here
#' \eqn{C} is BIC \eqn{= h \log(S/h) + p \log h} or AIC
#' \eqn{= h \log(S/h) + 2p}, \eqn{S} the sum of squared errors of the
#' max-normalized spectrum, \eqn{h} the number of points and \eqn{p = 4K}
#' plus the number of background parameters.
#'
#' Each peak count is an independent \code{spectralem()} run from scratch.
#'
#' @param x the \code{x} coordinates of the signal
#' @param y the function values of the signal at \code{x}
#' @param criterion \code{"bic"} or \code{"aic"}. AIC penalizes each extra
#'   parameter less, so it tends to select more peaks.
#' @param criterion_absolute_tolerance minimum improvement of the criterion,
#'   in criterion units; \code{NULL} disables this condition
#' @param criterion_relative_tolerance minimum improvement of the criterion as
#'   a fraction of its magnitude; \code{NULL} disables this condition
#' @param sse_tolerance minimum drop of the sum of squared errors;
#'   \code{NULL} disables this condition
#' @param sse_tolerance_mode \code{"absolute"} (drop in normalized SSE) or
#'   \code{"relative"} (drop as a fraction of the previous SSE)
#' @param patience number of further peak counts tried after an
#'   insignificant one before the search stops; 0 stops at the first
#' @param window half-width of the initial window around \eqn{K_0}; the
#'   paper's "minimum of five spectral analyses" is \code{window = 2}
#' @param min_peaks,max_peaks range of peak counts that may be tested
#' @param spar,curvature_threshold passed to \code{estimate_peak_count()}
#' @param n_cores fit the initial window's peak counts in parallel with
#'   \code{parallel::mclapply} (forking; ignored on Windows)
#' @param ... further arguments to \code{spectralem()} (not \code{K})
#'
#' @return The \code{spectralem()} result of the selected model, plus
#' \itemize{
#'   \item n_peaks - the selected number of peaks
#'   \item coarse_estimate - \eqn{K_0}
#'   \item selection - a data frame with one row per tested peak count:
#'         \code{n_peaks}, \code{phase} (window / up / down), \code{sse}
#'         (normalized), \code{sse_raw} (in the units of \code{y}),
#'         \code{bic}, \code{aic}, \code{compared_to} (the other model in
#'         the comparison that decided this row), \code{failed} (the
#'         conditions the model with MORE peaks failed in that comparison,
#'         comma separated; empty if it passed all; \code{NA} for window
#'         fits), \code{accepted} (whether the search moved to this peak
#'         count) and \code{selected}. Upward rows are the larger model;
#'         downward rows the smaller one.
#'   \item decision - why the search stopped at \code{n_peaks}: a list with
#'         \code{reason} (\code{"conditions"}, \code{"max_peaks"} or
#'         \code{"min_peaks"}) and
#'         \code{limiting_conditions}, a named count of the conditions failed by
#'         the rejected models with more peaks than the selected one -- the
#'         condition(s) that limited the number of peaks
#' }
#' @export
spectralem_select <- function(x, y,
                              criterion = c("bic", "aic"),
                              criterion_absolute_tolerance = NULL,
                              criterion_relative_tolerance = 1e-3,
                              sse_tolerance = 0.01,
                              sse_tolerance_mode = c("absolute", "relative"),
                              patience = 3,
                              window = 2,
                              min_peaks = 1,
                              max_peaks = 60,
                              spar = NULL,
                              curvature_threshold = 0.05,
                              n_cores = 1,
                              ...) {
  criterion <- match.arg(criterion)
  sse_tolerance_mode <- match.arg(sse_tolerance_mode)
  if (is.null(criterion_absolute_tolerance) && is.null(criterion_relative_tolerance) &&
      is.null(sse_tolerance)) {
    stop("enable at least one condition (criterion_absolute_tolerance, ",
         "criterion_relative_tolerance or sse_tolerance); with none, every extra peak is accepted")
  }
  args <- list(...)
  if ("K" %in% names(args)) {
    stop("spectralem_select() chooses K itself; do not pass it")
  }
  start_peaks <- args$start_peaks
  min_peaks <- max(1L, as.integer(min_peaks), length(start_peaks$pos))
  max_peaks <- as.integer(max_peaks)
  if (max_peaks < min_peaks) {
    stop("max_peaks (", max_peaks, ") is below min_peaks (", min_peaks, ")")
  }

  h <- length(x)
  scale <- max(abs(y))
  background <- if (is.null(args$background_model)) list(linear = TRUE) else args$background_model
  n_background <- 2L * isTRUE(background$linear) +
    sum(vapply(background, function(b) mode(b) == "numeric", logical(1)))

  fits <- list()
  rows <- list()
  fit_k <- function(k) do.call(spectralem, c(list(x = x, y = y, K = k), args))
  record <- function(k, res, phase) {
    sse_raw <- sum((y - res$fit)^2)
    sse <- sse_raw / scale^2
    p <- 4L * k + n_background
    fits[[as.character(k)]] <<- res
    rows[[length(rows) + 1]] <<- data.frame(
      n_peaks = k, phase = phase, sse = sse, sse_raw = sse_raw,
      bic = information_criterion(sse, h, p, "bic"),
      aic = information_criterion(sse, h, p, "aic"),
      compared_to = NA_integer_, failed = NA_character_,
      accepted = NA, selected = FALSE
    )
  }
  row_of <- function(k) {
    tab <- do.call(rbind, rows)
    r <- tab[tab$n_peaks == k, ]
    list(sse = r$sse, criterion = r[[criterion]])
  }
  mark <- function(k, accepted, compared_to, failed) {
    for (i in seq_along(rows)) {
      if (rows[[i]]$n_peaks == k) {
        rows[[i]]$accepted <<- accepted
        rows[[i]]$compared_to <<- as.integer(compared_to)
        rows[[i]]$failed <<- paste(failed, collapse = ",")
      }
    }
  }
  failures <- function(more, fewer) {
    failed_conditions(row_of(more), row_of(fewer), criterion_absolute_tolerance,
                      criterion_relative_tolerance, sse_tolerance, sse_tolerance_mode)
  }
  # Models with more peaks than the one finally selected that the conditions
  # rejected, with what they failed: the record of what limited the count.
  rejected <- list()
  reject <- function(k, failed) rejected[[length(rejected) + 1]] <<- list(n_peaks = k, failed = failed)

  # 1) Coarse estimate and initial window.
  k0 <- estimate_peak_count(x, y, spar = spar, curvature_threshold = curvature_threshold)
  k0 <- min(max(k0, min_peaks), max_peaks)
  lo <- max(min_peaks, k0 - window)
  hi <- min(max_peaks, lo + 2 * window)
  lo <- max(min_peaks, hi - 2 * window)
  ks <- lo:hi
  window_fits <- if (n_cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(ks, fit_k, mc.cores = n_cores)
  } else {
    lapply(ks, fit_k)
  }
  for (i in seq_along(ks)) {
    if (inherits(window_fits[[i]], "try-error")) stop(window_fits[[i]])
    record(ks[i], window_fits[[i]], "window")
  }
  # 2) Walk up from the bottom of the window -- through the counts already
  #    fitted, then beyond -- accepting a count only if it passes every
  #    enabled condition against the last accepted one, and looking
  #    `patience` counts past one that does not. (The paper takes the window's
  #    lowest criterion instead; walking it with the conditions means the
  #    enabled conditions always decide, wherever the best model lies.)
  best <- lo
  k <- lo
  misses <- 0
  while (k < max_peaks && misses <= patience) {
    k <- k + 1
    if (is.null(fits[[as.character(k)]])) record(k, fit_k(k), "up")
    failed <- failures(k, best)
    mark(k, length(failed) == 0, best, failed)
    if (length(failed) == 0) {
      best <- k
      misses <- 0
    } else {
      reject(k, failed)
      misses <- misses + 1
    }
  }
  # Ran out of room before `patience` was used up: the cap, not the
  # conditions, ended the search.
  reason <- if (misses <= patience) "max_peaks" else "conditions"

  # 3) Downward: nothing above the bottom of the window was worth it, so
  #    check the bottom itself: drop a peak while the larger model does not
  #    earn it.
  if (best == lo) {
    k <- best
    if (k <= min_peaks && reason == "conditions") reason <- "min_peaks"
    while (k > min_peaks) {
      record(k - 1, fit_k(k - 1), "down")
      failed <- failures(k, k - 1)
      # this row is the smaller model; it is moved to when the larger one fails
      mark(k - 1, length(failed) > 0, k, failed)
      if (length(failed) == 0) break
      reject(k, failed)
      reason <- "conditions"
      k <- k - 1
      best <- k
      if (k <= min_peaks) reason <- "min_peaks"
    }
  }

  limiting <- unlist(lapply(rejected, function(r) if (r$n_peaks > best) r$failed))
  decision <- list(
    reason = reason,
    limiting_conditions = as.list(table(factor(limiting, levels = unique(limiting))))
  )

  selection <- do.call(rbind, rows)
  selection$selected <- selection$n_peaks == best
  selection <- selection[order(selection$n_peaks), ]
  rownames(selection) <- NULL

  out <- fits[[as.character(best)]]
  out$n_peaks <- best
  out$coarse_estimate <- k0
  out$selection <- selection
  out$decision <- decision
  out
}
