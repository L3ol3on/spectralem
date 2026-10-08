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

# Is model `cand` (more peaks) a significant improvement over `ref`?
# Both the criterion and the SSE must improve by more than their tolerance
# (the paper stops on EITHER a negligible criterion gain OR an insignificant
# SSE drop, so continuing needs both).
is_significant_improvement <- function(cand, ref, criterion_tolerance,
                                       sse_tolerance, sse_tolerance_mode) {
  crit_gain <- (ref$criterion - cand$criterion) / max(abs(ref$criterion), .Machine$double.eps)
  sse_gain <- ref$sse - cand$sse
  if (sse_tolerance_mode == "relative") {
    sse_gain <- sse_gain / max(ref$sse, .Machine$double.xmin)
  }
  crit_gain > criterion_tolerance && sse_gain > sse_tolerance
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
#'         within [\code{min_peaks}, \code{max_peaks}]) and take the one with
#'         the lowest criterion.
#'   \item If that model is at the upper edge of the window, keep adding one
#'         peak at a time. A peak count is accepted only if it improves
#'         significantly on the last accepted model (see below). After an
#'         insignificant one, up to \code{patience} further peak counts are
#'         tried; the first significant one is accepted and the search goes
#'         on, otherwise the last accepted model is selected.
#'   \item If it is at the lower edge, remove one peak at a time as long as
#'         the larger model is NOT a significant improvement on the smaller
#'         one -- the mirror image of the upward rule.
#' }
#' Model \eqn{B} is a significant improvement on model \eqn{A} (fewer peaks)
#' when both
#' \deqn{(C_A - C_B) / |C_A| > \code{criterion\_tolerance}}
#' and the drop in the sum of squared errors, \eqn{S_A - S_B} (absolute) or
#' \eqn{(S_A - S_B) / S_A} (relative), exceeds \code{sse_tolerance}. Here
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
#' @param criterion_tolerance minimum relative improvement of the criterion
#' @param sse_tolerance minimum drop of the sum of squared errors
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
#'         \code{bic}, \code{aic}, \code{accepted} (whether the search moved
#'         to this peak count; \code{NA} for window fits) and \code{selected}
#' }
#' @export
spectralem_select <- function(x, y,
                              criterion = c("bic", "aic"),
                              criterion_tolerance = 1e-3,
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
      accepted = NA, selected = FALSE
    )
  }
  row_of <- function(k) {
    tab <- do.call(rbind, rows)
    r <- tab[tab$n_peaks == k, ]
    list(sse = r$sse, criterion = r[[criterion]])
  }
  mark <- function(k, accepted) {
    for (i in seq_along(rows)) {
      if (rows[[i]]$n_peaks == k) rows[[i]]$accepted <<- accepted
    }
  }
  significant <- function(more, fewer) {
    is_significant_improvement(row_of(more), row_of(fewer),
                               criterion_tolerance, sse_tolerance, sse_tolerance_mode)
  }

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
  tab <- do.call(rbind, rows)
  best <- tab$n_peaks[which.min(tab[[criterion]])]

  # 2) Upward: accept only significant improvements over the last accepted
  #    model, looking `patience` peak counts past an insignificant one.
  if (best == hi) {
    k <- best
    misses <- 0
    while (k < max_peaks && misses <= patience) {
      k <- k + 1
      record(k, fit_k(k), "up")
      if (significant(k, best)) {
        mark(k, TRUE)
        best <- k
        misses <- 0
      } else {
        mark(k, FALSE)
        misses <- misses + 1
      }
    }
  }

  # 3) Downward: drop a peak while the larger model does not earn it.
  if (best == lo) {
    k <- best
    while (k > min_peaks) {
      record(k - 1, fit_k(k - 1), "down")
      keep_larger <- significant(k, k - 1)
      mark(k - 1, !keep_larger)
      if (keep_larger) break
      k <- k - 1
      best <- k
    }
  }

  selection <- do.call(rbind, rows)
  selection$selected <- selection$n_peaks == best
  selection <- selection[order(selection$n_peaks), ]
  rownames(selection) <- NULL

  out <- fits[[as.character(best)]]
  out$n_peaks <- best
  out$coarse_estimate <- k0
  out$selection <- selection
  out
}
