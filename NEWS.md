# spectralem 0.10.0

* New `spectralem_select()`: chooses the number of peaks automatically, following
  Kasterke et al., "An expectation-maximization algorithm for spectral
  reconstruction under the spectral hard model", Chemometr. Intell. Lab. Syst.
  267 (2025) 105518, Section 5. A coarse count from the curvature of a smoothing
  spline (`estimate_peak_count()`, also exported) seeds a window of fits; the
  search then adds peaks while each is a significant improvement in BIC (or AIC)
  and in the sum of squared errors of the max-normalized spectrum.
* Extensions beyond the paper: AIC as an alternative criterion, a downward search
  when the coarse count is too high, a `patience` look-ahead over successive
  insignificant peak additions, an absolute or relative error tolerance, and
  parallel window fits (`n_cores`).

# spectralem 0.9

* Almost CRAN ready

# spectralem 0.1

* Added a `NEWS.md` file to track changes to the package.
