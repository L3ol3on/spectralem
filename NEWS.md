# spectralem 0.10.0

* New `spectralem_select()`: chooses the number of peaks automatically, following
  Kasterke et al., "An expectation-maximization algorithm for spectral
  reconstruction under the spectral hard model", Chemometr. Intell. Lab. Syst.
  267 (2025) 105518, Section 5. A coarse count from the curvature of a smoothing
  spline (`estimate_peak_count()`, also exported) seeds a window of fits; the
  search then adds peaks while each passes every enabled condition: an absolute
  and a relative improvement in BIC (or AIC), and a drop in the sum of squared
  errors of the max-normalized spectrum. `decision` records which conditions
  limited the number of peaks.
* Extensions beyond the paper: AIC as an alternative criterion, switchable
  conditions (the paper combines a relative criterion gain with an absolute error
  drop), the window walked with the conditions instead of taking its lowest
  criterion, a downward search when the coarse count is too high, a `patience`
  look-ahead over successive insignificant peak additions, and parallel window
  fits (`n_cores`).

# spectralem 0.9

* Almost CRAN ready

# spectralem 0.1

* Added a `NEWS.md` file to track changes to the package.
