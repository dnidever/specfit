# specfit
Spectral Fitting Program

Specfit is an automated spectral fitting program using a library of synthetic spectra.  It's very
generic but does requires fully reducted and wavelength calibrated 1-D spectra with uncertainties.

Here are the inputs, output and usage for the main specfit.pro program:

```
; INPUTS:
;  wave       The observed wavelength array
;  spec       The observed spectrum array
;  errspec    The error spectrum array
;  synstr     The synthetic grid structure.  It should have TEFF, LOGG, METAL, ALPHA and FILE
;               tags.  TEFF, LOGG, METAL and ALPHA should be 1D arrays of the UNIQUE stellar
;               parameters of the synthetic grid.  FILE should be a 4D string array (or 3D if
;               ALPHA only has one element) of dimensions [Nteff,Nlogg,Nmetal,Nalpha] with
;               the absolute path of the synthetic spectrum FITS file.  It should be an empty
;               string if the spectrum does not exist.
;  specpars   The spectrum parameters, specpars = [w0, w1, dw, res]
;  =mask      The mask to use with wave/spec.  1-for points to use,
;               0-for points to ignore.
;               This can also be an array of weights.
;  /zerovel   Do not solve for velocity, should already be set to rest frame.
;  /fitvsini  Fit vsini.
;  =fixvsini  Fix vsini at a certain value (in km/s).
;  /monte     Get errors via Monte Carlo simulation.  This takes 50x longer.
;  =nmonte    The number of Monte Carlo simulations (minimum is 10).
;  /plot      Plot the spectra.
;  /silent    Don't print anything to the screen
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  fitstr     Structure that contains all of the output information, including:
;     bestpars   The array of best-fitting parameters, Teff, logg, metal and alpha.
;     chisq      The chi squared of the best fit.
;     vrel       The relative velocity of the best fit
;  =fitstr_monte  The structure of Monte Carlo best fits.
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>specfit,wave,spec,errspec,synstr,[4000.,4500.,1.4,2.7],fitstr,dir=specdir
```