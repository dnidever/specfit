pro specfit_setres,wave,spec,res,spec2,error=error,stp=stp

;+
;
; SPECFIT_SETRES
;
; This program sets the resolution of a spectrum.  This can be a constant
; or vary with wavelength.  The wavelength sampling/dispersion needs to be
; linear if this is to work properly.  This must be run before SPECFIT_SETDISP.PRO
;
; INPUTS:
;  wave    The wavelength array (in A).  Needs to be linear otherwise the
;            resolution will vary with wavelength.
;  spec    The spectrum array
;  res     The resolution in A (Gaussian FWHM).  If res is a scalar then the
;            resolution will be a constant.  If res is an array then these
;            are polynomial coefficients, i.e. resarr = poly(wave,res)
;  /stp    Stop at the end of the program
;
; OUTPUTS:
;  spec2   New spectrum array.
;  =error  The error message if there was one.
;
; USAGE:
;  IDL>specfit_setres,wave,spec,3.0,spec2
;
; By D.Nidever  Oct 2008
;-

undefine,error,spec2

nwave = n_elements(wave)
nspec = n_elements(spec)
nres = n_elements(res)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nres eq 0 then begin
  print,'Syntax - specfit_setres,wave,spec,res,spec2,error=error,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Wave and spec don't match
if nwave ne nspec then begin
  print,'WAVE and SPEC must have same number of elements'
  error = 'WAVE and SPEC must have same number of elements'
  return
endif

; Convert res in A to pix
dw = median(slope(wave))
rdw = range(slope(wave))

; Not linear wavelength sampling
if rdw/dw gt 0.01 then begin
  print,'NOT linear wavelength sampling'
  error = 'NOT linear wavelength sampling'
  return
endif


; Constant resolution
if nres eq 1 then begin

  ; Convert res in A to pix
  respix = res/dw

  ; Use GSMOOTH
  ;spec2 = GSMOOTH(spec,respix)
  spec2 = GSMOOTH(spec,respix,widfwhm=3.0)


; Varying resolution
endif else begin

  resarr = POLY(wave,res)

  ; Don't know how to do this yet
  print,'CANNOT DO VARYING RESOLUTION YET'
  stop

endelse

if keyword_set(stp) then stop

end
