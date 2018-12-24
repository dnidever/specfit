pro specfit_setdisp,wave,spec,w0,w1,dw,wave2,spec2,log=log,binsmooth=binsmooth,error=error,stp=stp

;+
;
; SPECFIT_SETDISP
;
; This program sets the dispersion of a spectrum.  This must be run AFTER
; SPECFIT_SETRES.PRO.
;
; INPUTS:
;  wave        The wavelength array
;  spec        The spectrum array
;  w0          The minimum wavelength in Angstroms.
;  w1          The maximum wavelength in Angstroms.  MAX(WAVE) will
;                probably be slightly different than W1 because the exact
;                values are used for W1 and DW.
;  dw          The dispersion to use in Angstroms.  If /log then this
;                should be the DW in log Angstroms.
;  /log        Use a constant log(wave) sampling.  W0 and W1 should still
;                be in Angstroms.
;  /binsmooth  Use boxcar smoothing to reproduce the smoothing effects of binning.
;                This normally has a minimal effect compared to the
;                "resolution" smoothing by specfit_setres.pro
;  /stp        Stop at the end of the program
;
; OUTPUTS:
;  wave2       New wavelength array.  This will always be in A (i.e. linear).
;  spec2       New spectrum array.
;  =error      The error message if there was one.
;
; USAGE:
;  IDL>specfit_setdisp,wave,spec,4000.,4500.,1.0,wave2,spec2
;
; By D.Nidever  Oct 2008
;-

; This doesn't check that WAVE is actually within W0/W1.
; Can cause splining problems at the end.


undefine,error,wave2,spec2

nwave = n_elements(wave)
nspec = n_elements(spec)
nw0 = n_elements(w0)
nw1 = n_elements(w1)
ndw = n_elements(dw)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nw0 eq 0 or nw1 eq 0 or ndw eq 0 then begin
  print,'Syntax - specfit_setdisp,wave,spec,w0,w1,dw,wave2,spec2,log=log,binsmooth=binsmooth,error=error,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Wave and spec don't match
if nwave ne nspec then begin
  print,'WAVE and SPEC must have same number of elements'
  error = 'WAVE and SPEC must have same number of elements'
  return
endif

dwcur = median(slope(wave))  ;current dispersion

; Make new wavelength array
;--------------------------

; Linear scale
if not keyword_set(log) then begin

  nw = (w1-w0)/dw + 1.0d0
  if abs(nw-(long(nw)+1L)) lt 1d-6 then nw = long(nw)+1L   ; round up if VERY close
  nw = long(floor(nw))

  ;dw1 = (double(w1)-double(w0))/(nw-1.0d0)
  ;wave2 = dindgen(nw)*dw1 + double(w0)

  wave2 = dindgen(nw)*double(dw) + double(w0)

  ; Do we need to spline
  ;maxdiff = max(abs(wave2-wave)/wave)
  ;if maxdiff lt 1d-10 then dospline=0 else dospline=1

  ; Smoothing length in pixels
  ;smlen = round(dw1/dwcur)
  smlen = round(dw/dwcur)

; Log scale
endif else begin

  ; W0 and W1 are in Angstroms, DW in log(Angstroms)

  ; w1 = 10.^((nw-1)*dwlog+w0log)
  ; nw = ( log(w1)-log(w0) ) / dwlog
  nw = ( alog10(double(w1)) - alog10(double(w0)) ) / double(dw) + 1.0d0
  if abs(nw-(long(nw)+1L)) lt 1d-8 then nw = long(nw)+1L   ; round up if VERY close
  nw = long(floor(nw))

  wave2 = (10.0d0)^( dindgen(nw)*double(dw) + alog10(double(w0)) )

  ;nw = floor((w1-w0)/dw)
  ;; w1 = 10.^((nw-1)*dwlog+w0log)
  ;dwlog = (alog10(double(w1))-alog10(double(w0)))/(nw-1.0d0)
  ;wave2 = (10.0d0)^(dindgen(nw)*dwlog+alog10(double(w0)))

  ; Do we need to spline
  ;maxdiff = max(abs(wave2-wave)/wave)
  ;if maxdiff lt 1d-10 then dospline=0 else dospline=1

  ; Smoothing length in pixels
  dwlinear = 10.0d0^dw
  smlen = round(dwlinear/dwcur)

endelse
nwave2 = n_elements(wave2)

;; Already on the correct dispersion scale
;if nwave eq nwave2 then begin
;  diffwave = abs(wave2-wave)
;  if max(diffwave)/median(wave) lt 1e-7 then begin
;    wave2 = wave
;    spec2 = spec
;  endif
;endif

; Boxcar smoothing the spectrum to "simulate" binning
;-------------------------------------------------------
if keyword_set(binsmooth) and smlen ge 1.0 then begin
  spec1 = SMOOTH(spec,smlen,/edge_truncate,/nan)
endif else begin
  spec1 = spec
endelse

; Spline to the desired dispersion scale
;-----------------------------------------
;spec2 = spline(wave,spec,wave2)
spec2 = cspline(wave,spec1,wave2)  ; cspline is ~25x faster

if keyword_set(stp) then stop

end
