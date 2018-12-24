pro specfit_prepspec,wave0,spec0,errspec0,mask0,specpars,wave2,spec2,errspec2,mask2,object=object,$
                     vsini=vsini0,error=error,stp=stp,silent=silent,prepared=prepared,nonorm=nonorm,$
                     log=log,binsmooth=binsmooth0

;+
;
; SPECFIT_PREPSPEC
;
; This program preps spectra for comparison.  The spectrum will be
; smoothed to the proper resolution (unless /object is set),
; put on a new dispersion scale (with logarithmic spacing), and
; normalized.
;
; INPUTS:
;  wave       The wavelength array
;  spec       The spectrum array
;  errspec    The error array.  Set this to 0 if none exists (i.e. a
;               synthetic spectrum).
;  mask       The mask to use with wave/spec.  1-for points to use,
;               0-for points to ignore.  Set this to 0 to use all points.
;  specpars   A four element array of the spectral parameters: w0, w1, dw and res
;               w0  - The minimum wavelength in Angstroms.
;               w1  - The maximum wavelength in Angstroms.
;               dw  - The dispersion to use.  In log(wave) if /log is set.
;               res - The resolution in Angstroms (Gaussian FWHM).  If res is a scalar
;                       then the resolution will be a constant.  If res is an array then
;                       these are polynomial coefficients, i.e. resarr = poly(wave,res)
;                       If res<=0 then no smoothing is performed.
;  /object    This is an object spectrum, do NOT smooth.
;  /prepared  The spectrum has already been prepared, just trim if necessary.
;  =vsini     The rotational broadening (vsini) in km/s.
;  /nonorm    Don't normalize the spectrum (i.e. divide by the continuum)
;  =log       Logarithmic wavelength scale by default.  Use log=0 to
;               use linear wavelenth scale.
;  /binsmooth  Boxcar smooth (binsize=dispersion) the spectrum to simulate the binning
;               effects of a CCD observation.  The default is to perform this binning 
;               if /object is set and NOT perfrom this binning otherwise (i.e for synthetic
;               spectra).
;  /silent    Don't print anything to the screen
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  wave2      New wavelength array in Angstroms.
;  spec2      New spectrum array.
;  errspec2   New errpr array
;  mask2      New mask array.  If 0 was input for "mask" then mask2
;               will be all 1s
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>specfit_prepspec,wave,spec,errspec,mask,[4000.,4500.,1.0,3.7],wave2,spec2,errspec2,mask2
;
; By D.Nidever  Oct 2008
;-

undefine,error,wave,spec,errspec,mask,wave2,spec2,errspec2,mask2

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)
nmask = n_elements(mask0)
nspecpars = n_elements(specpars)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_prepspec,wave,spec,errspec,mask,specpars,wave2,spec2,errspec2,mask2,object=object'
  print,'                          prepared=prepared,vsini=vsini,nonorm=nonorm,log=log,binsmooth=binsmooth,'
  print,'                          error=error,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Logarithmic wavelength scale by default
if n_elements(log) eq 0 then log=1
if log[0] eq 0 then log=0 else log=1

; Internal arrays
wave = wave0
spec = spec0
errspec = errspec0
mask = mask0

; Wave and spec don't match
if nwave ne nspec then begin
  error = 'WAVE and SPEC must have same number of elements'
  print,error
  return
endif

; Errspec does not match SPEC/WAVE
; only exception is if errspec=0
if (nerrspec gt 1 and nerrspec ne nspec) or (nerrspec eq 1 and errspec[0] ne 0.0) then begin
  error = 'ERRSPEC must have the same number of elements as SPEC/WAVE or ERRSPEC=0'
  print,error
  return
endif

; Mask does not match SPEC/WAVE
; only exception is if errspec=0
if (nmask gt 1 and nmask ne nspec) or (nmask eq 1 and mask[0] ne 0.0) then begin
  error = 'MASK must have the same number of elements as SPEC/WAVE or MASK=0'
  print,error
  return
endif

; Mask=0 -> use all points
if (nmask eq 1 and mask[0] eq 0) then mask = wave*0.+1.0


; Not enough spectral parameters
if nspecpars lt 4 then begin
  error = 'SPECPARS must have at least 4 elements: W0, W1, DW, RES'
  print,error
  return
endif

; Spectral parameters
w0 = specpars[0]
w1 = specpars[1]
dw = specpars[2]
res = specpars[3:*]

; Vsini
if n_elements(vsini0) gt 0 then vsini=vsini0 else vsini=0.0

;mch = machar(/double)
;eps = mch.eps

; Out of range
if w0 lt (min(wave)-1e-7) then begin
  error = 'W0 lower than MIN(WAVE)'
  print,error
  return
endif
;if w1 gt (max(wave)+1d-7) then begin
;  error = 'W1 greater than MAX(WAVE)'
;  print,error
;  return
;endif


;------------------------
; PREPARING the spectrum
;------------------------
if not keyword_set(prepared) then begin

  ; Initialize temporary arrays
  undefine,ww,ss,ee,mm,ss2,ee2,mm2
  undefine,ss3,ee3,mm3,spec2a,errspec2a,mask2a


  ; MABYE IF THE SYNTHETIC SPECTRUM IS *VERY* HIGH RESOLUTION
  ; THEN IT SHOULD BE BINNED FIRST TO IMPROVE THE SPEED!!!!
  ; especially for specfit_setres and specfit_norm (DONE), vsini??

  ; MAYBE TRIM the spectrum if it is very long and we only need
  ; a small wavelength range.  Keep at least +/-5*res to help
  ; with the smoothing.

  ;-------------------------------------------------------
  ; Step 1: Make sure they are on LINEAR wavelength scale
  ;          For smoothing purposes
  ;-------------------------------------------------------
  disp = median(slope(wave))
  rdw = range(slope(wave))
  ; object spectra are NOT smoothed and therefore don't need to be on
  ;   linear scale
  if (rdw/disp gt 0.01) and not keyword_set(object) then begin   ;CHANGE THIS BACK!!!
  ;if (rdw/disp gt 0.01) then begin

    if not keyword_set(silent) then $
      print,'Spectrum NOT on linear wavelength scale.  Splining to linear wavelength scale'
    nww = range(wave)/disp
    nww = floor(nww)  ; can't go beyond end or spline will freak
    ww = findgen(nww)*disp+min(wave)
    ;ss = spline(wave,spec,ww)
    ss = cspline(wave,spec,ww)  ; cspline is ~25x faster

    ; Spline the error spectrum
    if nerrspec gt 1 then ee=cspline(wave,errspec,ww)

    ; Spline mask
    mm = cspline(wave,mask,ww)
    mm = float(mm gt 0.5)    ; if >0.5 then set to 1, else 0

  ; Already on a linear wavelength scale
  endif else begin
    ww = wave
    ss = spec
    if nerrspec gt 1 then ee = errspec
    mm = mask
  endelse


  ; ---Bin if RES_FINAL/DISP_ORIG and DISP_FINAL/DISP_ORIG is TOO HIGH---
  ;binflag = 1
  ;if (res/disp gt 10.0) and (10.d0^dw/disp gt 10.0) and (binflag eq 1) then begin
  ;
  ;  ww_orig = ww
  ;  ss_orig = ss
  ;  if nerrspec gt 1 then ee_orig=ee
  ;  mm_orig = mm
  ;
  ;  final_disp = (res/5.0) < (10.0^dw/5.0)
  ;  nbin = round(final_disp/disp)
  ;  nfinal = n_elements(ss)/nbin
  ;  ss = FREBIN(ss_orig,nfinal)
  ;  ww = FREBIN(ww_orig,nfinal)
  ;  if nerrspec gt 1 then ee=FREBIN(ee_orig,nfinal)
  ;  mm = FREbIN(mm_orig,nfinal)
  ;
  ;  ; ** THIS IS NOT GOOD ENOUGH!!!!!  IT CAUSES SHIFTS IN X AND Y.**
  ;  ; Could always use take every Nth point instead of binning.
  ;  ;  i.e. ss[0:*:3]
  ;
  ;  ; ** THIS IS VERY, VERY, VERY BAD!!! ***
  ;  ;ss = ss_orig[0:*:nbin]
  ;  ;ww = ww_orig[0:*:nbin]
  ;  ;if nerrspec gt 1 then ee = ee_orig[0:*:nbin]
  ;  ;mm = mm_orig[0:*:nbin]
  ;
  ;endif



  ;--------------------------------------
  ; Step 2: Rotational Broadening
  ;--------------------------------------
  if vsini[0] gt 0.0 then begin
    VSINI,ww,ss,vsini[0],ss2
  endif else begin
    ss2 = ss
  endelse
  if nerrspec gt 1 then ee2 = ee
  mm2 = mm


  ;--------------------------------------
  ; Step 3: Smooth to proper resolution
  ;--------------------------------------
  if not keyword_set(object) or res[0] le 0.0 then begin

    ; If the spectrum is VERY HIGH resolution, and the desired
    ; resolution is LOW, then maybe BIN first (maybe to 1/5 of the
    ;   new resolution).

    SPECFIT_SETRES,ww,ss2,res,ss3
    ; Smoothing reduces error.  Add gaussian-weighted errors in quadrature
    if nerrspec gt 1 then begin
      SPECFIT_SETRES,ww,ee2^2.0,res,ee3
      ee3 = sqrt(ee3)
    endif
  endif else begin
    if keyword_set(object) and not keyword_set(silent) then print,'Object spectrum - NO smoothing'
    if res[0] le 0.0 and not keyword_set(silent) then print,'Resolution <= 0.0 - NO smoothing'
    ss3 = ss2
    if nerrspec gt 1 then ee3 = ee2
  endelse
  mm3 = mm2      ; don't smooth the mask


  ;-------------------------------------
  ; Step 4: Put on new dispersion scale
  ;-------------------------------------
  if not keyword_set(object) then binsmooth=1 else binsmooth=0  ; only smooth syn spectra
  if keyword_set(binsmooth0) then binsmooth=1  ; /binsmooth set
  SPECFIT_SETDISP,ww,ss3,w0,w1,dw,wave2,spec2a,log=log,binsmooth=binsmooth
  undefine,dum,dum2
  if nerrspec gt 1 then SPECFIT_SETDISP,ww,ee3,w0,w1,dw,dum,errspec2a,log=log
  SPECFIT_SETDISP,ww,mm3,w0,w1,dw,dum2,mask2a,log=log
  mask2a = float(mask2a gt 0.5)    ; if >0.5 then set to 1, else 0



  ;--------------------
  ; Step 5: Normalize
  ;--------------------
  if nerrspec eq 1 then errspec2a=0
  if not keyword_set(nonorm) then begin
    SPECFIT_NORM,wave2,spec2a,errspec2a,spec2,cont
    if nerrspec gt 1 then errspec2=errspec2a/cont
  endif else begin
    spec2 = spec2a
    if nerrspec gt 1 then errspec2=errspec2a
  endelse
  mask2 = mask2a

  ;stop

  ; Remove temporary arrays
  undefine,ww,ss,ee,mm,ss2,ee2,mm2
  undefine,ss3,ee3,mm3,spec2a,errspec2a,mask2a



;-----------------------------
; ALREADY PREPARED, just TRIM
;-----------------------------
endif else begin

  ;------------------------------------------------
  ; Rotational Broadening
  ;   It doesn't matter if rotational broadening
  ;   is done before or after setres/setdisp/norm
  ;-------------------------------------------------
  if vsini[0] gt 0.0 then begin
    VSINI,wave,spec,vsini[0],temp
    spec = temp
    undefine,temp
  endif


  ; Check that the spectrum parameters match
  ; Log wavelength, same dispersion
  syndwlog = median(slope(alog10(wave)))
  diffdw = dw - syndwlog

  ; Log dispersion
  if range(slope(alog10(wave)))/median(slope(alog10(wave))) gt 1d-5 then begin
    error = 'Prepared spectrum DISPERSION is NOT sampled on LOG(WAVELENGTH) steps'
    print,error
    return
  end

  ; Comparing to input dispersion
  if diffdw/dw gt 1d-5 then begin
    error = 'Prepared spectrum DISPERSION '+strtrim(syndwlog,2)+$
            ' does NOT match the input dispersion '+strtrim(dw,2)
    print,error
    return
  endif

  ; Do we have W0
  w0closest = closest(w0,wave,ind=w0ind)
  if (w0closest-w0) gt 1d-5 then begin
    error = 'Prepared spectrum WAVELENGTH does NOT include W0 '+strtrim(w0,2)
    print,error
    return
  endif

  wave2 = wave
  spec2 = spec
  mask2 = mask
  if nerrspec gt 1 then errspec2 = errspec

  ; Trim if necessary
  if abs(min(wave2)-w0) gt 1d-5 or abs(max(wave2)-w1) gt 1d-5 then begin
    gdind = where(wave2 ge (w0-1d-5) and wave2 le (w1+1d-5),ngdind)
    wave2 = wave2[gdind]
    spec2 = spec2[gdind]
    mask2 = mask2[gdind]
    if nerrspec gt 1 then errspec2 = errspec2[gdind]
  endif

endelse

if keyword_set(stp) then stop

end
