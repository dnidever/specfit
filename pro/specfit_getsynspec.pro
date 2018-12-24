pro specfit_getsynspec,synstr,temp0,logg0,metal0,alpha0,wave,spec,interp=interp,nearest=nearest,$
                       head=head,error=error,silent=silent,stp=stp

;+
;
; SPECFIT_GETSYNSPEC
;
; This program gets a synthetic spectrum for the desired stellar
; parameters.
;
; If the desired stellar parameters do not fall on a grid point then
; the spectrum can either be interpolated (/interp) or the nearest grid point
; can be used (/nearest).
;
; **NOTE** The stellar parameters will be rounded to the nearest 0.001.
;
; INPUTS:
;  synstr    The synthetic grid structure.  It should have TEFF, LOGG, METAL, ALPHA and FILE
;              tags.  TEFF, LOGG, METAL and ALPHA should be 1D arrays of the UNIQUE stellar
;              parameters of the synthetic grid.  FILE should be a 4D string array (or 3D if
;              ALPHA only has one element) of dimensions [Nteff,Nlogg,Nmetal,Nalpha] with
;              the absolute path of the synthetic spectrum FITS file.  It should be an empty
;              string if the spectrum does not exist.
;  teff      The effective temperature in Kelvin.
;  logg      The surface gravity.
;  metal     The metallicity, [M/H].
;  alpha     The alpha abundance, [alpha/Fe]
;  /interp . Interpolate between grid points if necessary.
;  /nearest  If the stellar parameters are between grid points then
;              use the closest one instead of interpolating.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  wave    The wavelength array
;  spec    The spectrum array
;  =head   The header of the synthetic spectrum
;  =error  The error message if there was one.
;
; USAGE:
;  IDL>specfit_getsynspec,synstr,5050.,2.50,-1.5,0.0,wave,spec
;
; By D.Nidever  Oct 2008
;-

undefine,error,wave,spec,head

nsynstr = n_elements(synstr)
ntemp = n_elements(temp0)
nlogg = n_elements(logg0)
nmetal = n_elements(metal0)
nalpha = n_elements(alpha0)

; Not enough inputs
if ntemp eq 0 or nlogg eq 0 or nmetal eq 0 or nalpha eq 0 then begin
  print,'Syntax - specfit_getsynspec,synstr,teff,logg,metal,alpha,wave,spec,interp=interp,nearest=nearest,error=error,'
  print,'                            head=head,silent=silent,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Checking the stellar parameters
if (VALID_NUM(temp0,ftemp) eq 0) then begin
  error = temp+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(logg0,flogg) eq 0) then begin
  error = logg+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(metal0,fmetal) eq 0) then begin
  error = metal+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(alpha0,falpha) eq 0) then begin
  error = alpha+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif

; Local variables, round to nearest 0.001
temp = round(double(temp0)*1000.0d0)/1000.0d0
logg = round(double(logg0)*1000.0d0)/1000.0d0
metal = round(double(metal0)*1000.0d0)/1000.0d0
alpha = round(double(alpha0)*1000.0d0)/1000.0d0

; Synthetic grid structure
; SYNSTR must have five tags: file, TEFF, LOGG, METAL, ALPHA
; TEFF, LOGG, METAL and ALPHA should be the arrays of UNIQUE stellar values
; FILE should be a NteffxNloggxNmetalxNalpha string array with the model filenames
;  If the model doesn't exist then file should be an empty string.
type = size(synstr,/type)
; Wrong type
if type ne 8 then begin
  error = 'SYNSTR must be the synthetic grid structure'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'FILE') eq 0 then begin
  error = 'SYNSTR does not have FILE tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'TEFF') eq 0 then begin
  error = 'SYNSTR does not have TEFF tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'LOGG') eq 0 then begin
  error = 'SYNSTR does not have LOGG tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'METAL') eq 0 then begin
  error = 'SYNSTR does not have METAL tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'ALPHA') eq 0 then begin
  error = 'SYNSTR does not have ALPHA tag'
  if not keyword_set(silent) then print,error
  return
endif


; Get the unique stellar parameter arrays
; They should already be unique
teffarr = synstr.teff
loggarr = synstr.logg
metalarr = synstr.metal
alphaarr = synstr.alpha



;---------------
; Normal
;---------------
if not keyword_set(interp) and not keyword_set(nearest) then begin

  mch = machar(/double)
  eps = mch.eps

  tongrid = where(abs(teffarr-temp) lt 10.*eps,ntongrid)
  longrid = where(abs(loggarr-logg) lt 10.*eps,nlongrid)
  mongrid = where(abs(metalarr-metal) lt 10.*eps,nmongrid)
  aongrid = where(abs(alphaarr-alpha) lt 10.*eps,naongrid)

  if ntongrid eq 0 or nlongrid eq 0 or nmongrid eq 0 or naongrid eq 0 then begin
    error = 'NO Spectrum for Teff='+strtrim(temp,2)+' logg='+strtrim(logg,2)+' [M/H]='+strtrim(metal,2)+' [alpha/Fe]='+strtrim(alpha,2)
    if not keyword_set(silent) then print,error
    return
  endif

  filename = synstr.file[tongrid[0],longrid[0],mongrid[0],aongrid[0]]

  ; No file
  if strtrim(filename,2) eq '' then begin
    error = 'NO Spectrum for Teff='+strtrim(temp,2)+' logg='+strtrim(logg,2)+' [M/H]='+strtrim(metal,2)+' [alpha/Fe]='+strtrim(alpha,2)
    if not keyword_set(silent) then print,error
    return
  endif

  ; Does the spectrum exist
  test = FILE_TEST(filename)
  if test eq 0 then begin
    error = filename+' NOT FOUND'
    if not keyword_set(silent) then print,error
    return
  endif

  ; Restore the spectrum
  undefine,wave,spec,head
  RDFITS_SPEC,filename,wave,spec,head=head,error=error
  if n_elements(error) gt 0 then begin
    if not keyword_set(silent) then print,error
    return
  endif

endif


;---------------
; Interpolating
;---------------
if keyword_set(interp) then begin
  SPECFIT_MODELINTERP,synstr,temp,logg,metal,alpha,wave,spec,error=error,head=head,silent=silent
  return
endif


;---------------
; Nearest Model
;---------------
if keyword_set(nearest) then begin
  SPECFIT_MODELINTERP,synstr,temp,logg,metal,alpha,wave,spec,error=error,head=head,/nearest,silent=silent
  return
endif


if keyword_set(stp) then stop

end
