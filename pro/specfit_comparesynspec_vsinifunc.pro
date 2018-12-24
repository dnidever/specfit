function specfit_comparesynspec_vsinifunc,x,par,wsyn0=wsyn0,ssyn0=ssyn0,specpars,$
          prepared=prepared,wave2=wave2,spec2=spec2,errspec2=errspec2,mask2=mask2

; This compares a rotationally broadened synthetic spectrum
; with an observational spectrum using MPFITFUN


; I DON'T THINK THIS IS GOING TO WORK BECAUSE YOU NEED TO SPLINE EVERY
; TIME AND IT WILL BE SLOW!!!!
stop

ivsini = par[0]

; Prepare the synthetic spectrum
;-------------------------------
SPECFIT_PREPSPEC,wsyn0,ssyn0,0,0,specpars,wsyn,ssyn,prepared=prepared,vsini=ivsini
nwsyn = n_elements(wsyn)

; The two spectra are not of the same length
if nwsyn ne nwave then begin
  error = 'WAVE and WSYN not of the same length'
  if not keyword_set(silent) then print,error
  return
endif


;--------------------------
; Compare the two spectra
;--------------------------
SPECFIT_COMPARESPEC,wave2,spec2,wsyn,ssyn,specpars,chisq,vrel,vrelerr,errspec=errspec2,mask=mask2,$
                    zerovel=zerovel,error=error,silent=silent

if n_elements(error) ne 0 then begin
  if not keyword_set(silent) then print,'ERROR - ',error
  return
endif

wsynth = wsyn0
ssynth = CSPLINE(wsyn0*(vrel/cspeed+1.0d0),ssyn0,wsyn0)

return,(spec2-ssynth)/errspec2

stop

end
