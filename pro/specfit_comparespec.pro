pro specfit_comparespec,wave1,spec1,wave2,spec2,specpars,chisq,vrel,vrelerr,zerovel=zerovel,$
                        error=error,stp=stp,plot=pl,errspec=errspec0,mask=mask0,silent=silent

;+
;
; SPECFIT_COMPARESPEC
;
; This program compares two spectra and returns the best fitting
; chisq and relative velocity.  The two spectra are assumed to both
; have the same resolution and will NOT be smoothed.
;
; INPUTS:
;  wave1     The first wavelength array
;  spec1     The first spectrum array
;  wave2     The second wavelength array
;  spec2     The second spectrum array
;  specpars  The spectrum parameters, specpars = [w0, w1, dw, res]
;              NO smoothing is performed.
;  =errspec  The error spectrum for spec1.
;  =mask     The mask to use with wave1/spec1.  1-for points to use,
;              0-for points to ignore.
;              This can also be an array of weights.
;  /zerovel  Do not solve for velocity, should already be set to rest frame.
;  /plot     Plot the spectra.
;  /silent   Don't print anything to the screen
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  chisq   The chi squared of the best fit.
;  vrel    The relative velocity of the best fit
;  vrelerr The 1-sigma uncertainty of Vrel derived from the ChiSq
;            distribution.  This is only meaningful if "errspec" was input.
;  =error  The error message if there was one.
;
; USAGE:
;  IDL>specfit_comparespec,wave1,spec1,wave2,spec2,[4000.,4500.,1.4,2.7],chisq,vrel,vrelerr
;
; By D.Nidever  Oct 2008
;-

undefine,error,chisq,vrel,vrelerr

nwave1 = n_elements(wave1)
nspec1 = n_elements(spec1)
nwave2 = n_elements(wave2)
nspec2 = n_elements(spec2)
nspecpars = n_elements(specpars)

; Not enough inputs
if nwave1 eq 0 or nspec1 eq 0 or nwave2 eq 0 or nspec2 eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_comparespec,wave1,spec1,wave2,spec2,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,zerovel=zerovel,'
  print,'                             plot=plot,error=error,silent=silent,stp=stp'
  error = 'Not enough inputs'
  return
endif

cspeed = 2.99792458d5  ; speed of light in km/s

; Wave1 and spec1 not of the same length
if nwave1 ne nspec1 then begin
  error = 'WAVE1 and SPEC1 not of the same length'
  if not keyword_set(silent) then print,error
  return
endif

; Wave2 and spec2 not of the same length
if nwave2 ne nspec2 then begin
  error = 'WAVE2 and SPEC2 not of the same length'
  if not keyword_set(silent) then print,error
  return
endif

; Not enough spectrum parameters
if nspecpars lt 4 then begin
  error = 'Not enough spectrum parameters.  Need specpars=[w0,w1,dw,res]'
  if not keyword_set(silent) then print,error
  return
endif

; No error spectrum input
if n_elements(errspec0) eq 0 then errspec=0 else errspec=errspec0


; Error spectrum not of the right size
if n_elements(errspec0) gt 0 and n_elements(errspec0) ne nspec1 then begin
  error = 'ERRSPEC must have same number of elements as SPEC1/WAVE1'
  if not keyword_set(silent) then print,error
  return
endif


; Mask input
nmask = n_elements(mask0)
if nmask gt 0 then begin

  mask = mask0

  ; Not the right size
  if nmask ne nspec1 then begin
    error = 'MASK and SPEC1 not of the same length'
    if not keyword_set(silent) then print,error
    return
  endif

  ; No good points
  gdmask = where(mask gt 0.0,ngdmask)
  if ngdmask lt 1 then begin
  ;if total(mask) lt 1 then begin
    error = 'NO good points in mask'
    if not keyword_set(silent) then print,error
    return
  endif

  ;bd = where(mask ne 1 and mask ne 0,nbd)
  ;if nbd gt 0 then begin
  ;  error = 'MASK must be 0s and 1s ONLY'
  ;  if not keyword_set(silent) then print,error
  ;  return
  ;endif

endif else begin
  mask = wave1*0.+1.0
endelse

; Spectrum parameters
w0 = specpars[0]
w1 = specpars[1]
dw = specpars[2]
res = specpars[3:*]  ; might be an array


; Prepare the spectra if necessary
;----------------------------------
if abs(min(wave1)-w0) gt 1e-4 or abs(max(wave1)-w1) gt 1e-4 then begin
  SPECFIT_PREPSPEC,wave1,spec1,errspec,mask,specpars,ww1,ss1,ee1,mm1,/object
  nww1 = n_elements(ww1)

endif else begin
  ww1 = wave1
  ss1 = spec1
  ee1 = errspec
  mm1 = mask
endelse
nww1 = n_elements(ww1)

if abs(min(wave2)-w0) gt 1e-4 or abs(max(wave2)-w1) gt 1e-4 then begin
  SPECFIT_PREPSPEC,wave2,spec2,0,0,specpars,ww2,ss2,/object
  nww2 = n_elements(ww2)
endif else begin
  ww2 = wave2
  ss2 = spec2
endelse
nww2 = n_elements(ww2)

; The two spectra are not of the same length
if nww1 ne nww2 then begin
  error = 'WW1 and WW2 not of the same length'
  if not keyword_set(silent) then print,error
  return
endif


;--------------------------
; Compare the two spectra
;--------------------------


; If I divide by the synthetic spectrum then the 
; chisq is not going to be the same thing from
; syn spectrum to syn spectrum.

; No error spectrum input
if n_elements(errspec0) eq 0 then ee1=fltarr(nspec1)+1.0


; NO SHIFT
; Just do chisq
;-----------
if keyword_set(zerovel) then begin

  ; Pearson's chi square test
  ; chisq = total( (O-E)^2/E )

  ; Chi squared
  ; chisq = total( (O-E)^2/sig^2)
  dof = total(mm1)       ; degrees of freedom
  chisq = total( ((ss1-ss2)*mm1)^2.0/ee1^2.0 )
  ;chisq = total( (ss1-ss2)^2.0/ss2 )/nww1

  ; No shift
  vrel = 0.0
  vrelerr = 0.0

; SHIFTING
;-----------
endif else begin

  dwlog = median(slope(alog10(ww1)))
  dv = ( 10.d0^(dwlog)-1.0d0 )*cspeed
  maxshift = ceil(1000.0/dv)
  
  XCORLB,ss1,ss2,maxshift,shift,chisq,shifterr,errspec=ee1,mask=mm1

  ;XCORLB,ss1,ss2,10,shift,chisq,shifterr,errspec=ee1,mask=mm1
  ;
  ;; If the range was not big enough then widen it
  ;if abs(shift) gt 9.0 then XCORLB,ss1,ss2,25,shift,chisq,shifterr,errspec=ee1,mask=mm1
  ;if abs(shift) gt 24.0 then XCORLB,ss1,ss2,50,shift,chisq,shifterr,errspec=ee1,mask=mm1
  ;if abs(shift) gt 45.0 then XCORLB,ss1,ss2,100,shift,chisq,shifterr,errspec=ee1,mask=mm1
  ;if abs(shift) gt 90.0 then XCORLB,ss1,ss2,200,shift,chisq,shifterr,errspec=ee1,mask=mm1

  ; Convert shift to velocity
  ; delta log(wave) = log(v/c+1)
  ; v = (10^(delta log(wave))-1)*c
  dwlog = median(slope(alog10(ww1)))
  vrel = ( 10.d0^(shift*dwlog)-1.0d0 )*cspeed

  ; Vrel uncertainty
  dvreldshift = alog(10.0)*(10.d0^(shift*dwlog))*dwlog*cspeed  ; derivative wrt shift
  vrelerr = dvreldshift * shifterr

  dof = total(mm1)-1  ; degrees of freedom

endelse

rchisq = chisq/dof  ; reduced ChiSq


; Plot the spectra
if keyword_set(pl) then begin
  ;plot,[0],[0],xtit='Wavelength (A)',xr=[w0,w1],yr=[0.0,max([spec2,ssyn])>1.2],xs=1,ys=1
  plot,[0],[0],xtit='Wavelength (A)',xr=[w0,w1],yr=[0.0,max([ss1,ss2])>1.2],xs=1,ys=1
  oplot,ww1,ss1
  if nmask gt 0 then begin
    binarymask2 = float(mm1 gt 0.0)
    breaks1 = where(slope([0,binarymask2]) eq 1,nbreaks1)
    breaks2 = where(slope([binarymask2,0]) eq -1,nbreaks2)
    ;breaks1 = where(slope([0,mm1]) eq 1,nbreaks1)
    ;breaks2 = where(slope([mm1,0]) eq -1,nbreaks2)
    for i=0,nbreaks1-1 do begin
      oplot,ww1[breaks1[i]:breaks2[i]],ss1[breaks1[i]:breaks2[i]],co=150
      oplot,ww1[breaks1[i]:breaks2[i]],(ss2-s1)[breaks1[i]:breaks2[i]],co=150
    end
    ;gd = where(mm1 gt 0.0,ngd)
    ;oplot,ww1[gd],ss1[gd],ps=1,co=150
  endif
  if n_elements(errspec0) gt 1 then oplot,ww1,ee1
  oplot,ww2*(vrel/cspeed+1.0),ss2,co=250,linestyle=5
  xyouts,mean(ww1),0.2,'  R.ChiSq='+strtrim(string(rchisq,format='(F7.4)'),2)+' Vrel='+strtrim(string(vrel,format='(F9.3)'),2)+$
         ' +/-'+strtrim(string(vrelerr,format='(F9.3)'),2)+' km/s',align=0.5,charsize=1.2
  ;xyouts,mean(ww1),0.2,'Teff='+strtrim(string(teff,format='(I5)'),2)+' Logg='+string(logg,format='(F4.2)')+' '+$
  ;       '[M/H]='+strtrim(string(metal,format='(F5.2)'),2)+' [alpha/Fe]='+string(alpha,format='(F5.2)')+$
  ;       '  R.ChiSq='+strtrim(string(rchisq,format='(F7.4)'),2)+' Vrel='+strtrim(string(vrel,format='(F9.3)'),2)+$
  ;       ' +/-'+strtrim(string(vrelerr,format='(F9.3)'),2)+' km/s',align=0.5,charsize=1.2
endif


if keyword_set(stp) then stop

end
