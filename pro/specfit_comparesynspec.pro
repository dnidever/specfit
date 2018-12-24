pro specfit_comparesynspec,wave0,spec0,synstr,synpars,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec0,$
                           mask=mask0,zerovel=zerovel,error=error,silent=silent,stp=stp,plot=pl,$
                           wsynth=wsynth,ssynth=ssynth,synhead=synhead,interp=interp,nearest=nearest,$
                           fitvsini=fitvsini,fixvsini=fixvsini

;+
;
; SPECFIT_COMPARESYNSPEC
;
; This program compares a synthetic spectrum to an observed spectrum.
;
; INPUTS:
;  wave      The observed wavelength array
;  spec      The observed spectrum array
;  synstr    The synthetic grid structure.  It should have TEFF, LOGG, METAL, ALPHA and FILE
;              tags.  TEFF, LOGG, METAL and ALPHA should be 1D arrays of the UNIQUE stellar
;              parameters of the synthetic grid.  FILE should be a 4D string array (or 3D if
;              ALPHA only has one element) of dimensions [Nteff,Nlogg,Nmetal,Nalpha] with
;              the absolute path of the synthetic spectrum FITS file.  It should be an empty
;              string if the spectrum does not exist.
;  synpars   The stellar parameters for the synthetic spectrum
;              synpars = [Teff, logg, metal, alpha]
;  specpars  The spectrum parameters, specpars = [w0, w1, dw, res]. DW
;              should be in log(Angstroms)
;  =errspec  The error spectrum for the observed spectrum.
;  =mask     The mask to use with wave/spec.  1-for points to use,
;              0-for points to ignore.
;              This can also be an array of weights.
;  /interp   Interpolate synthetic spectrum between grid points.
;  /nearest  Use the closest synthetic spectrum to the input stellar
;              parameters instead of interpolating.
;  /zerovel  Do not solve for velocity, should already be set to rest frame.
;  /fitvsini  Fit vsini.
;  =fixvsini  Fix vsini at a certain value (in km/s).
;  /plot     Plot the spectra.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  chisq     The chi squared of the best fit.
;  vrel      The relative velocity of the best fit
;  vrelerr   The 1-sigma uncertainty of Vrel derived from the ChiSq
;              distribution.  This is only meaningful if "errspec" was input.
;  vsini     The best fitting vsini.
;  =wsynth   The final synthetic spectrum wavelength array.
;  =ssynth   The final synthetic spectrum (shifted to the best-fitting
;              velocity)
;  =synhead  The header of the synthetic spectrum.
;  =error    The error message if there was one.
;
; USAGE:
;  IDL>specfit_comparesynspec,synstr,wave,spec,synstr,[5000.,2.5,-1.5,0.0],[4000.,4500.,1.4,2.7],chisq,vrel,vrelerr
;
; By D.Nidever  Oct 2008
;-

undefine,error,chisq,vrel,vrelerr,wsynth,ssynth

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nsynstr = n_elements(synstr)
nsynpars = n_elements(synpars)
nspecpars = n_elements(specpars)
ndir = n_elements(dir)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nsynstr eq 0 or nsynpars eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_comparesynspec,wave,spec,synstr,synpars,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec,'
  print,'                                mask=mask,nearest=nearest,zerovel=zerovel,fitvsini=fitvsini,'
  print,'                                plot=plot,error=error,wsynth=wsynth,ssynth=ssynth,synhead=synhead,'
  print,'                                fixvsini=fixvsini,silent=silent,stp=stp'
  error = 'Not enough inputs'
  return
endif

cspeed = 2.99792458d5  ; speed of light in km/s

; Internal arrays
wave = wave0
spec = spec0

; Wave and spec not of the same length
if nwave ne nspec then begin
  error = 'WAVE and SPEC not of the same length'
  if not keyword_set(silent) then print,error
  return
endif

; Mask input
nmask = n_elements(mask0)
if nmask gt 0 then begin

  mask = mask0

  ; Not the right size
  if nmask ne nspec then begin
    error = 'MASK and SPEC not of the same length'
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
  mask = wave*0.+1.0
endelse


; Checking the synthetic spectrum grid structure
type = size(synstr,/type)
if type ne 8 then begin
  error = 'SYNSTR must be a structure'
  if not keyword_set(silent) then print,error
  return
endif

tags = tag_names(synstr)
indteff = max(strpos(tags,'TEFF'))
indlogg = max(strpos(tags,'LOGG'))
indmetal = max(strpos(tags,'METAL'))
indalpha = max(strpos(tags,'ALPHA'))
if tag_exist(synstr,'TEFF') eq 0 or tag_exist(synstr,'LOGG') eq 0 or tag_exist(synstr,'METAL') eq 0 or $
  tag_exist(synstr,'ALPHA') eq 0 or tag_exist(synstr,'FILE') eq 0 then begin
  error = 'SYNSTR must have TEFF, LOGG, METAL, ALPHA and FILE tags'
  if not keyword_set(silent) then print,error
  return
endif


; Not enough stellar parameters
if nsynpars lt 4 then begin
  error = 'Not enough stellar parameters.  Need synpars=[Teff,logg,metal,alpha]'
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
if n_elements(errspec0) gt 0 and n_elements(errspec0) ne nspec then begin
  error = 'ERRSPEC must have same number of elements as SPEC/WAVE'
  if not keyword_set(silent) then print,error
  return
endif

; Stellar parameters
teff = synpars[0]
logg = synpars[1]
metal = synpars[2]
alpha = synpars[3]

; Spectrum parameters
w0 = specpars[0]
w1 = specpars[1]
dw = specpars[2]
res = specpars[3:*]  ; might be an array


; Prepare the observed spectrum if necessary
;-------------------------------------------
if abs(min(wave0)-w0) gt 1e-4 or abs(max(wave0)-w1) gt 1d-4 or $
   median(spec0) lt 0.2 or median(spec0) gt 3.0 then begin

  SPECFIT_PREPSPEC,wave,spec,errspec,mask,specpars,wave2,spec2,errspec2,mask2,/object

; Already prepared
endif else begin

  wave2 = wave0
  spec2 = spec0
  if n_elements(errspec0) gt 0 then errspec2=errspec0
  mask2 = mask

endelse
nwave = n_elements(wave2)

; No error spectrum input
if n_elements(errspec0) eq 0 then undefine,errspec2



; Get the synthetic spectrum
;---------------------------
SPECFIT_GETSYNSPEC,synstr,teff,logg,metal,alpha,wsyn0,ssyn0,interp=interp,nearest=nearest,head=synhead,$
                   error=error,silent=silent

if n_elements(error) gt 0 then begin
  if not keyword_set(silent) then print,error
  return
end


; Is this synthetic spectrum prepared already
prepared = SXPAR(synhead,'PREPARED',count=nprepared)

;; Already prepared CANNOT fit vsini
;if keyword_set(fitvsini) and prepared ne 0 then begin
;  error = 'CANNOT fit vsini because synthetic spectra are already PREPARED'
;  print,error
;  retall        ; FATAL error
;endif
;; Already prepared CANNOT fix vsini
;if keyword_set(fixvsini) and prepared ne 0 then begin
;  error = 'CANNOT fix vsini because synthetic spectra are already PREPARED'
;  print,error
;  retall        ; FATAL error
;endif


; Tried using MPFIT, NOT GOING TO BE EFFICIENT!!!
;parinfo = {limited:[1,1],limits:[0.0,1000.0]}
;fa = {wsyn0:wsyn0,ssyn0:ssyn0,specpars:specpars,prepared:prepared,wave2:wave2,spec2:spec2,$
;      errspec2:errspec2,mask2:mask2}
;func = 'specfit_comparesynspec_vsinifunc'
;fpar = MPFIT(func,[50.0],parinfo=parinfo,functargs=fa,status=status,$
;                bestnorm=bestnorm,maxiter=10,niter=niter,perror=perror)  ;/quiet


; Solving for vsini
;-------------------
;if keyword_set(fitvsini) and not keyword_set(fixvsini) then nvsini = 10 else nvsini=1
;vsiniarr = findgen(nvsini)*100
;vsiniarr = findgen(nvsini)*50    ; 0-450 km/s, in steps of 50 km/s
;vsiniarr = [0.0, 10.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350., 400.]
vsiniarr = [0.0, 10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350., 400.]
if not keyword_set(fitvsini) then vsiniarr=0.0
nvsini = n_elements(vsiniarr)

chisqarr = fltarr(nvsini)
vrelarr = fltarr(nvsini)
vrelerrarr = fltarr(nvsini)

For i=0,nvsini-1 do begin

  ; Broaden with vsini
  ivsini=0.0
  if keyword_set(fitvsini) then begin
    ;; vsini will be the doppler broadening, Gaussian sigma
    ;; convert to A and FWHM
    ;vsini_sigma_ang = vsiniarr[i]*median(wsyn0)/cspeed
    ;vsini_fwhm_ang = vsini_sigma_ang *2.35
    ;specpars1 = specpars
    ;specpars1[3] = sqrt(specpars[3]^2.0 + vsini_fwhm_ang^2.0)  ; add in quadrature
    ivsini = vsiniarr[i]
  endif
  if keyword_set(fixvsini) then ivsini=fixvsini  ; fixing vsini at input value

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

  chisqarr[i] = chisq
  vrelarr[i] = vrel
  vrelerrarr[i] = vrelerr

  ;stop

Endfor ; vsini loop


; Find best vsini fit
vsini = 0.0    ; default
if keyword_set(fitvsini) then begin

  vsiniarr2 = scale_vector(dindgen(100),0.0,max(vsiniarr))
  chisqarr2 = cspline(vsiniarr,chisqarr,vsiniarr2)
  bestind = first_el(minloc(chisqarr2))
  vsini = vsiniarr2[bestind[0]]

  ; Redo the "best" fit
  ;----------------------

  ;; Broaden with vsini
  ;vsini_sigma_ang = vsini*median(wsyn0)/cspeed
  ;vsini_fwhm_ang = vsini_sigma_ang *2.35
  ;specpars1 = specpars
  ;specpars1[3] = sqrt(specpars[3]^2.0 + vsini_fwhm_ang^2.0)  ; add in quadrature

  ; Prepare the synthetic spectrum
  ;-------------------------------
  SPECFIT_PREPSPEC,wsyn0,ssyn0,0,0,specpars,wsyn,ssyn,prepared=prepared,vsini=vsini
  nwsyn = n_elements(wsyn)

  ;--------------------------
  ; Compare the two spectra
  ;--------------------------
  SPECFIT_COMPARESPEC,wave2,spec2,wsyn,ssyn,specpars,chisq,vrel,vrelerr,errspec=errspec2,mask=mask2,$
                      zerovel=zerovel,error=error,silent=silent

  ;stop

end  ; vsini best fit

; vsini input
if keyword_set(fixvsini) then vsini=fixvsini    ; set to input value

; Get final synthetic spectrum
;  redshift the spectrum, but keep same wavelength array
wsynth = wsyn
;ssynth = SPLINE(wsyn*(vrel/cspeed+1.0d0),ssyn,wsyn)
ssynth = CSPLINE(wsyn*(vrel/cspeed+1.0d0),ssyn,wsyn)
;ssynth = CSPLINE(wsyn/(vrel/cspeed+1.0d0),ssyn,wsyn)
;ssynth = CSPLINE(wsyn,ssyn,wsyn/(vrel/cspeed+1.0))
;ssynth = CSPLINE(wsyn,ssyn,wsyn*(vrel/cspeed+1.0))


; Reduced ChiSq
dof = total(mask2)
if not keyword_set(zerovel) then dof=dof-1    ; fit one parameter
rchisq = chisq/dof

;pl=1

; Plot the spectra
if keyword_set(pl) then begin
  ;plot,[0],[0],/nodata,xtit='Wavelength (A)',xr=[w0,w1],yr=[0.0,max([spec2,ssyn])>1.2],xs=1,ys=1
  plot,[0],[0],/nodata,xtit='Wavelength (A)',xr=[w0,w1],yr=[-0.15,1.4],xs=1,ys=1
  oplot,wave2,spec2
  if nmask gt 0 then begin
    binarymask2 = float(mask2 gt 0.0)
    breaks1 = where(slope([0,binarymask2]) eq 1,nbreaks1)
    breaks2 = where(slope([binarymask2,0]) eq -1,nbreaks2)
    ;breaks1 = where(slope([0,mask2]) eq 1,nbreaks1)
    ;breaks2 = where(slope([mask2,0]) eq -1,nbreaks2)
    for i=0,nbreaks1-1 do begin
      oplot,wave2[breaks1[i]:breaks2[i]],spec2[breaks1[i]:breaks2[i]],co=150
      oplot,wave2[breaks1[i]:breaks2[i]],(spec2-ssynth)[breaks1[i]:breaks2[i]],co=150
    end
    ;gd = where(mask2 gt 0.0,ngd)
    ;oplot,wave2[gd],spec2[gd],ps=1,co=150
  endif
  ;oplot,wsyn*(vrel/cspeed+1.0),ssyn,co=250,linestyle=5
  oplot,wsynth,ssynth,co=250,linestyle=5
  oplot,[0,1e4],[0,0],linestyle=2
  xyouts,mean(wave2),0.2,'Teff='+strtrim(string(teff,format='(I5)'),2)+' Logg='+string(logg,format='(F4.2)')+' '+$
         '[M/H]='+strtrim(string(metal,format='(F5.2)'),2)+' [alpha/Fe]='+string(alpha,format='(F5.2)'),align=0.5,charsize=1.2
  out2 = 'R.ChiSq='+strtrim(string(rchisq,format='(F7.4)'),2)+' Vrel='+strtrim(string(vrel,format='(F9.3)'),2)+$
         ' +/-'+strtrim(string(vrelerr,format='(F9.3)'),2)+' km/s'
  if keyword_set(fitvsini) or keyword_set(fixvsini) then out2=out2+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+' km/s'
  xyouts,mean(wave2),0.15,out2,align=0.5,charsize=1.2
endif

;wait,1
if keyword_set(stp) then stop

end
