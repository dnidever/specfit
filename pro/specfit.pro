pro specfit,wave0,spec0,errspec0,synstr,specpars,fitstr,mask=mask0,$
                 zerovel=zerovel,monte=monte,nmonte=nmonte0,error=error,$
                 stp=stp,plot=pl,silent=silent,fitstr_monte=monte_fitstr,$
                 fitvsini=fitvsini,fixvsini=fixvsini,psfile=psfile

;+
;
; SPECFIT
;
; This program fits synthetic spectra to an observed spectrum
;
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
;  =psfile    EPS filename to save the final plot to.
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
;
; By D.Nidever  Nov 2008
;-


; FUTURE SPECFIT IMPROVEMENTS:
; Add specpars to the "synstr" structures so it's clear what
; the parmeters are and so it can easily be used in
; specfit_driver.pro and similar programs.  Also, whether
; it's log or not.  Probably should change specfit so
; that specpars has the same format as in prepare_syngrid
;
; add another option to specfit to do an iterative search for
; velocity.  Use a rough grid of 3x3x3 (Teff, logg, [M/H]) to
; get a first guess of velocity.  Then use that velocity to set
; the spectrum to the rest frame and use the hybrid model
; with /zerovel to find best parameters.  Then get second
; guess of velocity using the best-fit spectrum.  Iterate
; one more time and get third velocit at the end.  Call
; this method /itervel.  Make the current default option
; /fitvel, but make /itervel the new default.
;
; Of course also add the "postsmooth" option.
;
; Also allow there be a "quickgrid", a small grid that
; gets compared all at once to a spectrum.


t0 = systime(1)

undefine,error,chisq,vrel,bestpars

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)
nsynstr = n_elements(synstr)
nspecpars = n_elements(specpars)
ndir = n_elements(dir)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 or nsynstr eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit,wave,spec,errspec,synstr,specpars,fitstr,mask=mask,'
  print,'                      zerovel=zerovel,fitvsini=fitvsini,monte=monte,'
  print,'                      nmonte=nmonte,fixvsini=fixvsini,plot=plot,error=error,'
  print,'                      psfile=psfile,stp=stp'
  error = 'Not enough inputs'
  return
endif

cspeed = 2.99792458d5  ; speed of light in km/s

; Wave and spec not of the same length
if nwave ne nspec then begin
  error = 'WAVE and SPEC not of the same length'
  print,error
  return
endif

; Spec and errspec not of the same length
if nspec ne nerrspec then begin
  error = 'SPEC and ERRSPEC not of the same length'
  print,error
  return
endif

; Checking that we have some finite elements
dum = where(finite(wave0) eq 1,ngdwave)
if ngdwave eq 0 then begin
  error = 'WAVE had NOT finite elements'
  print,error
  return
endif
dum = where(finite(spec0) eq 1,ngdspec)
if ngdspec eq 0 then begin
  error = 'SPEC had NOT finite elements'
  print,error
  return
endif
if nerrspec gt 0 then begin
  dum = where(finite(errspec0) eq 1,ngderrspec)
  if ngderrspec eq 0 then begin
    error = 'ERRSPEC had NOT finite elements'
    print,error
    return
  endif
endif
if n_elements(mask0) gt 0 then begin
  dum = where(finite(mask0) eq 1,ngdmask)
  if ngdmask eq 0 then begin
    error = 'MASK had NOT finite elements'
    print,error
    return
  endif
endif


; Internal arrays
wave = wave0
spec = spec0

; Mask input
nmask = n_elements(mask0)
if nmask gt 0 then begin

  mask = mask0

  ; Not the right size
  if nmask ne nspec then begin
    error = 'MASK and SPEC not of the same length'
    print,error
    return
  endif

  ; No good points
  gdmask = where(mask gt 0.0,ngdmask)
  if ngdmask lt 1 then begin
  ;if total(mask) lt 1 then begin
    error = 'NO good points in mask'
    print,error
    return
  endif

  ;bd = where(mask ne 1 and mask ne 0,nbd)
  ;if nbd gt 0 then begin
  ;  error = 'MASK must be 0s and 1s ONLY'
  ;  print,error
  ;  return
  ;endif

endif else begin
  mask = wave0*0.0+1.0
endelse


; Checking the synthetic spectrum grid structure
type = size(synstr,/type)
if type ne 8 then begin
  error = 'SYNSTR must be a structure'
  print,error
  return            
endif

if tag_exist(synstr,'TEFF') eq 0 or tag_exist(synstr,'LOGG') eq 0 or tag_exist(synstr,'METAL') eq 0 or $
  tag_exist(synstr,'ALPHA') eq 0 or tag_exist(synstr,'FILE') eq 0 then begin
  error = 'SYNSTR must have TEFF, LOGG, METAL, ALPHA and FILE tags'
  print,error
  return
endif


; Not enough spectrum parameters
if nspecpars lt 4 then begin
  error = 'Not enough spectrum parameters.  Need specpars=[w0,w1,dw,res]'
  print,error
  return
endif

; Compile PUSH.PRO
push,0

; Spectrum parameters
w0 = specpars[0]
w1 = specpars[1]
dw = specpars[2]
res = specpars[3:*]  ; might be an array


wavein = wave
specin = spec
errin = errspec0
maskin = mask

;#########################################
; Use SPECFIT_HYBRID to get the solution
;#########################################
SPECFIT_HYBRID,wavein,specin,errin,synstr,specpars,hyb_fitstr,zerovel=zerovel,mask=maskin,$
               fitvsini=fitvsini,fixvsini=fixvsini,silent=silent,plot=plot

bestpars = hyb_fitstr.bestpars
vrel = hyb_fitstr.vrel
vsini = hyb_fitstr.vsini
chisq = hyb_fitstr.chisq
dof = hyb_fitstr.dof
rchisq = hyb_fitstr.rchisq
npts = hyb_fitstr.npts
nused = hyb_fitstr.nused
; SNR
snr = median(spec/errspec0)


;#####################
; Monte Carlo errors
;#####################
if keyword_set(monte) then begin

  print,''
  print,'Using Monte Carlo simulation to obtain parameter uncertainties'

  ; Number of Monte Carlo simulations
  if keyword_set(nmonte0) then nmonte=nmonte0 else nmonte=50
  if nmonte lt 10 then nmonte=10   ; minimum of 10 simulations

  ; Get best model spectrum
  SPECFIT_GETSYNSPEC,synstr,bestpars[0],bestpars[1],bestpars[2],bestpars[3],wsyn0,ssyn0,head=synhead

  ; Is this synthetic spectrum prepared already
  prepared = SXPAR(synhead,'PREPARED',count=nprepared)

  ; Prepare the synthetic spectrum - if necessary
  SPECFIT_PREPSPEC,wsyn0,ssyn0,0,0,specpars,wsyn,ssyn2,prepared=prepared,vsini=vsini

  ; Shift to best Vrel
  wsyn2 = wsyn*(vrel/cspeed+1.0d0)
  ; Spline onto wave
  ssyn3 = CSPLINE(wsyn2,ssyn2,wave)
  ; Get object continuum
  SPECFIT_NORM,wavein,specin,errin,dum,cont
  ssyn4 = ssyn3*cont

  wmonte = wave
  nwave = n_elements(wave)
  errmonte = errspec0
  maskmonte = mask

  undefine,monte_fitstr
  ; Monte Carlo loop
  out =''
  for j=0,nmonte-1 do begin

    ; Progress bar
    len = strlen(out)  ; length of previous string
    out = strtrim(j+1,2)+'/'+strtrim(nmonte,2)
    if j eq 0 then writeu,-1,out else writeu,-1,strjoin(replicate(string(8B),len))+out

    ; Add noise to the spectrum
    smonte = ssyn4 + RANDOMN(seed,nwave)*errmonte
    
    ; Fit it
    undefine,fitstr2
    SPECFIT_HYBRID,wmonte,smonte,errmonte,synstr,specpars,fitstr2,zerovel=zerovel,mask=maskmonte,$
                   fitvsini=fitvsini,fixvsini=fixvsini,/silent

    ; Add to structure
    PUSH,monte_fitstr,fitstr2

  end
  print,''

  ; Calculate parameter uncertainties
  ;  Need to use RMS instead of STDDEV in case the mean is not
  ;  around the right value
  teffpars = monte_fitstr.bestpars[0]
  tefferr = sqrt(mean( (teffpars-bestpars[0])^2 ))
  ;tefferr = STDDEV(teffpars)
  loggpars = monte_fitstr.bestpars[1]
  loggerr = sqrt(mean( (loggpars-bestpars[1])^2 ))
  ;loggerr = STDDEV(loggpars)
  metalpars = monte_fitstr.bestpars[2]
  metalerr = sqrt(mean( (metalpars-bestpars[2])^2 ))
  ;metalerr = STDDEV(metalpars)
  alphapars = monte_fitstr.bestpars[3]
  alphaerr = sqrt(mean( (alphapars-bestpars[3])^2 ))
  ;alphaerr = STDDEV(alphapars)
  vrelpars = monte_fitstr.vrel
  vrelerr = sqrt(mean( (vrelpars-vrel)^2 ))
  ;vrelerr = STDDEV(vrelpars)
  vsinipars = monte_Fitstr.vsini
  vsinierr = sqrt(mean( (vsinipars-vsini)^2 ))
  ;vsinierr = STDDEV(vsinipars)

  ; Print Monte Carlo errors
  if not keyword_set(silent) then begin
    print,'Monte Carlo Teff error       = ',strtrim(string(tefferr,format='(I5)'),2)+' K'
    print,'Monte Carlo logg error       = ',strtrim(string(loggerr,format='(F4.2)'),2)
    print,'Monte Carlo [M/H] error      = ',strtrim(string(metalerr,format='(F5.2)'),2)
    print,'Monte Carlo [alpha/Fe] error = ',strtrim(string(alphaerr,format='(F5.2)'),2)
    print,'Monte Carlo Vrel error       = ',strtrim(string(vrelerr,format='(F5.2)'),2)
    if keyword_set(fitvsini) then $
      print,'Monte Carlo Vsini error       = ',strtrim(string(vsinierr,format='(F5.2)'),2)
  endif

  ; Error must be at least 0.5*step
  if n_elements(synstr.teff) gt 1 then dteff = synstr.teff[1]-synstr.teff[0] else dteff=0.0
  if n_elements(synstr.logg) gt 1 then dlogg = synstr.logg[1]-synstr.logg[0] else dlogg=0.0
  if n_elements(synstr.metal) gt 1 then dmetal = synstr.metal[1]-synstr.metal[0] else dmetal=0.0
  if n_elements(synstr.alpha) gt 1 then dalpha = synstr.alpha[1]-synstr.alpha[0] else dalpha=0.0
  uncteff = tefferr > 0.5*dteff
  unclogg = loggerr > 0.5*dlogg
  uncmetal = metalerr > 0.5*dmetal
  uncalpha = alphaerr > 0.5*dalpha
  uncvrel = vrelerr
  parerrors = [uncteff, unclogg, uncmetal, uncalpha]

  ; Final uncertainties
  if dteff eq 0.0 then uncteff = 0.0
  if dlogg eq 0.0 then unclogg = 0.0
  if dmetal eq 0.0 then uncmetal = 0.0
  if dalpha eq 0.0 then uncalpha = 0.0


  ; Print final output
  if not keyword_set(silent) then begin
    print,''
    print,'Final Fitted parameters with Monte Carlo errors:'
    print,'Teff = ',string(bestpars[0],format='(I5)')+' +/- '+strtrim(string(uncteff,format='(I5)'),2)+' K'
    print,'Logg = ',string(bestpars[1],format='(F4.2)')+' +/- '+strtrim(string(unclogg,format='(F6.2)'),2)
    print,'[M/H] = ',strtrim(string(bestpars[2],format='(F5.2)'),2)+' +/- '+strtrim(string(uncmetal,format='(F5.2)'),2)
    print,'[alpha/Fe] = ',strtrim(string(bestpars[3],format='(F5.2)'),2)+' +/- '+strtrim(string(uncalpha,format='(F5.2)'),2)
    print,'Vrel = ',strtrim(string(vrel,format='(F9.1)'),2)+' +/- '+strtrim(string(uncvrel,format='(F9.1)'),2)+' km/s'
    if keyword_set(fitvsini) then print,'Best vsini = ',strtrim(string(vsini,format='(F7.1)'),2)
    if keyword_set(fixvsini) then print,'Fixed vsini = ',strtrim(string(fixvsini,format='(F7.1)'),2)
    print,'R.ChiSq = ',string(rchisq,format='(F7.3)')
    print,'S/N = ',strtrim(string(snr,format='(F6.1)'),2)
  endif

  ; Make final FITSTR
  if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
  if keyword_set(zerovel) then fixzero=1 else fixzero=0
  fitstr = {specpars:double(specpars),zerovel:fix(fixzero),mask:fix(maskinput),npts:long(npts),nused:long(nused),snr:double(snr),bestpars:double(bestpars),$
            parerrors:double(parerrors),teff:double(bestpars[0]),tefferr:double(uncteff),logg:double(bestpars[1]),loggerr:double(unclogg),metal:double(bestpars[2]),$
            metalerr:double(uncmetal),alpha:double(bestpars[3]),alphaerr:double(uncalpha),chisq:float(chisq),dof:long(dof),$
            rchisq:float(rchisq),vrel:double(vrel),vrelerr:double(uncvrel),vsini:double(vsini),vsinierr:double(0.0)}

; NO monte carlo errors
endif else begin

  ; Print final output
  if not keyword_set(silent) then begin
    print,''
    print,'Final Fitted parameters:'
    print,'Teff = ',string(bestpars[0],format='(I5)')
    print,'Logg = ',string(bestpars[1],format='(F4.2)')
    print,'[M/H] = ',strtrim(string(bestpars[2],format='(F5.2)'),2)
    print,'[alpha/Fe] = ',strtrim(string(bestpars[3],format='(F5.2)'),2)
    print,'Vrel = ',strtrim(string(vrel,format='(F9.1)'),2)
    if keyword_set(fitvsini) then print,'Vsini = ',strtrim(string(vsini,format='(F9.1)'),2)
    if keyword_set(fixvsini) then print,'Fixed Vsini = ',strtrim(string(fixvsini,format='(F9.1)'),2)
    print,'R.ChiSq = ',string(rchisq,format='(F7.3)')
    print,'S/N = ',strtrim(string(snr,format='(F6.1)'),2)
  endif

  ; Make final FITSTR
  if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
  if keyword_set(zerovel) then fixzero=1 else fixzero=0
  fitstr = {specpars:double(specpars),zerovel:fix(fixzero),mask:fix(maskinput),npts:long(npts),nused:long(nused),snr:double(snr),bestpars:double(bestpars),$
            teff:double(bestpars[0]),logg:double(bestpars[1]),metal:double(bestpars[2]),alpha:double(bestpars[3]),chisq:float(chisq),dof:long(dof),$
            rchisq:float(rchisq),vrel:double(vrel),vrelerr:double(hyb_fitstr.vrelerr),vsini:double(vsini),vsinierr:double(0.0)}

endelse


; Plot
if keyword_set(pl) or keyword_set(psfile) then begin

  ; Get prepared object spectrum
  SPECFIT_PREPSPEC,wave,spec,errspec0,mask,specpars,wave2,spec2,errspec2,mask2,/object
  ; Get the best fitting, velocity-shifted synthetic spectrum
  SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,bestpars,specpars,errspec=errspec2,mask=mask2,$
                         zerovel=zerovel,fitvsini=fitvsini,fixvsini=fixvsini,error=error,$
                         /silent,wsynth=wsynth,ssynth=ssynth

  if keyword_set(psfile) then begin
    setdisp
    loadct,39,/silent
    !p.font = 0
    ps_open,psfile,/color,thick=4,/encap
  endif

  plot,[0],[0],/nodata,xtit='Wavelength (A)',xr=[w0,w1],yr=[-0.15,1.4],xs=1,ys=1
  oplot,wave2,spec2
  oplot,wave2,spec2-ssynth  ; difference
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

    ;for i=0,n_elements(wave2)-1 do begin
    ;  if i gt 1 then if mask2[i] eq 1.0 and mask2[i-1] eq 1.0 then $
    ;    oplot,wave2[i-1:i],spec2[i-1:i],co=150
    ;end
    ;gd = where(mask2 gt 0.0,ngd)
    ;oplot,wave2[gd],spec2[gd],ps=1,co=150
  endif

  oplot,wsynth,ssynth,co=250,linestyle=5
  oplot,[w0,w1],[1,1],co=200,linestyle=2
  oplot,[w0,w1],[0,0],co=200,linestyle=2

  if keyword_set(monte) then begin
    out = 'Teff='+strtrim(string(bestpars[0],format='(I5)'),2)+'+/-'+strtrim(string(uncteff,format='(I5)'),2)+' K'+$
          ' Logg='+strtrim(string(bestpars[1],format='(F4.2)'),2)+'+/-'+strtrim(string(unclogg,format='(F6.2)'),2)+$
          ' [M/H]='+strtrim(string(bestpars[2],format='(F5.2)'),2)+'+/-'+strtrim(string(uncmetal,format='(F5.2)'),2)+$
          ' [alpha/Fe]='+strtrim(string(bestpars[3],format='(F5.2)'),2)+'+/-'+strtrim(string(uncalpha,format='(F5.2)'),2)
    xyouts,mean([w0,w1]),0.2,out,align=0.5,charsize=1.2
    out = ' Vrel='+strtrim(string(vrel,format='(F9.1)'),2)+'+/-'+strtrim(string(uncvrel,format='(F9.1)'),2)+' km/s'+$
          ' R.ChiSq='+strtrim(string(rchisq,format='(F7.3)'),2)+' S/N='+strtrim(string(snr,format='(F6.1)'),2)
    if keyword_set(fitvsini) then out=out+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+$
              '+/-'+strtrim(string(vsinierr,format='(F9.1)'),2)+' km/s'
    if keyword_set(fixvsini) then out=out+' Fixed vsini='+strtrim(string(fixvsini,format='(F7.1)'),2)
    xyouts,mean([w0,w1]),0.15,out,align=0.5,charsize=1.2
  endif else begin
    out = 'Teff='+strtrim(string(bestpars[0],format='(I5)'),2)+' Logg='+string(bestpars[1],format='(F4.2)')+' '+$
          '[M/H]='+strtrim(string(bestpars[2],format='(F5.2)'),2)+' [alpha/Fe]='+strtrim(string(bestpars[3],format='(F5.2)'),2)+$
          ' Vel='+strtrim(string(vrel,format='(F9.1)'),2)+' R.ChiSq='+strtrim(string(rchisq,format='(F6.3)'),2)+$
          ' S/N='+strtrim(string(snr,format='(F6.1)'),2)
   if keyword_set(fitvsini) then out=out+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+' km/s'
   if keyword_set(fixvsini) then out=out+' Fixed vsini='+strtrim(string(fixvsini,format='(F7.1)'),2)+' km/s'
   xyouts,mean([w0,w1]),0.2,out,align=0.5,charsize=1.2
  endelse

  if keyword_set(psfile) then ps_close

endif

dt = systime(1)-t0
;print,'dt=',strtrim(dt,2)

;stop

if keyword_set(stp) then stop

end
