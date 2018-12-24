pro specfit_grid,wave0,spec0,errspec0,synstr,gridstr,specpars,fitstr,mask=mask0,$
                 zerovel=zerovel,nearest=nearest,interp=interp,error=error,$
                 fitvsini=fitvsini,fixvsini=fixvsini,stp=stp,plot=pl,silent=silent

;+
;
; SPECFIT_GRID
;
; This program compares an observed spectrum to a grid of synthetic spectra
;
; INPUTS:
;  wave       The observed wavelength array
;  spec       The observed spectrum array
;  errspec    The error spectrum array
;  synstr    The synthetic grid structure.  It should have TEFF, LOGG, METAL, ALPHA and FILE
;              tags.  TEFF, LOGG, METAL and ALPHA should be 1D arrays of the UNIQUE stellar
;              parameters of the synthetic grid.  FILE should be a 4D string array (or 3D if
;              ALPHA only has one element) of dimensions [Nteff,Nlogg,Nmetal,Nalpha] with
;              the absolute path of the synthetic spectrum FITS file.  It should be an empty
;              string if the spectrum does not exist.
;  gridstr    The grid structure of stellar parameters to search.  This should have
;               teff, logg, metal and alpha tags, each being a 3-element
;               array with minimum, step and nsteps.
;  specpars   The spectrum parameters, specpars = [w0, w1, dw, res]
;  =mask      The mask to use with wave/spec.  1-for points to use,
;               0-for points to ignore.
;               This can also be an array of weights.
;  /interp . Interpolate between grid points if necessary.
;  /nearest  If the stellar parameters are between grid points then
;              use the closest one instead of interpolating.
;  /zerovel   Do not solve for velocity, should already be set to rest frame.
;  /fitvsini  Fit vsini.
;  =fixvsini  Fix vsini at a certain value (in km/s).
;  /plot      Plot the spectra.
;  /silent    Don't print anything to the screen
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  fitstr     Structure that contains all of the output information, including:
;     bestpars   The array of best-fitting parameters, Teff, logg, metal and alpha.
;     chisq      The chi squared of the best fit.
;     vrel       The relative velocity of the best fit
;     vrelerr    The uncertainty of the relative velocity of the best fit
;     chisqarr   The entire chisq array of the grid
;     vrelarr    The entire vrel array of the grid
;     vrelerrarr The entire vrelerr array of the grid
;     teffarr    The 1D array of Teff values searched
;     loggarr    The 1D array of logg values searched
;     metalarr   The 1D array of [M/H] values searched
;     alphaarr   The 1D array of [alpha/Fe] values searched
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>specfit_grid,wave,spec,errspec,synstr,gridstr,[4000.,4500.,1.4,2.7],fitstr,dir=specdir
;
; By D.Nidever  Oct 2008
;-

undefine,error,chisq,vrel,chisqarr,vrelarr,vrelerrarr,teffarr,loggarr,metalarr,alphaarr

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)
nsynstr = n_elements(synstr)
ngridstr = n_elements(gridstr)
nspecpars = n_elements(specpars)
ndir = n_elements(dir)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 or nsynstr eq 0 or ngridstr eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_grid,wave,spec,errspec,synstr,gridstr,specpars,fitstr,mask=mask,'
  print,'                      zerovel=zerovel,inter=interp,nearest=nearest,fitvsini=fitvsini,'
  print,'                      fixvsini=fixvsini,plot=plot,error=error,stp=stp'
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



; Checking the grid structure
type = size(gridstr,/type)
if type ne 8 then begin
  error = 'GRIDSTR must be a structure'
  print,error
  return
endif

tags = tag_names(gridstr)
indteff = max(strpos(tags,'TEFF'))
indlogg = max(strpos(tags,'LOGG'))
indmetal = max(strpos(tags,'METAL'))
indalpha = max(strpos(tags,'ALPHA'))
if indteff[0] eq -1 or indlogg[0] eq -1 or indmetal[0] eq -1 or indalpha[0] eq -1 then begin
  error = 'GRIDSTR must have TEFF, LOGG, METAL and ALPHA tags'
  print,error
  return
endif

teffsz = size(gridstr.teff)
loggsz = size(gridstr.logg)
metalsz = size(gridstr.metal)
alphasz = size(gridstr.alpha)
if (teffsz[0] ne 1 or loggsz[0] ne 1 or metalsz[0] ne 1 or alphasz[0] ne 1) or $
   (teffsz[1] ne 3 or loggsz[1] ne 3 or metalsz[1] ne 3 or alphasz[1] ne 3)  then begin
  error = 'The TEFF, LOGG, METAL and ALPHA tags in GRIDSTR must be 3-element arrays: minimum, step, nsteps'
  print,error
  return
endif


; Getting grid parameters
teffpars = dblarr(3)
loggpars = dblarr(3)
metalpars = dblarr(3)
alphapars = dblarr(3)
for i=0,2 do begin

  tefftest = VALID_NUM(gridstr.teff[i],iteff)
  if tefftest eq 0 then begin
    error = 'TEFF['+strtrim(i,2)+']='+strtrim(gridstr.teff[i],2)+' is NOT A NUMBER'
    print,error
    return
  endif
  teffpars[i] = iteff

  loggtest = VALID_NUM(gridstr.logg[i],ilogg)
  if loggtest eq 0 then begin
    error = 'LOGG['+strtrim(i,2)+']='+strtrim(gridstr.logg[i],2)+' is NOT A NUMBER'
    print,error
    return
  endif
  loggpars[i] = ilogg

  metaltest = VALID_NUM(gridstr.metal[i],imetal)
  if metaltest eq 0 then begin
    error = 'METAL['+strtrim(i,2)+']='+strtrim(gridstr.metal[i],2)+' is NOT A NUMBER'
    print,error
    return
  endif
  metalpars[i] = imetal

  alphatest = VALID_NUM(gridstr.alpha[i],ialpha)
  if alphatest eq 0 then begin
    error = 'ALPHA['+strtrim(i,2)+']='+strtrim(gridstr.alpha[i],2)+' is NOT A NUMBER'
    print,error
    return
  endif
  alphapars[i] = ialpha

end

; Steps must be positive, unless Nsteps=1
if (teffpars[1] le 0.0 and teffpars[2] ne 1.0) or (loggpars[1] le 0.0 and loggpars[2] ne 1.0) or $
   (metalpars[1] le 0.0 and metalpars[2] ne 1.0) or (alphapars[1] le 0.0 and alphapars[2] ne 1.0) then begin
  print,'TEFF/LOGG/METAL/ALPHA steps MUST BE POSITIVE'
  error = 'TEFF/LOGG/METAL/ALPHA steps MUST BE POSITIVE'
  return
endif
; Nsteps must be integers
if (long(teffpars[2])-teffpars[2] ne 0.0) or (long(loggpars[2])-loggpars[2] ne 0.0) or $
   (long(metalpars[2])-metalpars[2] ne 0.0) or (long(alphapars[2])-alphapars[2] ne 0.0) then begin
  print,'TEFF/LOGG/METAL/ALPHA Nsteps MUST BE INTEGERS'
  error = 'TEFF/LOGG/METAL/ALPHA Nsteps MUST BE INTEGERS'
  return
endif
; Nsteps must be >=1
if teffpars[2] le 0.0 or loggpars[2] le 0.0 or metalpars[2] le 0.0 or alphapars[2] le 0.0 then begin
  print,'TEFF/LOGG/METAL/ALPHA Nsteps be >=1'
  error = 'TEFF/LOGG/METAL/ALPHA Nsteps be >=1'
  return
endif

; Not enough spectrum parameters
if nspecpars lt 4 then begin
  print,'Not enough spectrum parameters.  Need specpars=[w0,w1,dw,res]'
  error = 'Not enough spectrum parameters.  Need specpars=[w0,w1,dw,res]'
  return
endif


; Spectrum parameters
w0 = specpars[0]
w1 = specpars[1]
dw = specpars[2]
res = specpars[3:*]  ; might be an array


; Prepare the observed spectrum and error spectrum
;--------------------------------------------------
SPECFIT_PREPSPEC,wave,spec,errspec0,mask,specpars,wave2,spec2,errspec2,mask2,/object,/silent
nwave = n_elements(wave2)


; Stellar parameter pars are: minimum, stepsize, Nsteps
teffarr = dindgen(teffpars[2])*teffpars[1]+teffpars[0]
loggarr = dindgen(loggpars[2])*loggpars[1]+loggpars[0]
metalarr = dindgen(metalpars[2])*metalpars[1]+metalpars[0]
alphaarr = dindgen(alphapars[2])*alphapars[1]+alphapars[0]

; Setup the chisq and vrel grid
chisqarr = fltarr(teffpars[2],loggpars[2],metalpars[2],alphapars[2])+999999.
vrelarr = fltarr(teffpars[2],loggpars[2],metalpars[2],alphapars[2])+999999.
vrelerrarr = fltarr(teffpars[2],loggpars[2],metalpars[2],alphapars[2])+999999.
vsiniarr = fltarr(teffpars[2],loggpars[2],metalpars[2],alphapars[2])+999999.

t0 = systime(1)

; TEFF LOOP
;-------------------
FOR i=0,teffpars[2]-1 do BEGIN

  iteff = teffpars[0] + i*teffpars[1]

  ; LOGG LOOP
  ;---------------
  FOR j=0,loggpars[2]-1 do BEGIN

    ilogg = loggpars[0] + j*loggpars[1]

    ; METAL LOOP
    ;----------------
    FOR k=0,metalpars[2]-1 do BEGIN

      imetal = metalpars[0] + k*metalpars[1]

      ; ALPHA LOOP
      ;---------------
      FOR l=0,alphapars[2]-1 do BEGIN

        ialpha = alphapars[0] + l*alphapars[1]

        synpars = double([iteff,ilogg,imetal,ialpha])

        ; Compare the spectra
        undefine,chisq,vrel,error
        SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,synpars,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec2,$
                               mask=mask2,zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,$
                               fixvsini=fixvsini,error=error,/silent


        if n_elements(error) eq 0 then begin

          ; Put into the arrays
          chisqarr[i,j,k,l] = chisq
          vrelarr[i,j,k,l] = vrel
          vrelerrarr[i,j,k,l] = vrelerr
          vsiniarr[i,j,k,l] = vsini

          ;stop

        endif

        ;stop

      ENDFOR ; alpha loop
    ENDFOR ; metal loop
  ENDFOR ; logg loop
ENDFOR ; teff loop

dt = systime(1)-t0

; Find the best fitting synthetic spectrum

; Find best fit
best = first_el(minloc(chisqarr))
best2 = array_indices(chisqarr,best)
; one-element dimensions on the end get dropped
undefine,bestind_teff,bestind_logg,bestind_metal,bestind_alpha
if alphapars[2] eq 1 then bestind_alpha=0
if alphapars[2] eq 1 and metalpars[2] eq 1 then begin
  bestind_metal = 0
  bestind_alpha = 0
endif
if alphapars[2] eq 1 and metalpars[2] eq 1 and loggpars[2] eq 1 then begin
  bestind_logg = 0
  bestind_metal = 0
  bestind_alpha = 0
endif
bestind_teff = best2[0]
if n_elements(bestind_logg) eq 0 then bestind_logg = best2[1]
if n_elements(bestind_metal) eq 0 then bestind_metal = best2[2]
if n_elements(bestind_alpha) eq 0 then bestind_alpha = best2[3]

best_teff = teffarr[bestind_teff]
best_logg = loggarr[bestind_logg]
best_metal = metalarr[bestind_metal]
best_alpha = alphaarr[bestind_alpha]
best_chisq = min(chisqarr)
;best_chisq = chisqarr[bestind_teff,bestind_logg,bestind_metal,bestind_alpha]
best_vrel = vrelarr[bestind_teff,bestind_logg,bestind_metal,bestind_alpha]
best_vrelerr = vrelerrarr[bestind_teff,bestind_logg,bestind_metal,bestind_alpha]
best_vsini = vsiniarr[bestind_teff,bestind_logg,bestind_metal,bestind_alpha]
bestpars = [best_teff, best_logg, best_metal, best_alpha]
chisq = best_chisq
vrel = best_vrel
vrelerr = best_vrelerr

; Reduced ChiSq
dof = total(mask2)
if not keyword_set(zerovel) then dof=dof-1    ; fit one parameter
rchisq = chisq/dof

; Use REBIN to interpolate between points
; it uses bilinear interpolation
;sz = size(chisqarr)


; Error estimates
; Use the contour at min(chisq)+1, NOT REDUCED chisq
; Loop through stellar parameters, Teff, logg, metal, alpha
pars = fltarr(4,3)
pars[0,*] = teffpars
pars[1,*] = loggpars
pars[2,*] = metalpars
pars[3,*] = alphapars
parerrors = fltarr(4)
for i=0,3 do begin

  ipars = reform(pars[i,*])

  ; This parameter was fit
  if (ipars[2] gt 1) then begin

    ; Step through each parameter value and find the best ChiSq value
    ; in the rest of the chisq "cube"
    ichisqarr = fltarr(ipars[2])
    for j=0,ipars[2]-1 do begin

      case i of
      0: slice_chisq = chisqarr[j,*,*,*]
      1: slice_chisq = chisqarr[*,j,*,*]
      2: slice_chisq = chisqarr[*,*,j,*]
      3: slice_chisq = chisqarr[*,*,*,j]
      endcase

      best_slice_chisq = min(slice_chisq)
      ichisqarr[j] = best_slice_chisq

      ;stop

    end

    ; Spline
    pararr = findgen(ipars[2])*ipars[1]+ipars[0]
    pararr2 = scale_vector(findgen(ipars[2]*100.),min(pararr),max(pararr))
    ichisqarr2 = cspline(pararr,ichisqarr,pararr2)
    bestind = first_el(minloc(ichisqarr2))

    ; Find min(chisq)+1
    lo = (bestind[0]-100) > 0
    hi = (bestind[0]+100) < (n_elements(ichisqarr2)-1)
    coef = poly_fit(pararr2[lo:hi],ichisqarr2[lo:hi],2)
    ; coef = [c,b,a]
    ; (y-k)=a*(x-h)^2
    ;  k=4ac-b^2/4a  h=-b/2a (axis of symmetry)
    ;  The shape us determined by a alone.
    ;  Want x (relative to h) where y=k+1
    ;   y=a*x^2=1 -> x=1/sqrt(a)
    iparerror = 1.0/sqrt(coef[2])
    parerrors[i] = iparerror

  ; This parameter was NOT fit
  endif else begin
    parerrors[i] = 0
  endelse

  case i of
  0: tefferr = parerrors[i]
  1: loggerr = parerrors[i]
  2: metalerr = parerrors[i]
  3: alphaerr = parerrors[i]
  endcase

  ; THESE ERRORS SEEM **MUCH** TOO LOW!!

  ;stop

end

; Vrelerr estimate based on the global ChiSq distribution
; how much does Vrel vary within the 1-sigma error ellipse


; Output the best fit parameters
if not keyword_set(silent) then begin
  print,'Best Teff = ',string(best_teff,format='(I5)')+' +/- '+string(tefferr,format='(I5)')+' K'
  print,'Best Logg = ',string(best_logg,format='(F4.2)')+' +/- '+string(loggerr,format='(F6.4)')
  print,'Best [M/H] = ',strtrim(string(best_metal,format='(F5.2)'),2)+' +/- '+string(metalerr,format='(F5.2)')
  print,'Best [alpha/Fe] = ',strtrim(string(best_alpha,format='(F5.2)'),2)+' +/- '+string(alphaerr,format='(F5.2)')
  print,'Best R.ChiSq = ',string(rchisq,format='(F7.4)')
  print,'Best Vrel = ',strtrim(string(best_vrel,format='(F9.3)'),2)+' +/- '+strtrim(string(best_vrelerr,format='(F9.3)'),2)+' km/s'
  if keyword_set(fitvsini) then print,'Best vsini = ',strtrim(string(best_vsini,format='(F7.1)'),2)
  if keyword_set(fixvsini) then print,'Fixed vsini = ',strtrim(string(fixvsini,format='(F7.1)'),2)
endif


; Put everything in a structure
; specpars, bestpars, rchisq, rvel, rvelerr, chisqarr, rvelarr, rvelerrarr
; best model name, mask used, dof, # of points
if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
npts = nwave
nused = total(mask)
bestmodel = mksynthname(best_teff,best_logg,best_metal,best_alpha)
if keyword_set(zerovel) then fixzero=1 else fixzero=0
fitstr = {specpars:specpars,zerovel:fixzero,mask:maskinput,npts:npts,nused:nused,bestpars:bestpars,$
          parerrors:parerrors,teff:best_teff,tefferr:tefferr,logg:best_logg,loggerr:loggerr,metal:$
          best_metal,metalerr:metalerr,alpha:best_alpha,alphaerr:alphaerr,chisq:chisq,dof:dof,$
          rchisq:rchisq,vrel:vrel,vrelerr:vrelerr,vsini:best_vsini,bestmodel:bestmodel,$
          chisqarr:chisqarr,vrelarr:vrelarr,vrelerrarr:vrelerrarr,vsiniarr:vsiniarr,teffarr:teffarr,$
          loggarr:loggarr,metalarr:metalarr,alphaarr:alphaarr}

; Plot
if keyword_set(pl) then begin

  ; Get the best fitting, velocity-shifted synthetic spectrum
  SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,bestpars,specpars,errspec=errspec2,mask=mask2,zerovel=zerovel,$
                         fixvsini=fixvsini,fitvsini=fitvsini,nearest=nearest,interp=interp,error=error,$
                         /silent,wsynth=wsynth,ssynth=ssynth
  ;wset,2
  ;gg = where(chisqarr lt 500.)
  ;maxchisq = max(chisqarr[gg])
  ;displayc,reform(chisqarr),teffarr,loggarr,/log,max=maxchisq+1
  ;oplot,[best_teff],[best_logg],ps=1,co=250

  ;bestmodel = mksynthname(best_teff,best_logg,best_metal,best_alpha)
  ;SPECFIT_GETSYNSPEC,synstr,best_teff,best_logg,best_metal,best_alpha,wsyn,ssyn
  ;SPECFIT_PREPSPEC,wsyn,ssyn,0,0,specpars,ww,ss
  ;ww2 = ww*(best_vrel/cspeed+1.0)  ; shift

  plot,[0],[0],/nodata,xtit='Wavelength (A)',xr=[w0,w1],yr=[-0.15,1.2],xs=1,ys=1
  oplot,wave2,spec2
  if nmask gt 0 then begin
    binarymask2 = float(mask2 gt 0.0)
    breaks1 = where(slope([0,binarymask2]) eq 1,nbreaks1)
    breaks2 = where(slope([binarymask2,0]) eq -1,nbreaks2)
    ;breaks1 = where(slope([0,mask2]) eq 1,nbreaks1)
    ;breaks2 = where(slope([mask2,0]) eq -1,nbreaks2)
    for i=0,nbreaks1-1 do begin
      oplot,wave2[breaks1[i]:breaks2[i]],spec2[breaks1[i]:breaks2[i]],co=150
    end
    ;gd = where(mask2 gt 0.0,ngd)
    ;oplot,wave2[gd],spec2[gd],ps=1,co=150
  endif
  ;oplot,ww2,ss,co=250,linestyle=5
  oplot,wsynth,ssynth,co=250,linestyle=5
  oplot,wave2,spec2-ssynth
  oplot,[0,1e5],[0,0],linestyle=2
  out = 'Teff='+strtrim(string(best_teff,format='(I5)'),2)+' Logg='+string(best_logg,format='(F4.2)')+' '+$
        '[M/H]='+strtrim(string(best_metal,format='(F5.2)'),2)+' [alpha/Fe]='+strtrim(string(best_alpha,format='(F5.2)'),2)
  xyouts,mean(wave2),0.2,out,align=0.5,charsize=1.2
  out2 = 'R.ChiSq='+string(rchisq,format='(F6.4)')+' Vel='+strtrim(string(best_vrel,format='(F9.3)'),2)
  if keyword_set(fitvsini) then out=out+' vsini='+strtrim(string(best_vsini,format='(F7.1)'),2)+' km/s'
  if keyword_set(fixvsini) then out=out+' Fixed vsini='+strtrim(string(fixvsini,format='(F7.1)'),2)+' km/s'
  xyouts,mean(wave2),0.15,out2,align=0.5,charsize=1.2

endif

;stop

if keyword_set(stp) then stop

end
