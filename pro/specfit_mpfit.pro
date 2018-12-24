;+
;
; SPECFIT_MPFIT
;
; This program finds the best fitting synthetic spectrum to an observed spectrum
; using Levenberg-Marquardt least-squares minimization with MPFIT.PRO
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
;  specpars   The spectrum parameters, specpars = [w0, w1, dw, res]
;  =mask      The mask to use with wave/spec.  1-for points to use,
;               0-for points to ignore.
;               This can also be an array of weights.
;  =initpars  Initial guess for the stellar parameters.  If not input, a rough grid of
;               Teff and logg will be searched.
;  /interp . Interpolate between grid points if necessary.
;  /nearest  If the stellar parameters are between grid points then
;              use the closest one instead of interpolating.
;  /zerovel   Do not solve for velocity, should already be set to rest frame.
;  /fitvsini  Fit vsini.  The synthetic spectra must *NOT* be prepared.
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
;     alphaarr   The 1D array of ]alpha/Fe] values searched
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>specfit_genetic,wave,spec,errspec,synstr,[4000.,4500.,1.4,2.7],fitstr,dir=specdir
;
; By D.Nidever  Nov 2008
;-

;###########################################################################

pro specfit_mpfit_dummy
forward_function specfit_mpfit_par2ind, specfit_mpfit_shift
end

;###########################################################################

function specfit_mpfit_par2ind,pars,synstr

; This gets the index for a given parameter set

npars = n_elements(pars)
nsynstr = n_elements(synstr)

if npars lt 4 or nsynstr eq 0 then begin
  print,'syntax - ind = specfit_mpfit_par2ind(pars,synstr)'
  return,-1
endif

; Find the closest grid point
teffind = first_el(minloc(abs(synstr.teff-pars[0])))
loggind = first_el(minloc(abs(synstr.logg-pars[1])))
metalind = first_el(minloc(abs(synstr.metal-pars[2])))
alphaind = first_el(minloc(abs(synstr.alpha-pars[3])))

ind = [teffind, loggind, metalind, alphaind]

return,ind

end

;###########################################################################

function specfit_mpfit_shift,indin,synstr,shift,parin,parout

; This returns the parameters and indices for certains shifts in the
; the parameters
;
; INPuTS:
;  indin   Input parameter indices
;  synstr  Synthetic spectrum grid structure
;  shift   The shift index
;
; OUTPUTS:
;  indout  The shifted parameter indices. -1 if outside the grid
;  parin   The input parameters
;  parout  The output parameters

nindin = n_elements(indin)
nsynstr = n_elements(synstr)
nshift = n_elements(shift)

undefine,indout,parin,parout

if nindin lt 4 or nsynstr eq 0 or nshift lt 4 then begin
  print,'syntax - specfit_mpfit_shift,indin,synstr,shift,indout,parin,parout'
  return,-1
endif

indout = indin + shift
parin = [ synstr.teff[indin[0]], synstr.logg[indin[1]], synstr.metal[indin[2]], synstr.alpha[indin[3]] ]
parout = -1

; Teff
nteff = n_elements(synstr.teff)
if indout[0] lt 0 then return,-1
if indout[0] gt (nteff-1) then return,-1

; logg
nlogg = n_elements(synstr.logg)
if indout[1] lt 0 then return,-1
if indout[1] gt (nlogg-1) then return,-1

; Metal
nmetal = n_elements(synstr.metal)
if indout[2] lt 0 then return,-1
if indout[2] gt (nmetal-1) then return,-1

; alpha
nalpha = n_elements(synstr.alpha)
if indout[3] lt 0 then return,-1
if indout[3] gt (nalpha-1) then return,-1

parout = [ synstr.teff[indout[0]], synstr.logg[indout[1]], synstr.metal[indout[2]], synstr.alpha[indout[3]] ]

return,indout

end

;###########################################################################

pro specfit_mpfit,wave0,spec0,errspec0,synstr,specpars,fitstr,mask=mask0,$
                  zerovel=zerovel,error=error,stp=stp,plot=pl,silent=silent,$
                  initpars=initpars0,nearest=nearest,interp=interp,fitvsini=fitvsini,$
                  fixvsini=fixvsini

;COMMON specfit_mpfit,mpstr,chigrid
COMMON specfit_mpfit,chigrid

undefine,error,chisq,vrel

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)
nsynstr = n_elements(synstr)
nspecpars = n_elements(specpars)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 or nsynstr eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_genetic,wave,spec,errspec,synstr,specpars,fitstr,mask=mask,'
  print,'                      fixvsini=fixvsini,zerovel=zerovel,plot=plot,error=error,stp=stp'
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

tags = tag_names(synstr)
indteff = max(strpos(tags,'TEFF'))
indlogg = max(strpos(tags,'LOGG'))
indmetal = max(strpos(tags,'METAL'))
indalpha = max(strpos(tags,'ALPHA'))
if tag_exist(synstr,'TEFF') eq 0 or tag_exist(synstr,'LOGG') eq 0 or tag_exist(synstr,'METAL') eq 0 or $
  tag_exist(synstr,'ALPHA') eq 0 or tag_exist(synstr,'FILE') eq 0 then begin
  error = 'SYNSTR must have TEFF, LOGG, METAL, ALPHA and FILE tags'
  print,error
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


;------------------------------
; DETERMINE FIRST GUESS
;------------------------------
if n_elements(initpars0) eq 0 then begin

  ntry = 5
  ; Try various Teff/logg values
  dteff = range(synstr.teff)/(ntry+1.0)
  teffarr = (dindgen(ntry)+1.0d0)*dteff + min(synstr.teff)
  dlogg = range(synstr.logg)/(ntry+1.0)
  loggarr = (dindgen(ntry)+1.0d0)*dlogg + min(synstr.logg)

  chisqarr = fltarr(ntry^2)
  teffvals = fltarr(ntry^2)
  loggvals = fltarr(ntry^2)

  ; Loop through Teff
  count=0
  for i=0,ntry-1 do begin
    ; Loop through logg
    for j=0,ntry-1 do begin

      trypars = [teffarr[i], loggarr[j], -0.50d0, 0.0d0]
      SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,trypars,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec2,mask=mask2,$
                             zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,fixvsini=fixvsini,error=error,/silent

      if n_elements(error) gt 0 then chisq=1e30
      chisqarr[count] = chisq
      teffvals[count] = teffarr[i]
      loggvals[count] = loggarr[j]
      count++

    end
  end

  ; Get best Teff/logg values
  best_ind = first_el(minloc(chisqarr))
  initpars = [teffvals[best_ind], loggvals[best_ind], -0.50d0, 0.0d0]

endif else initpars=initpars0



;------------------------------
; RUN MPFIT
;------------------------------
t0 = systime(1)

; Function to call
func = 'specfit_mpfit_func'

parinfo = REPLICATE({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0],step:0.03},4)
;parinfo = REPLICATE({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0],relstep:0.03},4)
;parinfo = REPLICATE({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0],relstep:0.05},4)
;parinfo = REPLICATE({fixed:0,limited:[0,0],limits:[0.0d0,0.0d0],step:0.0},4)
;parinfo = REPLICATE({fixed:0,limited:[1,1],limits:[0.0d0,0.0d0],relstep:0.10},4)
;parinfo = REPLICATE({fixed:0,limited:[1,1],limits:[0.0d0,0.0d0],relstep:0.02},4)
;parinfo[0].limits = minmax(synstr.teff)
;parinfo[1].limits = minmax(synstr.logg)
;parinfo[2].limits = minmax(synstr.metal)
;parinfo[0].step = 250.
;parinfo[1].step = 0.2
;parinfo[2].step = 0.1
if range(synstr.teff) eq 0.0 then parinfo[0].fixed = 1
if range(synstr.logg) eq 0.0 then parinfo[1].fixed = 1
if range(synstr.metal) eq 0.0 then parinfo[2].fixed = 1
if range(synstr.alpha) eq 0.0 then parinfo[3].fixed = 1
parinfo[0].step = 250
parinfo[0].step = 0.2
parinfo[0].step = 0.1
parinfo[0].step = 0

; Function arguments
functargs = {wave:wave2,spec:spec2,errspec:errspec2,mask:mask2,$
             synstr:synstr,specpars:specpars,zerovel:keyword_set(zerovel),$
             nearest:keyword_set(nearest),interp:keyword_set(interp),fitvsini:keyword_set(fitvsini),$
             fixvsini:0.0}
if n_elements(fixvsini) gt 0 then functargs.fixvsini=fixvsini[0]
;             par:dblarr(1e4,4),chisq:dblarr(1e4),vrel:dblarr(1e4),vrelerr:dblarr(1e4),count:0L}
;mpstr = {par:dblarr(1e4,4),chisq:dblarr(1e4),vrel:dblarr(1e4),vrelerr:dblarr(1e4),count:0L}

;par = [15000.0d0, 3.50d0, -0.5d0, 0.0d0]

maxiter = 20

;fpar = MPFITFUN(func,wave2,spec2,errspec2,initpars,bestnorm=chisq,dof=dof,functargs=functargs,$
;                perror=perror,status=status,yfit=yfit,parinfo=parinfo,maxiter=maxiter,$
;                nfev=nfev,niter=niter) ;,/quiet)



; I might need two methods: one for gird, one for /interp

; Keep track of which grid points we've already evaluated
chigrid = fltarr(n_elements(synstr.teff),n_elements(synstr.logg),n_elements(synstr.metal),$
                 n_elements(synstr.alpha))

par = initpars
ind = specfit_mpfit_par2ind(par,synstr)
count=0
flag=0
WHILE (flag eq 0) do begin

  ; Get chisq in neighboring pixels
  specfit_mpfit_neighbors,par,initchi,chisqarr,steparr,_Extra=functargs
  minchi = min(chisqarr)
  best = first_el(minloc(chisqarr))

  ; Calculate the gradient
  if minchi lt initchi then begin

    ; [npars, 2]
    best2 = array_indices(chisqarr,best)
    shft = lonarr(4)
    if best2[1] eq 0 then shft[best2[0]]=-1 else shft[best2[0]]=1
    indout = specfit_mpfit_shift(ind,synstr,shft,parin,parout)

    ; New parameters
    par = parout
    ind = indout
    chisq = minchi

  ; Nearest neighbors NOT better
  endif else begin

    ; Get chisq for ALL neighboring pixels
    specfit_mpfit_neighbors,par,dum,chisqarr_all,/all,_Extra=functargs
    minchi_all = min(chisqarr_all)
    best_all = first_el(minloc(chisqarr_all))
    best2_all = array_indices(chisqarr_all,best_all)
    
    ; Found a better solution
    ; 3*3*3*3=81 elements, the center element is 40
    if (best_all ne 40) then begin
    
      shft = best2_all-1
      indout = specfit_mpfit_shift(ind,synstr,shft,parin,parout)

      ; New parameters
      par = parout
      ind = indout
      chisq = minchi_all

    ; None are better, FINISHED
    endif else begin
      flag=1
    endelse

    ;stop

  endelse

  ;stop

  count++
ENDWHILE

; Maybe refine the solution if /interp

niter = count
dum = where(chigrid gt 0.0,nfev)
dof = n_elements(spec2)-4-1
fpar = par

; Get the final best solution
undefine,chisq,vrel,error
SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,fpar,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec2,mask=mask2,$
                       zerovel=zerovel,nearest=nearest,fitvsini=fitvsini,fixvsini=fixvsini,interp=interp,error=error,$
                       /silent,wsynth=wsynth,ssynth=ssynth

dt = systime(1)-t0

rchisq = chisq/dof
best_teff = fpar[0]
best_logg = fpar[1]
best_metal = fpar[2]
best_alpha = fpar[3]

best_vrel = vrel
;best_vrel = mpstr.vrel[nfev-1]
;vrelerr = mpstr.vrelerr[nfev-1]

bestpars = fpar
; CALCULATE PARAMETER ERRORS!
;PARERRORS = PERROR * SQRT(chisq / DOF)   ; scaled uncertainties
;tefferr = parerrors[0]
;loggerr = parerrors[1]
;metalerr = parerrors[2]
;alphaerr = parerrors[3]

; Output the best fit parameters
if not keyword_set(silent) then begin
  print,'Best Teff = ',string(best_teff,format='(I5)')+' K'
  print,'Best Logg = ',string(best_logg,format='(F4.2)')
  print,'Best [M/H] = ',strtrim(string(best_metal,format='(F5.2)'),2)
  print,'Best [alpha/Fe] = ',strtrim(string(best_alpha,format='(F5.2)'),2)
  print,'Best R.ChiSq = ',string(rchisq,format='(F7.4)')
  print,'Best Vrel = ',strtrim(string(best_vrel,format='(F9.3)'),2)+' km/s'
  if keyword_set(fitvsini) and not keyword_set(fixvsini) then print,'Best vsini = ',strtrim(string(vsini,format='(F7.1)'),2)
  if keyword_set(fixvsini) then print,'Fixed vsini = ',strtrim(string(vsini,format='(F7.1)'),2)
endif


; Put everything in a structure
; specpars, bestpars, rchisq, rvel, rvelerr, chisqarr, rvelarr, rvelerrarr
; best model name, mask used, dof, # of points
if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
npts = nwave
nused = total(mask)
if keyword_set(zerovel) then fixzero=1 else fixzero=0
;fitstr = {specpars:specpars,zerovel:fixzero,mask:maskinput,npts:npts,nused:nused,bestpars:bestpars,parerrors:parerrors,$
;          teff:best_teff,tefferr:tefferr,logg:best_logg,loggerr:loggerr,metal:best_metal,metalerr:metalerr,$ 
;          alpha:best_alpha,alphaerr:alphaerr,chisq:chisq,dof:dof,rchisq:rchisq,vrel:best_vrel,vrelerr:vrelerr,$
fitstr = {specpars:specpars,zerovel:fixzero,mask:maskinput,npts:npts,nused:nused,bestpars:bestpars,$
          teff:best_teff,logg:best_logg,metal:best_metal,$ 
          alpha:best_alpha,chisq:chisq,dof:dof,rchisq:rchisq,vrel:best_vrel,vrelerr:vrelerr,$
          vsini:vsini,niter:niter,nfev:nfev}


; Plot
if keyword_set(pl) then begin

  ;wset,2
  ;gg = where(chisqarr lt 500.)
  ;maxchisq = max(chisqarr[gg])
  ;displayc,reform(chisqarr),teffarr,loggarr,/log,max=maxchisq+1
  ;oplot,[best_teff],[best_logg],ps=1,co=250

  ;bestmodel = mksynthname(best_teff,best_logg,best_metal,best_alpha)
  ;SPECFIT_GETSYNSPEC,best_teff,best_logg,best_metal,best_alpha,wsyn,ssyn,dir=dir
  ;SPECFIT_MODELINTERP,synstr,best_teff,best_logg,best_metal,best_alpha,wsyn,ssyn
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
  oplot,wsynth,ssynth,co=250,linestyle=5
  oplot,wave2,spec2-ssynth
  oplot,[0,1e5],[0,0],linestyle=2

  ;oplot,ww2,ss,co=250,linestyle=5
  out = 'Teff='+string(best_teff,format='(I5)')+' Logg='+string(best_logg,format='(F4.2)')+' '+$
        '[M/H]='+strtrim(string(best_metal,format='(F5.2)'),2)+' [alpha/Fe]='+strtrim(string(best_alpha,format='(F5.2)'),2)+$
        ' R.ChiSq='+string(rchisq,format='(F6.4)')+' Vel='+strtrim(string(best_vrel,format='(F9.3)'),2)
  if keyword_set(fitvsini) or keyword_set(fixvsini) then out=out+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+' km/s'
  xyouts,mean(wave2),0.2,out,align=0.5,charsize=1.2

endif

;stop

if keyword_set(stp) then stop

end
