;+
;
; SPECFIT_GENETIC
;
; This program finds the best fitting synthetic spectrum to an observed spectrum
; using a genetic algorithm.
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
;     alphaarr   The 1D array of ]alpha/Fe] values searched
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>specfit_genetic,wave,spec,errspec,synstr,[4000.,4500.,1.4,2.7],fitstr,dir=specdir
;
; By D.Nidever  Nov 2008
;-

;pro specfit_genetic_dummy
;
;end

;###########################################################################
;
;function specfit_genetic_func,par,_Extra=extra
;
;; This is the function that the genetic algorithm calls
;
;;COMMON specfit,wave2,spec2,errspec2,mask2,specpars,zerovel
;
;synpars = reform(par)
;
;wave = extra.wave
;spec = extra.spec
;errspec = extra.errspec
;mask = extra.mask
;synstr = extra.synstr
;specpars = extra.specpars
;zerovel = extra.zerovel
;
;undefine,chisq,vrel,error
;SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,$
;                       zerovel=zerovel,error=error,/silent
;
;if n_elements(error) gt 0 then begin
;  chisq = 999999.
;  vrel = 999999.
;  vrelerr = 999999.
;endif
;
;;print,extra.count,synpars,chisq,vrel
;
; Keep track of the outputs
;extra.par[extra.count,*] = synpars
;extra.chisq[extra.count] = chisq
;extra.vrel[extra.count] = vrel
;extra.vrelerr[extra.count] = vrelerr
;extra.count++
;
;;stop
;
;return,chisq
;
;end
;
;###########################################################################

pro specfit_genetic,wave0,spec0,errspec0,synstr,specpars,fitstr,mask=mask0,$
                    inter=interp,nearest=nearest,zerovel=zerovel,error=error,$
                    fitvsini=fitvsini,fixvsini=fixvsini,stp=stp,plot=pl,silent=silent

;COMMON specfit,wave2,spec2,errspec2,mask2,specpars,zerovel
COMMON specfit,synstr2,functargs,genestr

undefine,error,chisq,vrel

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)
nsynstr = n_elements(synstr)
nspecpars = n_elements(specpars)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 or nsynstr eq 0 or nspecpars eq 0 then begin
  print,'Syntax - specfit_genetic,wave,spec,errspec,synstr,specpars,fitstr,mask=mask,'
  print,'                      inter=interp,nearest=nearest,zerovel=zerovel,fitvsini=fitvsini,'
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
  error = 'Not enough spectrum parameters.  Need specpars=[w0,w1,dw,res]'
  print,error
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
; RUN the GENETIC ALGORITHM
;------------------------------

; The parameter space to be explored
prange = fltarr(2,4)
prange[*,0] = minmax(synstr.teff)
prange[*,1] = minmax(synstr.logg)
prange[*,2] = minmax(synstr.metal)
prange[*,3] = minmax(synstr.alpha)

; Function to call
func = 'specfit_genetic_func'

; Function arguments
functargs = {wave:wave2,spec:spec2,errspec:errspec2,mask:mask2,$
             specpars:specpars,zerovel:keyword_set(zerovel),nearest:keyword_set(nearest),$
            ; synstr:synstr,specpars:specpars,zerovel:keyword_set(zerovel),nearest:keyword_set(nearest),$
             interp:keyword_set(interp),fitvsini:keyword_set(fitvsini),fixvsini:0.0}
if n_elements(fixvsini) gt 0 then functargs.fixvsini=fixvsini[0]
;             par:dblarr(1e4,4),chisq:dblarr(1e4),vrel:dblarr(1e4),vrelerr:dblarr(1e4),count:0L}
;genestr = {par:dblarr(1e4,4),chisq:dblarr(1e4),vrel:dblarr(1e4),vrelerr:dblarr(1e4),timearr:dblarr(1e4),count:0L}
synstr2 = synstr  ; for the common block


ftol = 1.e-2 
t0 = systime(1)
pcross = 0.90 ;0.95
gene_length = 25
pmutate = 0.01 ;0.05 ;0.01
itmax = 20 ;15 ;30  ;10
npop = 200 ;150 ;100 ;200 ;50
fpar = RMD_GA(    ftol,                               $ 
                  function_value = function_value,    $ 
                  function_name = func,               $ 
               ;   functargs = functargs,              $
                  prange = prange,                    $ 
                  ncalls = ncalls,                    $ 
                  quiet = quiet,                      $ 
                  iterproc = iterproc,                $ 
                  iterargs = iterargs,                $ 
                  pcross = pcross,                    $ 
                  gene_length = gene_length,          $ 
                  pmutate = pmutate,                  $ 
                  itmax = itmax,                      $ 
                  npop = npop) 

;  PCROSS:              Probability that two genes will be crossed (default: 0.9)
;  PMUTATE:             Probability that a single gene will have one of its
;                       chromosomes bit-flipped (default:0.01)
;  ITMAX:               Maximum number of generations to proceed in the evolution (default:10)
;  NPOP:                Number of chromosomes in the population (default: 20)
;  GENE_LENGTH:         Length in bits of each member of the population.  The more
;                       bits the better the numerical precision of the result (default: 5)
;  STRETCH_FACTOR:      Amount by which to *stretch* the range of accepted fitness
;                       values.


; Final results
chisq = function_value
dof = total(mask)
if not keyword_set(zerovel) then dof=dof-1    ; fit one parameter
rchisq = chisq/dof
best_teff = fpar[0]
best_logg = fpar[1]
best_metal = fpar[2]
best_alpha = fpar[3]
;best_vrel = genestr.vrel[genestr.count-1]
;vrelerr = genestr.vrelerr[genestr.count-1]
bestpars = reform(fpar)

; Get best model
undefine,chisq,vrel,error
SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,fpar,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec2,$
                       /nearest,zerovel=zerovel,fitvsini=fitvsini,fixvsini=fixvsini,error=error,/silent,$
                       wsynth=wsynth,ssynth=ssynth,synhead=synhead
best_vrel = vrel

; If NEAREST find the right parameters
if keyword_set(nearest) then begin
  hdteff = SXPAR(synhead,'TEFF',count=nhdteff)
  if nhdteff gt 0 then best_teff = hdteff
  hdlogg = SXPAR(synhead,'LOGG',count=nhdlogg)
  if nhdlogg gt 0 then best_logg = hdlogg
  hdmetal = SXPAR(synhead,'FEH',count=nhdmetal)
  if nhdmetal gt 0 then best_metal = hdmetal
  hdalpha = SXPAR(synhead,'ALPHA_FE',count=nhdalpha)
  if nhdalpha gt 0 then best_alpha = hdalpha

  bestpars = [best_teff, best_logg, best_metal, best_alpha]
endif


; Output the best fit parameters
if not keyword_set(silent) then begin
  print,'Best Teff = ',string(best_teff,format='(I5)')+' K'
  print,'Best Logg = ',string(best_logg,format='(F4.2)')
  print,'Best [M/H] = ',strtrim(string(best_metal,format='(F5.2)'),2)
  print,'Best [alpha/Fe] = ',strtrim(string(best_alpha,format='(F5.2)'),2)
  print,'Best R.ChiSq = ',string(rchisq,format='(F7.4)')
  print,'Best Vrel = ',strtrim(string(best_vrel,format='(F9.3)'),2)+' km/s'
  if keyword_set(fitvsini) then print,'Best vsini = ',strtrim(string(vsini,format='(F7.1)'),2)
  if keyword_set(fixvsini) then print,'Fixed vsini = ',strtrim(string(fixvsini,format='(F7.1)'),2)
endif


; Put everything in a structure
; specpars, bestpars, rchisq, rvel, rvelerr, chisqarr, rvelarr, rvelerrarr
; best model name, mask used, dof, # of points
if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
npts = nwave
nused = total(mask)
;bestmodel = mksynthname(best_teff,best_logg,best_metal,best_alpha)
if keyword_set(zerovel) then fixzero=1 else fixzero=0
fitstr = {specpars:specpars,zerovel:fixzero,mask:maskinput,npts:npts,nused:nused,bestpars:bestpars,$
          teff:best_teff,logg:best_logg,metal:best_metal,alpha:best_alpha,chisq:chisq,dof:dof,$
          rchisq:rchisq,vrel:best_vrel,vrelerr:vrelerr,vsini:vsini,pcross:pcross,$
          gene_length:gene_length,pmutate:pmutate,itmax:itmax,npop:npop,ncalls:ncalls}

; Plot
if keyword_set(pl) then begin

  ; Get the best fitting, velocity-shifted synthetic spectrum
  ;SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,bestpars,specpars,errspec=errspec2,mask=mask2,$
  ;                       zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent,wsynth=wsynth,ssynth=ssynth

  ;wset,2
  ;gg = where(chisqarr lt 500.)
  ;maxchisq = max(chisqarr[gg])
  ;displayc,reform(chisqarr),teffarr,loggarr,/log,max=maxchisq+1
  ;oplot,[best_teff],[best_logg],ps=1,co=250

  ;bestmodel = mksynthname(best_teff,best_logg,best_metal,best_alpha)
  ;SPECFIT_GETSYNSPEC,best_teff,best_logg,best_metal,best_alpha,wsyn,ssyn,dir=dir
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
        '[M/H]='+strtrim(string(best_metal,format='(F5.2)'),2)+' [alpha/Fe]='+strtrim(string(best_alpha,format='(F5.2)'),2)+$
        ' R.ChiSq='+string(rchisq,format='(F6.4)')+' Vel='+strtrim(string(best_vrel,format='(F9.3)'),2)
  if keyword_set(fitvsini) then out=out+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+' km/s'
  if keyword_set(fixvsini) then out=out+' Fixed vsini='+strtrim(string(fixvsini,format='(F7.1)'),2)+' km/s'
  xyouts,mean(wave2),0.2,out,align=0.5,charsize=1.2

endif

dt = systime(1)-t0
;print,'dt=',strtrim(dt,2)

;stop

if keyword_set(stp) then stop

end
