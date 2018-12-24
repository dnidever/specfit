pro specfit_hybrid,wave0,spec0,errspec0,synstr,specpars,fitstr,mask=mask0,$
                 zerovel=zerovel,error=error,stp=stp,plot=pl,silent=silent,$
                 fitvsini=fitvsini,fixvsini=fixvsini

;+
;
; SPECFIT_HYBRID
;
; This program fits synthetic spectra to an observed spectrum with
; SPECFIT_GENETIC.PRO and then SPECFIT_GRID.PRO to refine the solution.
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
;  IDL>specfit_hybrid,wave,spec,errspec,synstr,[4000.,4500.,1.4,2.7],fitstr,dir=specdir
;
; By D.Nidever  Nov 2008
;-

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
  print,'Syntax - specfit_hybrid,wave,spec,errspec,synstr,specpars,fitstr,mask=mask,'
  print,'                      zerovel=zerovel,fitvsini=fitvsini,fixvsini=fixvsini,'
  print,'                      plot=plot,error=error,stp=stp'
  error = 'Not enough inputs'
  return
endif

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


;; Prepare the observed spectrum and error spectrum
;;--------------------------------------------------
;SPECFIT_PREPSPEC,wave,spec,errspec0,mask,specpars,wave2,spec2,errspec2,mask2,/object
;nwave = n_elements(wave2)

wavein = wave
specin = spec
errin = errspec0

;#####################################################
; STEP 1:  SPECFIT_GENETIC to narrow down the region
;#####################################################
save_except = !EXCEPT & !EXCEPT=0   ; catch match errors

SPECFIT_GENETIC,wavein,specin,errin,synstr,specpars,gene_fitstr,/nearest,mask=mask,$
                zerovel=zerovel,fitvsini=fitvsini,fixvsini=fixvsini,/silent ;,plot=plot

dum = check_math()      ; catch math errors
!EXCEPT = save_except


;###############################################
; STEP 2:  SPECFIT_MPFIT to refine the solution
;###############################################

bestpars1 = gene_fitstr.bestpars

; Run mpfit
SPECFIT_MPFIT,wavein,specin,errin,synstr,specpars,mpfit_fitstr,mask=mask,zerovel=zerovel,$
              fitvsini=fitvsini,fixvsini=fixvsini,/silent,initpars=gene_fitstr.bestpars ;,plot=plot
bestpars2 = mpfit_fitstr.bestpars

; Final parameters
bestpars = bestpars2
best_teff = bestpars[0]
best_logg = bestpars[1]
best_metal = bestpars[2]
best_alpha = bestpars[3]
chisq = mpfit_fitstr.chisq
dof = mpfit_fitstr.dof
rchisq = mpfit_fitstr.rchisq
vrel = mpfit_fitstr.vrel
vrelerr = mpfit_fitstr.vrelerr
vsini = mpfit_fitstr.vsini
;bestpars = mpfit_fitstr.bestpars
npts = mpfit_fitstr.npts
nused = mpfit_fitstr.nused

;; Set up the grid to use
;grid = {teff:dblarr(3),logg:dblarr(3),metal:dblarr(3),alpha:dblarr(3)}
;
;; Teff
;nteffstep = 5
;indteff = where(abs(synstr.teff-bestpars1[0]) lt 1d-6,nindteff)
;nteffarr = n_elements(synstr.teff)
;loindteff = (indteff-nteffstep) > 0
;tefflo = synstr.teff[loindteff]
;hiindteff = (indteff+nteffstep) < (nteffarr-1)
;nteff = hiindteff-loindteff+1
;if nteffarr gt 1 then dteff = synstr.teff[1]-synstr.teff[0] else dteff=0.0
;grid.teff = [tefflo, dteff, nteff]
;
;; Logg
;nloggstep = 5
;indlogg = where(abs(synstr.logg-bestpars1[1]) lt 1d-6,nindlogg)
;nloggarr = n_elements(synstr.logg)
;loindlogg = (indlogg-nloggstep) > 0
;logglo = synstr.logg[loindlogg]
;hiindlogg = (indlogg+nloggstep) < (nloggarr-1)
;nlogg = hiindlogg-loindlogg+1
;if nloggarr gt 1 then dlogg = synstr.logg[1]-synstr.logg[0] else dlogg=0.0
;grid.logg = [logglo, dlogg, nlogg]
;
;; Metal
;nmetalstep = 5
;indmetal = where(abs(synstr.metal-bestpars1[2]) lt 1d-6,nindmetal)
;nmetalarr = n_elements(synstr.metal)
;loindmetal = (indmetal-nmetalstep) > 0
;metallo = synstr.metal[loindmetal]
;hiindmetal = (indmetal+nmetalstep) < (nmetalarr-1)
;nmetal = hiindmetal-loindmetal+1
;if nmetalarr gt 1 then dmetal = synstr.metal[1]-synstr.metal[0] else dmetal=0.0
;grid.metal = [metallo, dmetal, nmetal]
;
;; Alpha
;nalphastep = 5
;indalpha = where(abs(synstr.alpha-bestpars1[3]) lt 1d-6,nindalpha)
;nalphaarr = n_elements(synstr.alpha)
;loindalpha = (indalpha-nalphastep) > 0
;alphalo = synstr.alpha[loindalpha]
;hiindalpha = (indalpha+nalphastep) < (nalphaarr-1)
;nalpha = hiindalpha-loindalpha+1
;if nalphaarr gt 1 then dalpha = synstr.alpha[1]-synstr.alpha[0] else dalpha=0.0
;grid.alpha = [alphalo, dalpha, nalpha]
;
;
; Run the grid
;SPECFIT_GRID,wavein,specin,errin,synstr,grid,specpars,grid_fitstr,mask=mask,zerovel=zerovel,$
;             fitvsini=fitvsini,/silent ;,plot=plot
;bestpars2 = grid_fitstr.bestpars
;
;;; Final parameters
;;best_teff = grid_fitstr.bestpars[0]
;;best_logg = grid_fitstr.bestpars[1]
;;best_metal = grid_fitstr.bestpars[2]
;;best_alpha = grid_fitstr.bestpars[3]
;;chisq = grid_fitstr.chisq
;;dof = grid_fitstr.dof
;;rchisq = grid_fitstr.rchisq
;;vrel = grid_fitstr.vrel
;;vrelerr = grid_fitstr.vrelerr
;;bestpars = grid_fitstr.bestpars
;;npts = grid_fitstr.npts
;;nused = grid_fitstr.nused
;
;
;; Do a search along each axis
;;------------------------------
;gridbest = REPLICATE({specpars:dblarr(4),bestpars:dblarr(4),zerovel:0,mask:0,npts:0L,nused:0.,$
;                      chisq:0.0,dof:0.0,rchisq:0.0,vrel:0.0,vrelerr:0.0,vsini:0.0},3)
;for i=0,2 do begin
;
;  ; Set up the grid to use
;  grid1 = {teff:dblarr(3),logg:dblarr(3),metal:dblarr(3),alpha:dblarr(3)}
;
;  ; Teff
;  nteffstep = 1
;  indteff = where(abs(synstr.teff-bestpars2[0]) lt 1d-6,nindteff)
;  nteffarr = n_elements(synstr.teff)
;  loindteff = (indteff-nteffstep) > 0
;  tefflo = synstr.teff[loindteff]
;  hiindteff = (indteff+nteffstep) < (nteffarr-1)
;  nteff = hiindteff-loindteff+1
;  if nteffarr gt 1 then dteff = synstr.teff[1]-synstr.teff[0] else dteff=0.0
;  grid1.teff = [tefflo, dteff, nteff]
;  if i eq 0 then grid1.teff  = [min(synstr.teff),(synstr.teff[1]-synstr.teff[0]),nteffarr]
;
;  ; Logg
;  nloggstep = 1
;  indlogg = where(abs(synstr.logg-bestpars2[1]) lt 1d-6,nindlogg)
;  nloggarr = n_elements(synstr.logg)
;  loindlogg = (indlogg-nloggstep) > 0
;  logglo = synstr.logg[loindlogg]
;  hiindlogg = (indlogg+nloggstep) < (nloggarr-1)
;  nlogg = hiindlogg-loindlogg+1
;  if nloggarr gt 1 then dlogg = synstr.logg[1]-synstr.logg[0] else dlogg=0.0
;  grid1.logg = [logglo, dlogg, nlogg]
;  if i eq 1 then grid1.logg  = [min(synstr.logg),(synstr.logg[1]-synstr.logg[0]),nloggarr]
;
;  ; Metal
;  nmetalstep = 1
;  indmetal = where(abs(synstr.metal-bestpars2[2]) lt 1d-6,nindmetal)
;  nmetalarr = n_elements(synstr.metal)
;  loindmetal = (indmetal-nmetalstep) > 0
;  metallo = synstr.metal[loindmetal]
;  hiindmetal = (indmetal+nmetalstep) < (nmetalarr-1)
;  nmetal = hiindmetal-loindmetal+1
;  if nmetalarr gt 1 then dmetal = synstr.metal[1]-synstr.metal[0] else dmetal=0.0
;  grid1.metal = [metallo, dmetal, nmetal]
;  if i eq 2 then grid1.metal  = [min(synstr.metal),(synstr.metal[1]-synstr.metal[0]),nmetalarr]
;
;  ; Alpha
;  grid1.alpha = [bestpars2[3],0.0,1.0]
;
;
;  ; Run the grid
;  SPECFIT_GRID,wavein,specin,errin,synstr,grid1,specpars,grid1_fitstr,mask=mask,zerovel=zerovel,$
;               fitvsini=fitvsini,/silent ;,plot=plot
;
;  ; Best parameters
;  temp = gridbest[i]
;  STRUCT_ASSIGN,grid1_fitstr,temp
;  gridbest[i] = temp
;
;  ;stop
;
;end
;
;
;; Final parameters
;bestind = first_el(minloc(gridbest.rchisq))
;fgrid = gridbest[bestind[0]]
;best_teff = fgrid.bestpars[0]
;best_logg = fgrid.bestpars[1]
;best_metal = fgrid.bestpars[2]
;best_alpha = fgrid.bestpars[3]
;chisq = fgrid.chisq
;dof = fgrid.dof
;rchisq = fgrid.rchisq
;vrel = fgrid.vrel
;vrelerr = fgrid.vrelerr
;vsini = fgrid.vsini
;bestpars = fgrid.bestpars
;npts = fgrid.npts
;nused = fgrid.nused
;
;;; Final parameters
;;best_teff = grid_fitstr.bestpars[0]
;;best_logg = grid_fitstr.bestpars[1]
;;best_metal = grid_fitstr.bestpars[2]
;;best_alpha = grid_fitstr.bestpars[3]
;;chisq = grid_fitstr.chisq
;;dof = grid_fitstr.dof
;;rchisq = grid_fitstr.rchisq
;;vrel = grid_fitstr.vrel
;;vrelerr = grid_fitstr.vrelerr
;;bestpars = grid_fitstr.bestpars
;;npts = grid_fitstr.npts
;;nused = grid_fitstr.nused

; Output the best fit parameters
if not keyword_set(silent) then begin
  print,'Best Teff = ',string(best_teff,format='(I5)')+' K'
  print,'Best Logg = ',string(best_logg,format='(F4.2)')
  print,'Best [M/H] = ',strtrim(string(best_metal,format='(F5.2)'),2)
  print,'Best [alpha/Fe] = ',strtrim(string(best_alpha,format='(F5.2)'),2)
  print,'Best R.ChiSq = ',string(rchisq,format='(F7.4)')
  print,'Best Vrel = ',strtrim(string(vrel,format='(F9.3)'),2)+' +/- '+strtrim(string(vrelerr,format='(F9.3)'),2)+' km/s'
  if keyword_set(fitvsini) then print,'Best vsini = '+strtrim(string(vsini,format='(F7.1)'),2)
  if keyword_set(fixvsini) then print,'Fixed vsini = '+strtrim(string(fixvsini,format='(F7.1)'),2)
endif


; Make final FITSTR
if n_elements(mask0) gt 0 then maskinput=1 else maskinput=0
if keyword_set(zerovel) then fixzero=1 else fixzero=0
fitstr = {specpars:specpars,zerovel:fixzero,mask:maskinput,npts:npts,nused:nused,bestpars:bestpars,$
          teff:best_teff,logg:best_logg,metal:best_metal,alpha:best_alpha,chisq:chisq,dof:dof,$
          rchisq:rchisq,vrel:vrel,vrelerr:vrelerr,vsini:vsini}


; Plot
if keyword_set(pl) then begin

  ; Get prepared object spectrum
  SPECFIT_PREPSPEC,wave,spec,errspec0,mask,specpars,wave2,spec2,errspec2,mask2,/object
  ; Get the best fitting, velocity-shifted synthetic spectrum
  SPECFIT_COMPARESYNSPEC,wave2,spec2,synstr,bestpars,specpars,errspec=errspec2,mask=mask2,$
                         zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent,wsynth=wsynth,ssynth=ssynth

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
        ' R.ChiSq='+string(rchisq,format='(F6.4)')+' Vel='+strtrim(string(vrel,format='(F9.3)'),2)
  if keyword_set(fitvsini) then out=out+' vsini='+strtrim(string(vsini,format='(F7.1)'),2)+' km/s'
  if keyword_set(fitxsini) then out=out+' Fixed vsini='+strtrim(string(fixvsini,format='(F7.1)'),2)+' km/s'
  xyouts,mean(wave2),0.2,out,align=0.5,charsize=1.2

endif

dt = systime(1)-t0
;print,'dt=',strtrim(dt,2)

;stop

if keyword_set(stp) then stop

end
