pro prepare_synspec,synstr0,grid,specpars0,outdir

;+
;
; PREPARE_SYNSPEC
;
; This prepares a grid of synthetic spectra for a run of SPECFIT (fitting to
; observed spectra).  A rough grid of synthetic spectra (normalized) must
; already exist.
;
; These steps are taken:
; 1.) interpolate the synthethic spectra
; 2.) gaussian smooth to proper resolution
; 3.) set onto desires dispersion scale (logarithmic)
; 4.) normalize.  the same procedure as is used on the observed spectra
;        for consistency.
;
; A SYNSTR structure will be created for the prepared spectra and saved
; in "outdir" with the filename "synstr.dat". 
;
; INPUTS:
;  synstr    The structure describing the synthetic spectra grid (made by
;              MKSYNSTR.PRO)
;  grid      A structure with the grid parameters for the prepared spectral grid.
;              TEFF, LOGG, METAL and ALPHA tags with the 1D unique stellar parameter
;              values.
;  specpars  Structure with the spectrum parameters. W0, DW, NW, RES, LOG tags.
;  outdir    The directory to copy the prepared spectra to.
;
; OUTPUTS:
;  Prepared spectra (interpolated, smoothed, on proper dispersion scale, normalized).
;  A SYNSTR structure will be saved in "outdir".
;
; USAGE:
;  IDL>prepare_synspec,synstr0,grid,specpars,outdir
;
; By D.Nidever  Nov 2008
;- 

nsynstr = n_elements(synstr0)
ngrid = n_elements(grid)
nspecpars = n_elements(specpars0)
noutdir = n_elements(outdir)

if nsynstr eq 0 or ngrid eq 0 or nspecpars eq 0 or noutdir eq 0 then begin
  print,'prepare_synspec,synstr0,grid,specpars,outdir'
  return
endif

;;--------------------------------
;; ROUGH SYNTHETIC SPECTRAL GRID
;;-------------------------------
;restore,'/net/halo/dln5q/bstars/synspec/finegrid/finegrid_synstr.dat'
;
;;--------------------------------
;; PARAMETERS FOR "PREPARED" GRID
;;--------------------------------
;teff = dindgen(133)*250.0d0 + 7000.0d0
;logg = dindgen(26)*0.2d0 
;metal = [-2.0d0, -1.5d0, -1.0d0, -0.5d0, -0.25d0, 0.0d0]
;alpha = [0.0d0]
;grid = {teff:teff,logg:logg,metal:metal,alpha:alpha}
;;grid = {teff:dblarr(3),logg:dblarr(3),metal:dblarr(3),alpha:dblarr(3)}
;;grid.teff = [7000., 250., 133]  ; 7000-40000
;;grid.logg = [0.0, 0.2d0, 26]
;;grid.metal = [-1.0, 0.1d0, 11]
;;grid.alpha = [0.0, 0.0, 1]
;
;;specpars = [4000., 5200., 1.4, 2.7]
;specpars = {w0:3.5652573434202d0, dw:1.44903610837930d-4, nw:1070L, res:2.7d0, log:1}
;
;outdir = '/net/halo/dln5q/bstars/synspec/prepgrid/'

; Make the directory if necessary
if FILE_TEST(outdir,/directory) eq 0 then FILE_MKDIR,outdir

;-------------------------
; PREPARING THE SPECTRA
;-------------------------
;teffpars = grid.teff
;loggpars = grid.logg
;metalpars = grid.metal
;alphapars = grid.alpha

w0log = specpars0.w0   ; log(wave)
dwlog = specpars0.dw   ; dlog(wave)
nw = specpars0.nw
res = specpars0.res

lambda = (10.0d0)^(dindgen(nw)*dwlog+double(w0log))
specpars = [min(lambda),max(lambda),dwlog,res]


; Stellar parameter pars are: minimum, stepsize, Nsteps
;teffarr = dindgen(teffpars[2])*teffpars[1]+teffpars[0]
;loggarr = dindgen(loggpars[2])*loggpars[1]+loggpars[0]
;metalarr = dindgen(metalpars[2])*metalpars[1]+metalpars[0]
;alphaarr = dindgen(alphapars[2])*alphapars[1]+alphapars[0]
teffarr = grid.teff
loggarr = grid.logg
metalarr = grid.metal
alphaarr = grid.alpha

nteff = n_elements(teffarr)
nlogg = n_elements(loggarr)
nmetal = n_elements(metalarr)
nalpha = n_elements(alphaarr)

;nmodels = teffpars[2]*loggpars[2]*metalpars[2]*alphapars[2]
nmodels = nteff*nlogg*nmetal*nalpha
print,''
print,'PREPARING ',strtrim(string(nmodels,format='(I10)'),2),' SYNTHETIC SPECTRA'
print,''
wait,2


; TEFF LOOP
;-------------------
;FOR i=0,teffpars[2]-1 do BEGIN
FOR i=0,nteff-1 do BEGIN

  ;iteff = teffpars[0] + i*teffpars[1]
  iteff = teffarr[i]

  ; LOGG LOOP
  ;---------------
  ;FOR j=0,loggpars[2]-1 do BEGIN
  FOR j=0,nlogg-1 do BEGIN

    ;ilogg = loggpars[0] + j*loggpars[1]
    ilogg = loggarr[j]

    ; METAL LOOP
    ;----------------
    ;FOR k=0,metalpars[2]-1 do BEGIN
    FOR k=0,nmetal-1 do BEGIN

      ;imetal = metalpars[0] + k*metalpars[1]
      imetal = metalarr[k]

      ; ALPHA LOOP
      ;---------------
      ;FOR l=0,alphapars[2]-1 do BEGIN
      FOR l=0,nalpha-1 do BEGIN

        ;ialpha = alphapars[0] + l*alphapars[1]
        ialpha = alphaarr[l]

        ; PREPARE THE SPECTRUM

        synpars = double([iteff,ilogg,imetal,ialpha])


        ; STEP 1: Interpolate the spectra
        undefine,wave,spec,error
        SPECFIT_MODELINTERP,synstr0,iteff,ilogg,imetal,ialpha,ww,ss,error=error,head=head,/silent
        if n_elements(error) gt 0 then begin
          print,error
          goto,bomb
        endif

        ; PREP the spectra
        SPECFIT_PREPSPEC,ww,ss,0,0,specpars,wave2,spec2


        ;; STEP 2: Gaussian Smooth
        ;SPECFIT_SETRES,ww,ss,res,ss2
        ;
        ;; STEP 3: Set onto logarithmic dispersion scale
        ;;SPECFIT_SETDISP,ww,ss2,w0,w1,dw,wave2,spec2a,/log
        ;;nw = floor((w1-w0)/dw)
        ;; w1 = 10.^((nw-1)*dwlog+w0log)
        ;;dwlog = (alog10(double(w1))-alog10(double(w0)))/(nw-1.0d0)
        ;;wave2 = (10.0d0)^(dindgen(nw)*dwlog+alog10(double(w0)))
        ;wave2 = (10.0d0)^(dindgen(nw)*dwlog+double(w0log))
        ;spec2a = cspline(ww,ss2,wave2)  ; cspline is ~25x faster
        ;
        ;; STEP 4: Normalize
        ;SPECFIT_NORM,wave2,spec2a,spec2,cont
        ;spec2 = float(spec2)

        ; Update the header
        head2 = head
        SXADDPAR,head2,'TEFF',iteff
        SXADDPAR,head2,'LOGG',ilogg
        SXADDPAR,head2,'FEH',imetal
        SXADDPAR,head2,'ALPHA_FE',ialpha
        SXADDPAR,head2,'CRVAL1',w0log
        SXADDPAR,head2,'CDELT1',dwlog

        SXADDPAR,head2,'PREPARED',1
        SXADDPAR,head2,'RES',res


        ; STEP 5: Save
        name = MKSYNTHNAME(iteff,ilogg,imetal,ialpha)
        newname = outdir+name+'.fits'
        FITS_WRITE,newname,spec2,head2


        BOMB:

      ENDFOR ; alpha loop
    ENDFOR ; metal loop
  ENDFOR ; logg loop
ENDFOR ; teff loop

; Make a SYNSTR for the prepared spectra
MKSYNSTR,outdir,teffarr,loggarr,metalarr,alphaarr,synstr
print,'Saving SYNSTR of the PREPARED spectra to '+outdir+'/synstr.dat'
SAVE,synstr,file=outdir+'/synstr.dat'

;stop

end
