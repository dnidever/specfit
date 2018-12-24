function specfit_mpfit_dev,par,dp,noderiv=noderiv,_Extra=extra
  
; This is the function that MPFIT calls and returns the deviates
  
COMMON specfit_mpfit,mpstr

large_number = 1d30

synpars = reform(par)

wave = extra.wave
spec = extra.spec
errspec = extra.errspec
mask = extra.mask
synstr = extra.synstr
specpars = extra.specpars
zerovel = extra.zerovel
nearest = extra.nearest
interp = extra.interp

npts = n_elements(spec)
npar = n_elements(par)

; Get the model at this point
;-----------------------------
undefine,chisq,vrel,error
SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,$
                       zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent,wsynth=wsynth,ssynth=ssynth

if n_elements(error) gt 0 then begin
  dp = fltarr(npts,npar)+large_number
  return,spec*0.0+large_number
endif

resid = spec-ssynth
dev = resid/errspec


; Calculate the partial derivatives
; dp = [npts,npar]
if not keyword_set(noderiv) and n_params() gt 1 then begin

  dp = fltarr(npts,npar)

  ; Find the closest grid point
  teffind = first_el(minloc(abs(synstr.teff-par[0])))
  loggind = first_el(minloc(abs(synstr.logg-par[1])))
  metalind = first_el(minloc(abs(synstr.metal-par[2])))
  alphaind = first_el(minloc(abs(synstr.alpha-par[3])))

  bestind = [teffind, loggind, metalind, alphaind]
  ngrid = [n_elements(synstr.teff),n_elements(synstr.logg),n_elements(synstr.metal),n_elements(synstr.alpha)]
  tagind = [0,1,2,3]

  ; Loop through the parameters
  for i=0,npar-1 do begin

    ; We can find a partial derivative
    if ngrid[i] gt 1 then begin

      ; PLUS step
      ;------------
      ; Not at end
      if (bestind[i] lt (ngrid[i]-1)) then begin
        synpars_plus = par
        synpars_plus[i] = (synstr.(tagind[i]))[bestind[i]+1]

        undefine,chisq,vrel,error
        SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars_plus,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,$
                               zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent,wsynth=wsynth_plus,$
                               ssynth=ssynth_plus

        if n_elements(error) gt 0 then begin
          dp = fltarr(npts,npar)+large_number
          return,dev
        endif

        par_plus = synpars_plus[i]
        resid_plus = spec-ssynth_plus

      ; Already at the end
      endif else begin
        par_plus = par[i]
        resid_plus = resid
      endelse

      ; MINUS step
      ;-------------
      ; Not at beginning
      if (bestind[i] gt 0) then begin
        synpars_minus = par
        synpars_minus[i] = (synstr.(tagind[i]))[bestind[i]-1]

        undefine,chisq,vrel,error
        SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars_minus,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,$
                               zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent,wsynth=wsynth_minus,$
                               ssynth=ssynth_minus

        if n_elements(error) gt 0 then begin
          dp = fltarr(npts,npar)+large_number
          return,dev
        endif

        par_minus = synpars_minus[i]
        resid_minus = spec-ssynth_minus

      ; Already at the beginning
      endif else begin
        par_plus = par[i]
        resid_plus = resid
      endelse

      ; Calculate the partial derivative
      step = par_plus - par_minus
      ;step = 1.0
      dp[*,i] = ((ssynth_plus-ssynth_minus)/step) / errspec

      ;stop

    endif ; can find derivative

  end  ; parameter loop

  ;dp = -dp

endif  ; finding deriative


;if n_elements(error) gt 0 then begin
;
;  large_number = 1d30
;  chisq = large_number
;  vrel = large_number
;  vrelerr = large_number
;
;  ssynth = wave*0.0+large_number
;
;endif
;
; Keep track of the outputs
;mpstr.par[mpstr.count,*] = synpars
;mpstr.chisq[mpstr.count] = chisq
;mpstr.vrel[mpstr.count] = vrel
;mpstr.vrelerr[mpstr.count] = vrelerr
;mpstr.count++

;stop

return,dev

end
