function specfit_mpfit_func,x,par,_Extra=extra
  
; This is the function that MPFIT calls
  
COMMON specfit_mpfit,mpstr

synpars = reform(par)

wave = extra.wave
spec = extra.spec
errspec = extra.errspec
mask = extra.mask
synstr = extra.synstr
specpars = extra.specpars
zerovel = extra.zerovel
nearest = extra.nearest
fitvsini = extra.fitvsini
interp = extra.interp

undefine,chisq,vrel,error
SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec,mask=mask,$
                       zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,error=error,$
                       /silent,wsynth=wsynth,ssynth=ssynth

if n_elements(error) gt 0 then begin

  large_number = 1d30
  chisq = large_number
  vrel = large_number
  vrelerr = large_number

  ssynth = wave*0.0+large_number

endif

; Keep track of the outputs
mpstr.par[mpstr.count,*] = synpars
mpstr.chisq[mpstr.count] = chisq
mpstr.vrel[mpstr.count] = vrel
mpstr.vrelerr[mpstr.count] = vrelerr
mpstr.count++

;stop

return,ssynth

end
