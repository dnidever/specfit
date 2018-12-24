function specfit_genetic_func,par,_Extra=extra
  
; This is the function that the genetic algorithm calls
  
;COMMON specfit,wave2,spec2,errspec2,mask2,specpars,zerovel
COMMON specfit,synstr,ft,genestr

synpars = reform(par)

;extra = ft
;wave = extra.wave
;spec = extra.spec
;errspec = extra.errspec
;mask = extra.mask
;synstr = extra.synstr
;specpars = extra.specpars
;zerovel = extra.zerovel
;nearest = extra.nearest
;interp = extra.interp

undefine,chisq,vrel,error
SPECFIT_COMPARESYNSPEC,ft.wave,ft.spec,synstr,synpars,ft.specpars,chisq,vrel,vrelerr,vsini,errspec=ft.errspec,$
                       mask=ft.mask,zerovel=ft.zerovel,nearest=ft.nearest,interp=ft.interp,$
                       fitvsini=ft.fitvsini,fixvsini=ft.fixvsini,error=error,/silent
;SPECFIT_COMPARESYNSPEC,ft.wave,ft.spec,ft.synstr,synpars,ft.specpars,chisq,vrel,vrelerr,errspec=ft.errspec,$
;                       mask=ft.mask,zerovel=ft.zerovel,nearest=ft.nearest,interp=ft.interp,error=error,/silent
;SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars,specpars,chisq,vrel,vrelerr,errspec=errspec,mask=mask,$
;                       zerovel=zerovel,nearest=nearest,interp=interp,error=error,/silent

;if n_elements(error) eq 0 then stop

;chisq=99.
;vrel=99.
;vrelerr=99.

if n_elements(error) gt 0 then begin
  large_number = 1d30
  chisq = large_number
  vrel = large_number
  vrelerr = large_number
  vsini = large_number
endif

;print,extra.count,synpars,chisq,vrel

; Keep track of the outputs
;extra.par[extra.count,*] = synpars
;extra.chisq[extra.count] = chisq
;extra.vrel[extra.count] = vrel
;extra.vrelerr[extra.count] = vrelerr
;extra.count++

;genestr.par[genestr.count,*] = synpars
;genestr.chisq[genestr.count] = chisq
;genestr.vrel[genestr.count] = vrel
;genestr.vrelerr[genestr.count] = vrelerr
;genestr.timearr[genestr.count] = systime(1)
;genestr.count++

;stop

return,chisq

end
