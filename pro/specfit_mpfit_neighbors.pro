pro specfit_mpfit_neighbors,par,initchi,chisqarr,steparr,all=all,_Extra=extra
  
; Get Chisq at the neighboring pixels
  
;COMMON specfit_mpfit,mpstr,chigrid
COMMON specfit_mpfit,chigrid

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
fitvsini = extra.fitvsini

npts = n_elements(spec)
npar = n_elements(par)

bestind = specfit_mpfit_par2ind(par,synstr)

; Nearest neighbors
;---------------------
If not keyword_set(all) then begin

  ; Get the model at this point
  ;-----------------------------
  ichigrid = chigrid[bestind[0],bestind[1],bestind[2],bestind[3]]       
  if ichigrid eq 0.0 then begin

    undefine,chisq,vrel,error
    SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars,specpars,initchi,vrel,vrelerr,vsini,errspec=errspec,mask=mask,$
                           zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,error=error,$
                           /silent,wsynth=wsynth,ssynth=ssynth

    if n_elements(error) gt 0 then begin
      initchi = large_number
    endif
    chigrid[bestind[0],bestind[1],bestind[2],bestind[3]] = initchi

  endif else initchi=ichigrid


  ; Get neighbors
  ;----------------
  chisqarr = fltarr(4,2)+large_number
  steparr = fltarr(4,2)

  ngrid = [n_elements(synstr.teff),n_elements(synstr.logg),n_elements(synstr.metal),n_elements(synstr.alpha)]
  tagind = [0,1,2,3]

  ; Loop through the parameters
  for i=0,npar-1 do begin

    ; We can find a partial derivative
    if ngrid[i] gt 1 then begin

      ; PLUS step
      ;------------
      shft = lonarr(4)
      shft[i] = 1
      indout = specfit_mpfit_shift(bestind,synstr,shft,parin,parout)

      ; Not at end
      if (indout[0] ne -1) then begin
        ;synpars_plus = par
        ;synpars_plus[i] = (synstr.(tagind[i]))[bestind[i]+1]
        synpars_plus = parout
        step = abs(synpars_plus[i] - par[i])

        ichigrid = chigrid[indout[0],indout[1],indout[2],indout[3]] 
        if ichigrid eq 0.0 then begin

          undefine,chisq,vrel,error
          SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars_plus,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec,mask=mask,$
                                 zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,error=error,$
                                 /silent,wsynth=wsynth_plus,ssynth=ssynth_plus

          if n_elements(error) gt 0 then begin
            chisq = large_number
          endif
          chigrid[indout[0],indout[1],indout[2],indout[3]] = chisq

        endif else chisq=ichigrid

      ; Already at the end
      endif else begin
        chisq = large_number
        step = 0.0
      endelse
      chisqarr[i,1] = chisq
      steparr[i,1] = step

      ; MINUS step
      ;-------------
      shft = lonarr(4)
      shft[i] = -1
      indout = specfit_mpfit_shift(bestind,synstr,shft,parin,parout)

      ; Not at beginning
      if (bestind[i] gt 0) then begin
        ;synpars_minus = par
        ;synpars_minus[i] = (synstr.(tagind[i]))[bestind[i]-1]
        synpars_minus = parout
        step = abs(synpars_minus[i] - par[i])

        ichigrid = chigrid[indout[0],indout[1],indout[2],indout[3]] 
        if ichigrid eq 0.0 then begin

          undefine,chisq,vrel,error
          SPECFIT_COMPARESYNSPEC,wave,spec,synstr,synpars_minus,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec,mask=mask,$
                                 zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,error=error,$
                                 /silent,wsynth=wsynth_minus,ssynth=ssynth_minus

          if n_elements(error) gt 0 then begin
            chisq = large_number
          endif
          chigrid[indout[0],indout[1],indout[2],indout[3]] = chisq

        endif else chisq=ichigrid

      ; Already at the beginning
      endif else begin
        chisq = large_number
        step = 0.0
      endelse
      chisqarr[i,0] = chisq
      steparr[i,0] = step

    endif

  end  ; parameter loop

; ALL neighbors
;----------------
Endif Else Begin

  chisqarr = fltarr(3,3,3,3)+large_number

  ; Teff loop
  for i=-1,1 do begin

    ; logg loop
    for j=-1,1 do begin

      ; Metal loop
      for k=-1,1 do begin

        ; Alpha loop
        for l=-1,1 do begin

          shft = [i,j,k,l]
          indout = specfit_mpfit_shift(bestind,synstr,shft,parin,parout)

          ; Okay
          if (indout[0] ne -1) then begin

            ichigrid = chigrid[indout[0],indout[1],indout[2],indout[3]] 
            if ichigrid eq 0.0 then begin

              undefine,chisq,vrel,error
              SPECFIT_COMPARESYNSPEC,wave,spec,synstr,parout,specpars,chisq,vrel,vrelerr,vsini,errspec=errspec,mask=mask,$
                                     zerovel=zerovel,nearest=nearest,interp=interp,fitvsini=fitvsini,error=error,$
                                     /silent,wsynth=wsynth_plus,ssynth=ssynth_plus

              if n_elements(error) gt 0 then begin
                chisq = large_number
              endif
              chigrid[indout[0],indout[1],indout[2],indout[3]] = chisq

            endif else chisq=ichigrid

          ; Already at the end
          endif else begin
            chisq = large_number
          endelse
          chisqarr[i+1,j+1,k+1,l+1] = chisq


        end ; alpha loop
      end ; metal loop
    end ; logg loop
  end ; Teff loop


Endelse ; ALL neighbors

;stop

end
