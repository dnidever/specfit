pro specfit_multi,wave,spec,errspec,ssyn,wcenarr,widtharr,fitstr

; This compares an observed and sythetic spectrum in several separate
; regions and computes one grand total chisq value.
;
; INPUTS:
;  wave      Observed wavelength array
;  spec      Observed spectrum
;  errspec   Observed noise/error spectrum
;  ssyn      Synthetic spectrum
;  wcenarr   Array of central wavelength for the regions
;  widtharr  Width for each region (+/-).

; The observed spectrum should ALREADY have the velocity/doppler shift removed
;
; The spectra should ALREADY be on the SAME WAVELENGTH SCALE.

nwave = n_elements(wave)
nregions = n_elements(wcenarr)

; Loop through the regions
deviates = wave*0.0
; Deviates = (Y-Model)/error
count = 0L
FOR i=0,nregions-1 do begin

  ; Get region
  icen = wcenarr[i]
  iwidth = widtharr[i]
  ind = where(abs(wave-icen) le iwidth,nind)
  contwidth = (5*iwidth) > 10.0
  cind = where(abs(wave-icen) lt contwidth and abs(wave-icen) gt iwidth,ncind)

  ispec = spec[ind]
  ierr = errspec[ind]
  isyn = ssyn[ind]
  icspec = spec[cind]  ; continuum
  icsyn = ssyn[cind]

  ;-------------------------
  ; Fit continuum for each
  ;-------------------------
  npoly = 1 ;3
  x = scale(wave[cind],minmax(wave[cind]),[-1.0,1.0])
  y = icspec
  coef = robust_poly_fit(x,y,npoly)
  ;cont = poly(x,coef)*max(spec[ind])
  x2 = scale(wave[ind],minmax(wave[cind]),[-1.0,1.0])
  cont1 = poly(x2,coef)
  inspec = ispec/cont1
  inerr = ierr/cont1

  y = icsyn
  coef2 = robust_poly_fit(x,y,npoly)
  cont2 = poly(x2,coef2)
  insyn = isyn/cont2

  ;plot,wave[ind],inspec,yr=[0.9,1.1],xs=1,ys=1,xtit='Wavelength (A)',tit=strtrim(i+1,2)+' WAVE='+strtrim(icen,2)
  ;oplot,wave[ind],insyn,co=250
  ;oplot,wave[ind],wave[ind]*0.0+1.0,linestyle=2
  plot,wave,spec,xr=[-20,20]+icen,yr=[0.8,1.1],xs=1,ys=1
  oplot,wave,ssyn,co=250

  stop

  ; Calculate deviates
  dev = (inspec-insyn)/inerr
  deviates[count:count+nind-1] = dev

  ; Increment counter
  count = count+nind

END

stop

; Clip the end off deviates
deviates = deviates[0:count-1]

; Final chisq
chisq = total(deviates^2)
n = n_elements(deviates)
rchisq = chisq/n

stop

end
