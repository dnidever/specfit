pro specfit_norm,wave0,spec0,errspec0,spec2,cont,plot=pl,stp=stp,error=error
;pro specfit_norm,wave0,spec0,spec2,cont,plot=pl,stp=stp,error=error

;+
;
; SPECFIT_NORM
;
; This program normalizes a spectrum.  It uses the error spectrum
; to identify outliers (i.e. lines).  The "good" points are binned
; and then a cubic spline is used to obtain the continuum for the
; entire spectrum.
;
; INPUTS:
;  wave      The wavelength array
;  spec      The spectrum array
;  errspec   The error spectrum.   Set this to 0 if none exists (i.e. a
;               synthetic spectrum).
;  /plot     Plot spectra.
;  /stp      Stop at the end of the program
;
; OUTPUTS:
;  spec2     Normalized spectrum
;  cont      The continuum
;  =error    The error message if there was one.
;
; USAGE:
;  IDL>specfit_norm,wave,spec,errspec,spec2,cont
;
; By D.Nidever  Oct 2008
;-

undefine,error,spec2

nwave = n_elements(wave0)
nspec = n_elements(spec0)
nerrspec = n_elements(errspec0)

; Not enough inputs
if nwave eq 0 or nspec eq 0 or nerrspec eq 0 then begin
  print,'Syntax - specfit_norm,wave,spec,errspec,spec2,cont,stp=stp,error=error'
  error = 'Not enough inputs'
  return
endif

; Wave and spec don't match
if nwave ne nspec then begin
  error = 'WAVE and SPEC must have same number of elements'
  print,error
  return
endif

wave = wave0
spec = spec0

; Error spectrum
if n_elements(errspec0) gt 0 then errspec=errspec0 else errspec=dblarr(nspec)
if n_elements(errspec) eq 1 and errspec[0] eq 0.0 then errspec=dblarr(nspec)  ; synthetic spectrum
if n_elements(errspec) ne nwave then begin
  error = 'ERRSPEC and WAVE/SPEC must have same number of elements (or 0 for synthetic spectrum)'
  print,error
  return
endif

;goto,skiptohere

; I THINK THE SPECTRUM SHOULD BE BINNED TO ~1A PER PIXEL
; BEFORE SMOOTHING, AND THEN SPLINE TO THE FULL SPECTRUM
; OTHERWISE IT COULD BE *VERY* SLOW FOR HIGH-RES.

; Dispersion
dw = median(slope(wave))


;##########################
; NINTH ATTEMPT
;##########################

; Just use wide Gaussian smoothing
smlen = round(200./dw) < nwave
;smlen = round(150./dw) < nwave
;smlen = round(100./dw) < nwave

; If hig-resolution, bin first.
dwthresh = 5.0  ;1.0
if dw lt dwthresh then begin

  ; Bin the spectrum
  binsize = round(dwthresh/dw) > round(dw*2)
  BINDATA,wave,spec,wbin,sbin,binsize=binsize
  nbin = n_elements(wbin)

  ; Just use wide Gaussian smoothing
  dw = wbin[1]-wbin[0]
  smlen = round(200./dw) < nbin

  ;; add the ends
  ;wbin = [wave[0],wbin,wave[nwave-1]]
  ;sbin = [spec[0],sbin,spec[nwave-1]]

  next = 2*smlen < nbin ;200  ; 100
  npoly = next < nbin ;400  ;200
  ;smspec1 = gsmooth(spec,5)
  smspec1 = savgolsm(sbin,[5,5,4])

  ;coef1 = linfit(findgen(npoly),smspec1[0:npoly-1])
  ;coef2 = linfit(findgen(npoly),smspec1[nwave-npoly:nwave-1])
  coef1 = robust_linefit(findgen(npoly),smspec1[0:npoly-1])
  coef2 = robust_linefit(findgen(npoly),smspec1[nbin-npoly:nbin-1])
  ;coef1 = robust_poly_fit(findgen(npoly),smspec1[0:npoly-1],1)
  ;coef2 = robust_poly_fit(findgen(npoly),smspec1[nwave-npoly:nwave-1],1)
  f1 = poly(findgen(next)-next,coef1)
  f2 = poly(findgen(next)+npoly,coef2)
  w1 = findgen(next)*dw-next*dw+min(wbin)
  w2 = findgen(next)*dw+dw+max(wbin)

  ;plot,findgen(nwave),spec,xr=[-next,nwave+next],xs=1
  ;oplot,findgen(next)-next,f1,co=250
  ;oplot,findgen(next)+nwave,f2,co=250

  extspec = [f1,sbin,f2]
  extwave = [w1,wbin,w2]

  ; GAUSSIAN
  extsmspec = gsmooth(extspec,smlen,widfwhm=3.0)

  ; Spline to the original wavelength scale
  cont = CSPLINE(extwave,extsmspec,wave)

  spec2 = spec/cont


; Low resolution, no binning
endif else begin

  next = 2*smlen < nwave ;200  ; 100
  npoly = next < nwave ;400  ;200
  ;smspec1 = gsmooth(spec,5)
  smspec1 = savgolsm(spec,[5,5,4])

  ;coef1 = linfit(findgen(npoly),smspec1[0:npoly-1])
  ;coef2 = linfit(findgen(npoly),smspec1[nwave-npoly:nwave-1])
  coef1 = robust_linefit(findgen(npoly),smspec1[0:npoly-1])
  coef2 = robust_linefit(findgen(npoly),smspec1[nwave-npoly:nwave-1])
  ;coef1 = robust_poly_fit(findgen(npoly),smspec1[0:npoly-1],1)
  ;coef2 = robust_poly_fit(findgen(npoly),smspec1[nwave-npoly:nwave-1],1)
  f1 = poly(findgen(next)-next,coef1)
  f2 = poly(findgen(next)+npoly,coef2)

  ;plot,findgen(nwave),spec,xr=[-next,nwave+next],xs=1
  ;oplot,findgen(next)-next,f1,co=250
  ;oplot,findgen(next)+nwave,f2,co=250

  extspec = [f1,spec,f2]

  ; GAUSSIAN
  ;extsmspec = gsmooth(extspec,smlen)
  extsmspec = gsmooth(extspec,smlen,widfwhm=3.0)

  ; Remove extended ends
  gsmspec = extsmspec[next:n_elements(extsmspec)-next-1]

  ;gsmspec = gsmooth(spec,smlen)

  cont = gsmspec
  spec2 = spec/cont

endelse

if keyword_set(stp) then stop

return



;##########################
; EIGHTH ATTEMPT
;##########################
;
;; Use Gaussian smoothed spec
;; spline pieces
;
;smlen = round(200./dw) < nwave
;
;next = 2*smlen < nwave ;200  ; 100
;npoly = next < nwave ;400  ;200
;;smspec1 = gsmooth(spec,5)
;smspec1 = savgolsm(spec,[5,5,4])
;coef1 = robust_poly_fit(findgen(npoly),smspec1[0:npoly-1],1)
;coef2 = robust_poly_fit(findgen(npoly),smspec1[nwave-npoly:nwave-1],1)
;f1 = poly(findgen(next)-next,coef1)
;f2 = poly(findgen(next)+npoly,coef2)
;
;;plot,findgen(nwave),spec,xr=[-next,nwave+next],xs=1
;;oplot,findgen(next)-next,f1,co=250
;;oplot,findgen(next)+nwave,f2,co=250
;
;extspec = [f1,spec,f2]
;
;
;; Smooth the spectrum
;;--------------------
;; SAVGOL
;;smspec = savgolsm(spec,[30,30,4])
;;smspec = savgolsm(spec,[25,25,4])  ; good
;;smspec = savgolsm(spec,[25,25,3])
;;smspec = savgolsm(spec,[25,25,2])
;;smspec = savgolsm(spec,[15,15,4])
;;smspec = savgolsm(spec,[10,10,4])
;;smspec = savgolsm(spec,[5,5,4])
;; GAUSSIAN
;extsmspec = gsmooth(extspec,smlen)
;;smspec = gsmooth(extspec,50)   ; very good, except at edges
;;smspec = gsmooth(spec,40)
;;smspec = gsmooth(spec,30)
;;smspec = gsmooth(spec,20)
;;smspec = gsmooth(spec,10)
;;smspec = gsmooth(spec,5)
;; BOXCAR
;;smspec = smooth(spec,smlen,/edge_truncate)
;;smspec = smooth(spec,60,/edge_truncate)  ; problems at edges
;;smspec = smooth(spec,50,/edge_truncate)
;;smspec = smooth(spec,40,/edge_truncate)
;;smspec = smooth(spec,30,/edge_truncate)  ;good
;;smspec = smooth(spec,20,/edge_truncate)
;;smspec = smooth(spec,10,/edge_truncate)
;;smspec = smooth(spec,5,/edge_truncate)
;; MEDIAN
;;smspec = median(spec,50)
;;smspec = median(spec,40)
;;smspec = median(spec,30)
;;smspec = median(spec,20)
;;smspec = median(spec,10)
;;smspec = median(spec,5)
;; HANNING
;;han = hanning(20)
;;han = hanning(10)
;;han = hanning(5)
;;smspec = convol(spec,han,/edge_truncate,/norm)
;; LEE
;;smspec = leefilt(spec,15)
;;smspec = leefilt(spec,10)
;;smspec = leefilt(spec,5)
;; low pass filter
;;Coeff = DIGITAL_FILTER(0, 0.2, 50, nwave)
;;smspec = CONVOL(spec, Coeff,/edge_truncate) 
;
;; Remove extended ends
;gsmspec = extsmspec[next:n_elements(extsmspec)-next-1]
;
;; Smoothing reduces error.  Add gaussian-weighted errors in quadrature
;;errspec2 = sqrt( gsmooth(errspec^2.0,smlen) )
;;errspec = errspec2
;
;; First order continuum is wide gaussian smoothing
;cont1 = gsmspec
;spec1 = spec/cont1
;errspec1 = errspec/cont1
;;cont = gsmspec
;;spec2 = spec/cont
;;return
;
;;stop
;
;; plotting
;;pl = 1
;;midpl = 1
;;growpl = 1
;
;; Dispersion
;dw = median(slope(wave))
;
;nspec = n_elements(spec)
;x = findgen(nspec)
;;tspec = spec
;;tspec = smspec
;;tspec = spec1
;;tspec = savgolsm(spec1,[30,30,4])
;;tspec = gsmooth(spec1,30)
;tspec = smooth(spec1,30,/edge)
;;tspec = smooth(spec1,10,/edge)
;;tspec = smooth(spec1,20,/edge)
;tx = x
;smspec = tspec
;;smspec = smooth(spec1,30,/edge)
;smspec = savgolsm(tspec,[10,10,4])
;;smspec = gsmooth(tspec,10)
;;smspec = tspec
;
;; Initial continuum 
;;nbin = 10 >  round(200./dw)  ; 200A width
;nbin = 10 >  round(100./dw)  ; 200A width
;;nbin = 10 >  round(50./dw)  ; 200A width
;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;tcont = CSPLINE(xmnbin,specbin,x)
;
;; Outlier rejection using the error spectrum
;niterate = 10  ;20
;;order = 15
;
;count = 0
;flag = 0
;WHILE (flag eq 0) do begin
;
;  ;ratio = tspec/tcont
;  ratio = spec1/tcont
;  smratio = smspec/tcont
;  ;smslope = slope(smratio,/acc)
;
;  ; Find outliers using error spectrum
;  ;------------------------------------
;
;  ; Parameters
;  order =  30 ;15 ;20 ;10 ; 20.
;  low_rej = 3 ;4 ;3   ; 0.9 ;2
;  high_rej = 0 ;10 ;5  ; 0 ;1.1 ;0
;  ; Reject positive outliers once a decent continuum has been found
;  if count ge 5 then high_rej = 6 ;3
;  grow = 1 ;0 ;1 ;0
;
;  ;It doesn't make sense to use the error spectrum to find outliers if we're going to
;  ;smooth the spectrum a lot.
;  ;stop
;
;  ; Find outliers
;  sig = errspec1/tcont < 0.2      ; 0.2 maximum
;  sig = 0.01 > errspec1/tcont     ; 0.01 for synthetic spectrum
;  high_thresh = 1.0+high_rej*sig
;  low_thresh = 1.0-low_rej*sig
;  ;high_thresh = 1.2+fltarr(nspec)
;  ;low_thresh = 0.8+fltarr(nspec)
;  nbdhi = 0
;  if high_rej gt 0 then bdhi = where(ratio gt high_thresh,nbdhi)
;  if nbdhi gt 0 then tspec[bdhi] = !values.f_nan
;  nbdlo = 0
;  if low_rej gt 0 then bdlo = where(ratio lt low_thresh,nbdlo)
;  if nbdlo gt 0 then tspec[bdlo] = !values.f_nan
;
;  ; DON'T GROW!!  This causes problems by throwing out few remaining
;  ; good points b/w close, wide absorption lines.
;
;  ; Grow
;  if grow gt 0 then begin
;    badmask = x*0.
;    ;if nbdhi gt 0 then bad[bdhi]=1
;    if nbdlo gt 0 then badmask[bdlo]=1
;    gkernel = fltarr(long(grow)+2)+1.0
;    flag2 = 0
;    count2 = 0
;
;    ; Use shift instead of convolve
;    ; make a lftgtr and rtgrt arrays first before the loop 
;    ;   to check that points are getting larger
;    smratio_lftnbor = [ 0.0, smratio[0:nwave-2] ]  ; left neighbor
;    smratio_rtnbor = [ smratio[1:nwave-1], 0.0 ]   ; right neighbor
;    lftgtr = float(smratio gt smratio_lftnbor)   ; greater than left neighbor
;    rtgtr = float(smratio gt smratio_rtnbor)     ; greater than right neighbor
;
;    ;growpl = 1
;    if keyword_set(growpl) then begin
;      plot,smratio,ps=-1
;      b = where(badmask eq 1.0)
;      oplot,b,smratio[b],ps=4,co=250
;    endif
;
;    badmask0 = badmask
;    tbadmask = badmask
;    while (flag2 eq 0) do begin
;
;      badmask_lft = [ 0.0, badmask[0:nwave-2] ]  ; badmask to the left 1 element
;      badmask_rt = [ badmask[1:nwave-1], 0.0 ]   ; badmask to the right 1 element
;
;      ht_thresh = 0.98
;      pk_thresh = 0.90 ;0.95
;      newbadmask = float( badmask eq 0.0 and smratio lt ht_thresh and ratio lt ht_thresh and $
;                          ((badmask_lft eq 1.0 and lftgtr eq 1.0) OR $
;                          (badmask_rt eq 1.0 and rtgtr eq 1.0)) and $
;                          ((smratio lt pk_thresh and ratio lt pk_thresh) OR lftgtr eq 0.0 OR rtgtr eq 0.0) )
;      nnewbadmask = total(newbadmask)
;
;      ; This method still has problem in the blue part of noisy spectra
;      ; where the lines are so close and the isn't much continuum b/w them.
;
;      if nnewbadmask gt 0 then begin
;        badmask = badmask+newbadmask
;        badmask = badmask/(badmask>1.0)
;      endif
;      ;if nkeep gt 0 then badmask[keep] = 1.0
;
;      ; Finished
;      if nnewbadmask eq 0 then flag2 =1
;
;      if keyword_set(growpl) then begin
;        b = where(badmask eq 1.0)
;        oplot,b,smratio[b],ps=1,co=150
;
;        ;stop
;        wait,0.1
;      endif
;
;      ;stop
;
;      count2++
;    end
;
;    ;stop
;    bdgrow = where(badmask eq 1.0,nbdgrow)
;    if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;
;  endif ; grow
;
;  if keyword_set(growpl) then stop
;
;  ;if grow gt 0 then begin
;  ;  bad = x*0.
;  ;  if nbdhi gt 0 then bad[bdhi]=1
;  ;  if nbdlo gt 0 then bad[bdlo]=1
;  ;  gkernel = fltarr(long(grow)+2)+1.0
;  ;  bad2 = convol(bad,gkernel,/center,/edge_truncate)
;  ;  bdgrow = where(bad2 gt 0.5,nbdgrow)
;  ;  if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;  ;  if nbdgrow gt 0 then tx[bdgrow] = !values.f_nan
;  ;endif
;
;  ;bd = where(finite(tspec) eq 0,nbd)
;  ;if nbd gt 0 then tx[bd] = !values.f_nan
;
;  ; Get the continuum
;  ;-------------------
;
;  ; Bin the data
;  ;nbin = ceil( float(nspec)/float(order) )
;  ;nbin = 10 >  round(100./dw)  ; 100A width
;  ;nbin = 10 >  round(100./dw)  ; 100A width
;  BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;  ;perc = 0.5 ;1.0-median(errspec1)
;  ;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0,perc=perc
;
;  ; Spline
;  gd = where(finite(xmnbin) eq 1 and finite(specbin) eq 1,ngd)
;  if ngd eq 0 then begin
;    error = 'No more good points for splining the continuum'
;    print,error
;    stop
;  endif
;  ;spl_cont = CSPLINE(xbin,specbin,x)
;  spl_cont = CSPLINE(xmnbin[gd],specbin[gd],x)
;
;  tcont_last = tcont
;  tcont = spl_cont
;
;  ;; Normalized spectrum
;  ;ratio = tspec/spl_cont
;
;
;  ; Finished
;  if count ge niterate or max(abs(tcont_last-tcont)) lt 1e-6 then flag=1
;
;  ; Plot continuum normalization
;  ;midpl = 1
;  if keyword_set(midpl) then begin
;    !p.multi=[0,1,2]
;    plot,spec1,tit='Iteration='+strtrim(count,2)
;    oplot,tspec,co=200
;    gd = where(finite(tspec) eq 1,ngd)
;    oplot,gd,spec1[gd],ps=1,co=150
;    oplot,spl_cont,co=250,thick=2.0
;    ;oplot,xbin,specbin,ps=1,co=250,sym=2
;    oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;  
;    ;plot,spec/tcont,yr=[0,2.0],ys=1
;    plot,ratio,yr=[0,2.0],ys=1
;    oplot,gd,ratio[gd],ps=1,co=150
;    oplot,minmax(x),[1,1],co=200
;    if high_rej gt 0 then oplot,high_thresh,co=250
;    if low_rej gt 0 then oplot,low_thresh,co=250
;    !p.multi=0  
;
;    wait,0.1
;    stop
;  end
;
;  ;stop
;
;  count++
;
;ENDWHILE
;
;
;;pl = 1
;; plotting
;if keyword_set(pl) then begin
;  !p.multi=[0,1,2]
;  plot,spec1,xs=1
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec1[gd],ps=1,co=150
;  oplot,spl_cont,co=250,thick=2.0
;  ;oplot,xbin,specbin,ps=1,co=250,sym=2
;  oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;
;  ;!p.multi=[1,1,2]
;  plot,[0],/nodata,xr=minmax(x),yr=[0.0,1.5],xs=1,ys=1
;  oplot,spec1/spl_cont
;  oplot,gd,(spec1/spl_cont)(gd),ps=1,co=150
;  ;oplot,spec/cont4,co=200
;  ;oplot,spec/cont2
;  ;oplot,spec/cont1,co=200
;  oplot,minmax(x),[1,1],co=200
;  if high_rej gt 0 then oplot,high_thresh,co=250
;  if low_rej gt 0 then oplot,low_thresh,co=250
;
;  wait,0.2
;  !p.multi=0
;  stop
;endif
;
;;cont = spl_cont
;;spec2 = spec/cont
;
;cont = cont1*spl_cont
;spec2 = spec/cont
;
;;stop
;return
;
;
;
;##########################
; SEVENTH ATTEMPT
;##########################

; Same as FIFTH attempt but growing the deep outlier lines to the
; continuum, to remove the wings.
;
;; plotting
;;pl = 1
;;midpl = 1
;;growpl = 1
;
;; Dispersion
;dw = median(slope(wave))
;
;nspec = n_elements(spec)
;x = findgen(nspec)
;;tspec = spec
;tspec = smspec
;tx = x
;;smspec = savgolsm(spec,[10,10,4])
;;smspec = gsmooth(spec,10)
;;smspec = spec
;
;; Initial continuum 
;nbin = 10 >  round(200./dw)  ; 200A width
;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;tcont = CSPLINE(xmnbin,specbin,x)
;
;; Outlier rejection using the error spectrum
;niterate = 10  ;20
;order = 15
;
;count = 0
;flag = 0
;WHILE (flag eq 0) do begin
;
;  ;ratio = tspec/tcont
;  ratio = spec/tcont
;  smratio = smspec/tcont
;  ;smslope = slope(smratio,/acc)
;
;  ; Find outliers using error spectrum
;  ;------------------------------------
;
;  ; Parameters
;  order =  15 ;20 ;10 ; 20.
;  low_rej = 3 ;4 ;3   ; 0.9 ;2
;  high_rej = 0 ;10 ;5  ; 0 ;1.1 ;0
;  ; Reject positive outliers once a decent continuum has been found
;  if count ge 5 then high_rej = 6 ;3
;  grow = 0 ;1 ;0
;
;
;  ;It doesn't make sense to use the error spectrum to find outliers if we're going to
;  ;smooth the spectrum a lot.
;  ;stop
;
;  ; Find outliers
;  sig = errspec/tcont < 0.2      ; 0.2 maximum
;  sig = 0.01 > errspec/tcont     ; 0.01 for synthetic spectrum
;  high_thresh = 1.0+high_rej*sig
;  low_thresh = 1.0-low_rej*sig
;  ;high_thresh = 1.2+fltarr(nspec)
;  ;low_thresh = 0.8+fltarr(nspec)
;  nbdhi = 0
;  if high_rej gt 0 then bdhi = where(ratio gt high_thresh,nbdhi)
;  if nbdhi gt 0 then tspec[bdhi] = !values.f_nan
;  nbdlo = 0
;  if low_rej gt 0 then bdlo = where(ratio lt low_thresh,nbdlo)
;  if nbdlo gt 0 then tspec[bdlo] = !values.f_nan
;
;  ; DON'T GROW!!  This causes problems by throwing out few remaining
;  ; good points b/w close, wide absorption lines.
;
;  ; Grow
;  if grow gt 0 then begin
;    badmask = x*0.
;    ;if nbdhi gt 0 then bad[bdhi]=1
;    if nbdlo gt 0 then badmask[bdlo]=1
;    gkernel = fltarr(long(grow)+2)+1.0
;    flag2 = 0
;    count2 = 0
;
;    ; Use shift instead of convolve
;    ; make a lftgtr and rtgrt arrays first before the loop 
;    ;   to check that points are getting larger
;    smratio_lftnbor = [ 0.0, smratio[0:nwave-2] ]  ; left neighbor
;    smratio_rtnbor = [ smratio[1:nwave-1], 0.0 ]   ; right neighbor
;    lftgtr = float(smratio gt smratio_lftnbor)   ; greater than left neighbor
;    rtgtr = float(smratio gt smratio_rtnbor)     ; greater than right neighbor
;
;    ;growpl = 1
;    if keyword_set(growpl) then begin
;      plot,smratio,ps=-1
;      b = where(badmask eq 1.0)
;      oplot,b,smratio[b],ps=4,co=250
;    endif
;
;    badmask0 = badmask
;    tbadmask = badmask
;    while (flag2 eq 0) do begin
;
;      badmask_lft = [ 0.0, badmask[0:nwave-2] ]  ; badmask to the left 1 element
;      badmask_rt = [ badmask[1:nwave-1], 0.0 ]   ; badmask to the right 1 element
;
;      ht_thresh = 0.98
;      pk_thresh = 0.90 ;0.95
;      newbadmask = float( badmask eq 0.0 and smratio lt ht_thresh and ratio lt ht_thresh and $
;                          ((badmask_lft eq 1.0 and lftgtr eq 1.0) OR $
;                          (badmask_rt eq 1.0 and rtgtr eq 1.0)) and $
;                          ((smratio lt pk_thresh and ratio lt pk_thresh) OR lftgtr eq 0.0 OR rtgtr eq 0.0) )
;      nnewbadmask = total(newbadmask);
;
;      ; This method still has problem in the blue part of noisy spectra
;      ; where the lines are so close and the isn't much continuum b/w them.
;
;      if nnewbadmask gt 0 then begin
;        badmask = badmask+newbadmask
;        badmask = badmask/(badmask>1.0)
;      endif
;      ;if nkeep gt 0 then badmask[keep] = 1.0
;
;      ; Finished
;      if nnewbadmask eq 0 then flag2 =1
;
;      if keyword_set(growpl) then begin
;        b = where(badmask eq 1.0)
;        oplot,b,smratio[b],ps=1,co=150
;
;        ;stop
;        wait,0.1
;      endif
;
;      ;stop
;
;      count2++
;    end
;
;    ;stop
;    bdgrow = where(badmask eq 1.0,nbdgrow)
;    if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;
;  endif ; grow
;
;  if keyword_set(growpl) then stop;
;
;  ;if grow gt 0 then begin
;  ;  bad = x*0.
;  ;  if nbdhi gt 0 then bad[bdhi]=1
;  ;  if nbdlo gt 0 then bad[bdlo]=1
;  ;  gkernel = fltarr(long(grow)+2)+1.0
;  ;  bad2 = convol(bad,gkernel,/center,/edge_truncate)
;  ;  bdgrow = where(bad2 gt 0.5,nbdgrow)
;  ;  if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;  ;  if nbdgrow gt 0 then tx[bdgrow] = !values.f_nan
;  ;endif;
;
;  ;bd = where(finite(tspec) eq 0,nbd)
;  ;if nbd gt 0 then tx[bd] = !values.f_nan
;
;  ; Get the continuum
;  ;-------------------;
;
;  ; Bin the data
;  ;nbin = ceil( float(nspec)/float(order) )
;  ;nbin = 10 >  round(100./dw)  ; 100A width
;  nbin = 10 >  round(100./dw)  ; 100A width
;  BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;
;  ; Spline
;  gd = where(finite(xmnbin) eq 1 and finite(specbin) eq 1,ngd)
;  if ngd eq 0 then begin
;    error = 'No more good points for splining the continuum'
;    print,error
;    stop
;  endif
;  ;spl_cont = CSPLINE(xbin,specbin,x)
;  spl_cont = CSPLINE(xmnbin[gd],specbin[gd],x)
;
;  tcont_last = tcont
;  tcont = spl_cont
;
;  ;; Normalized spectrum
;  ;ratio = tspec/spl_cont
;
;
;  ; Finished
;  if count ge niterate or max(abs(tcont_last-tcont)) lt 1e-6 then flag=1  
;
;  ; Plot continuum normalization
;  ;midpl = 1
;  if keyword_set(midpl) then begin
;    !p.multi=[0,1,2]
;    plot,spec,tit='Iteration='+strtrim(count,2)
;    gd = where(finite(tspec) eq 1,ngd)
;    oplot,gd,spec[gd],ps=1,co=150
;    oplot,spl_cont,co=250,thick=2.0
;    ;oplot,xbin,specbin,ps=1,co=250,sym=2
;    oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;  
;    ;plot,spec/tcont,yr=[0,2.0],ys=1
;    plot,ratio,yr=[0,2.0],ys=1
;    oplot,gd,ratio[gd],ps=1,co=150
;    oplot,minmax(x),[1,1],co=200
;    if high_rej gt 0 then oplot,high_thresh,co=250
;    if low_rej gt 0 then oplot,low_thresh,co=250
;    !p.multi=0  
;
;    wait,0.1
;    stop
;  end
;
;  ;stop
;
;  count++
;
;ENDWHILE
;
;
;;pl = 1
;; plotting
;if keyword_set(pl) then begin
;  !p.multi=[0,1,2]
;  plot,spec,xs=1
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec[gd],ps=1,co=150
;  oplot,spl_cont,co=250,thick=2.0
;  ;oplot,xbin,specbin,ps=1,co=250,sym=2
;  oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;
;  ;!p.multi=[1,1,2]
;  plot,[0],/nodata,xr=minmax(x),yr=[0.0,1.5],xs=1,ys=1
;  oplot,spec/spl_cont
;  oplot,gd,(spec/spl_cont)(gd),ps=1,co=150
;  ;oplot,spec/cont4,co=200
;  ;oplot,spec/cont2
;  ;oplot,spec/cont1,co=200
;  oplot,minmax(x),[1,1],co=200
;  if high_rej gt 0 then oplot,high_thresh,co=250
;  if low_rej gt 0 then oplot,low_thresh,co=250
;
;  wait,0.2
;  !p.multi=0
;endif
;
;cont = spl_cont
;spec2 = spec/cont
;
;;stop
;return
;
;
;##########################
; SIXTH ATTEMPT
;##########################

; Same as FIFTH attempt but with line fitting and removal

;; Make error spectrum
;if n_elements(errspec0) gt 0 then errspec=errspec0 else errspec=dblarr(nspec)
;if n_elements(errspec) eq 1 then errspec=dblarr(nspec)
;
;wave = wave0
;spec = spec0
;
;x1 = findgen(nspec)
;x2 = scale_vector(findgen(nspec*10.0),min(x1),max(x1))
;wave = cspline(x1,wave0,x2)
;spec = cspline(x1,spec0,x2)
;errspec = cspline(x1,errspec,x2)
;
;
;; Dispersion
;dw = median(slope(wave))
;
;nspec = n_elements(spec)
;x = findgen(nspec)
;tspec = spec
;tx = x
;
;
; ; Initial continuum 
;nbin = 10 >  round(200./dw)  ; 200A width
;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;tcont = CSPLINE(xmnbin,specbin,x)
;
;; Outlier rejection using the error spectrum
;niterate = 10  ;20
;;low_rej = 3   ; 0.9 ;2
;;high_rej = 0 ;5  ; 0 ;1.1 ;0
;;order = 20 ;20 ;10 ; 20.
;;grow = ceil(5.0/dw) ; 3 ;1
;;grow = grow > 1
;order = 15
;
;count = 0
;flag = 0
;WHILE (flag eq 0) do begin
;
;  ; Find deep absorption lines, fit and remove them
;  ;-------------------------------------------------
;  normspec = spec/tcont
;  sm = savgolsm(normspec,[10,10,4])
;  maxima,sm,minarr,maxarr
;
;  ;diff = normspec
;  tspec = spec
;  diff = normspec
;  plot,wave,normspec
;  undefine,allstr
;
;  nmin = n_elements(minarr)
;  for i=0,nmin-1 do begin
;
;    ht_thresh = 0.9
;    if (normspec[minarr[i]] lt ht_thresh) then begin
;
;      speclinefit,wave,normspec,wave[minarr[i]],outstr,/gaussian,/constant,/silent
;      ;speclinefit,wave,normspec,wave[minarr[i]],outstr,/moffat,/constant,/silent
;
;      stop
;
;      ; Must have positive EW
;      if outstr.ew gt 0.0 and outstr.par[0] gt 0.1 and outstr.fwhm lt 100. then begin
;
;        gg = where(abs(wave-outstr.cen) lt 5.0*outstr.fwhm,ngg)
;        w = wave[gg]
;        yy = mpfitpeak_gauss(w,outstr.par)
;        ;yy = mpfitpeak_moffat(w,outstr.par)
;        line = yy-outstr.par[3]  ; remove baseline
;        oplot,w,1.0-yy,co=250
;
;        diff[gg] = diff[gg]+line
;        tspec[gg] = tspec[gg]+line*tcont[gg]
;
;        print,outstr.par
;        push,allstr,outstr
;
;      endif
;
;    end
;
;  end
;
;  stop
;  wait,0.2
;
;  ratio = diff
;
;
;  ; Find outliers using error spectrum
;  ;-------------------------------------
;
;  ; Parameters
;  order =  15 ;20 ;10 ; 20.
;  low_rej = 3 ;4 ;3   ; 0.9 ;2
;  high_rej = 0 ;5  ; 0 ;1.1 ;0
;  ; Reject positive outliers once a decent continuum has been found
;  ;if count ge 5 then high_rej = 3
;  if count ge 5 then high_rej = 6 ;3
;  grow = 0
;
;  ; DON'T GROW!!  This causes problems by throwing out few remaining
;  ; good points b/w close, wide absorption lines
;
;  ; Find outliers
;  sig = errspec/tcont < 0.2      ; 0.2 maximum
;  sig = 0.01 > errspec/tcont     ; 0.02 for synthetic spectrum
;  high_thresh = 1.0+high_rej*sig
;  low_thresh = 1.0-low_rej*sig
;  nbdhi = 0
;  if high_rej gt 0 then bdhi = where(ratio gt high_thresh,nbdhi)
;  if nbdhi gt 0 then tspec[bdhi] = !values.f_nan
;  nbdlo = 0
;  if low_rej gt 0 then bdlo = where(ratio lt low_thresh,nbdlo)
;  if nbdlo gt 0 then tspec[bdlo] = !values.f_nan
;
;  ; Grow
;  if grow gt 0 then begin
;    bad = x*0.
;    if nbdhi gt 0 then bad[bdhi]=1
;    if nbdlo gt 0 then bad[bdlo]=1
;    gkernel = fltarr(long(grow)+2)+1.0
;    bad2 = convol(bad,gkernel,/center,/edge_truncate)
;    bdgrow = where(bad2 gt 0.5,nbdgrow)
;    if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;    ;if nbdgrow gt 0 then tx[bdgrow] = !values.f_nan
;  endif
;
;  ;bd = where(finite(tspec) eq 0,nbd)
;  ;if nbd gt 0 then tx[bd] = !values.f_nan
;
;
;  ; Get the continuum
;  ;---------------------
;  ; Bin the data
;  ;;nbin = ceil( float(nspec)/float(order) )
;  ;;nbin = 10 >  round(100./dw)  ; 100A width
;  ;nbin = 10 >  round(200./dw)  ; 200A width
;  ;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;
;  ; At half steps  
;  ;BINDATA,tx[nbin/2:*],tspec[nbin/2:*],xbin2,specbin2,binsize=nbin,/nan,xmnbin=xmnbin2,min=tx[nbin/2]
;  ;xbinall = [xmnbin,xmnbin2]
;  ;specbinall = [specbin,specbin2]
;  ;si = sort(xbinall)
;  ;xbinall = xbinall[si]  ; put them in order
;  ;specbinall = specbinall[si]
;
;  ; Spline
;  ;;spl_cont = CSPLINE(xbin,specbin,x)
;  ;spl_cont = CSPLINE(xmnbin,specbin,x)
;  ;;spl_cont = CSPLINE(xbinall,specbinall,x)
;  ;tcont = spl_cont
;
;  ; Masked Gaussian smoothing, 100A wide
;  ;smlen = ceil( float(nspec)/float(order) )
;  ;smlen = 10 > ( round(nspec*0.01) > round(100./dw) )  ; 100A width
;  ;smlen = smlen < round(nspec*0.5)
;  ;tcont = GSMOOTH(tspec,smlen)
;  smlen = 10 >  round(100./dw)  ; 100A width
;  ;sm_cont = GSMOOTH(tspec/spl_cont,smlen)
;  tcont = GSMOOTH(tspec,smlen)
;
;  ;tcont = spl_cont * sm_cont
;
;  ;; Normalized spectrum
;  ;ratio = tspec/tcont
;
;
;  ; Finished
;  if count ge niterate then flag=1
;
;  ; Plot continuum normalization
;  !p.multi=[0,1,2]
;  plot,x,spec,tit='Iteration='+strtrim(count,2)
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec[gd],ps=1,co=150
;  oplot,tcont,co=250,thick=2.0
;  ;oplot,xbin,specbin,ps=1,co=250,sym=2
;  oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;  ;oplot,xbinall,specbinall,ps=1,co=250,sym=3,thick=2.5
;  ;oplot,sm_cont2,co=90,thick=2.0
;  ;oplot,spl_cont,co=90,thick=2.0
;  
;  plot,spec/tcont,yr=[0,2.0],ys=1
;  ;plot,tspec/sm_cont,yr=[0,2.0],ys=1
;  oplot,gd,ratio[gd],ps=1,co=150;  ;oplot,minmax(x),[1,1],co=200
;  if high_rej gt 0 then oplot,high_thresh,co=250
;  if low_rej gt 0 then oplot,low_thresh,co=250
;  
;  ;wait,0.1
;  stop
;
;  count++
;
;ENDWHILE
;
;;stop
;
;;pl = 1
;; plotting
;if keyword_set(pl) then begin
;  !p.multi=[0,1,2]
;  plot,spec,xs=1
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec[gd],ps=1,co=150
;  oplot,tcont,co=250,thick=2.0
;  ;oplot,xbin,specbin,ps=1,co=250,sym=2
;  oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;
;  ;!p.multi=[1,1,2]
;  plot,[0],/nodata,xr=minmax(x),yr=[0.0,1.5],xs=1,ys=1
;  oplot,spec/tcont
;  oplot,gd,(spec/tcont)(gd),ps=1,co=150
;  ;oplot,spec/cont4,co=200
;  ;oplot,spec/cont2
;  ;oplot,spec/cont1,co=200
;  oplot,minmax(x),[1,1],co=200
;  if high_rej gt 0 then oplot,high_thresh,co=250
;  if low_rej gt 0 then oplot,low_thresh,co=250
;
;  wait,0.2
;  !p.multi=0
;  ;stop
;endif
;
;cont = tcont
;spec2 = spec/cont
;
;stop
;return


;##########################
; FIFTH ATTEMPT
;##########################

; Hybrid method:
; heavy gaussian smoothing, outlier rejection using the error
; spectrum, and splining from the good points (maybe some binning)
;
; NOTES:
; maybe use the error spectrum (zeros for synthetic)
; to figure out at what fraction of the maximum the
; continuum should be.
;
; Use the error spectrum to set the outlier rejection
; threshold for each pixel.  Set lower and upper limits
; on what the threshold can be.  Input 0 for synthetic
; spectra.
;
; Maybe use upper and lower outlier rejection.
;
; maybe use a hybrid method that uses heavy gaussian
; smoothing and then splining.  The splining is good for
; getting the shape right over the wide lines, but it needs
; more points to work with.  So first use heavy gaussian
; smoothing (with outlier rejection?) to get rid of the
; big lines and then spline with the rest of the "good"
; lines (maybe binned).

;; Make error spectrum
;if n_elements(errspec0) gt 0 then errspec=errspec0 else errspec=dblarr(nspec)
;if n_elements(errspec) eq 1 then errspec=dblarr(nspec)
;
;; Dispersion
;dw = median(slope(wave))
;
;nspec = n_elements(spec)
;x = findgen(nspec)
;tspec = spec
;tx = x
;
;;; Heavy Gaussian smoothing
;;smlen = 10 > ( round(nspec*0.10) > round(200./dw) )
;;smlen = smlen < round(nspec*0.5)
;;smspec = GSMOOTH(tspec,smlen)
;;
;;; Bin the data
;;order = 10
;;nbin = ceil( float(nspec)/float(order) )
;;;BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;;BINDATA,tx,smspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;;
;;; Spline
;;spl_cont = CSPLINE(xbin,specbin,x)
;;spl_cont = CSPLINE(xmnbin,specbin,x)
;
;; Outlier rejection using the error spectrum
;niterate = 10  ;20
;;low_rej = 3   ; 0.9 ;2
;;high_rej = 0 ;5  ; 0 ;1.1 ;0
;;order = 20 ;20 ;10 ; 20.
;;grow = ceil(5.0/dw) ; 3 ;1
;;grow = grow > 1
;order = 15
;
;count = 0
;flag = 0
;WHILE (flag eq 0) do begin
;
;  ; Bin the data
;  nbin = ceil( float(nspec)/float(order) )
;  BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;
;  ; Spline
;  ;spl_cont = CSPLINE(xbin,specbin,x)
;  spl_cont = CSPLINE(xmnbin,specbin,x)
;
;  ; Normalized spectrum
;  ratio = tspec/spl_cont
;
;
;  ; Find outliers using error spectrum
;
;  ; Parameters
;  order =  15 ;20 ;10 ; 20.
;  low_rej = 4 ;3   ; 0.9 ;2
;  high_rej = 0 ;5  ; 0 ;1.1 ;0
;  ; Reject positive outliers once a decent continuum has been found
;  ;if count ge 5 then high_rej = 3
;  if count ge 5 then high_rej = 6 ;3
;  grow = 0
;
;  ; DON'T GROW!!  This causes problems by throwing out few remaining
;  ; good points b/w close, wide absorption lines.
;
;  ; Find outliers
;  sig = errspec/spl_cont < 0.2      ; 0.2 maximum
;  sig = 0.01 > errspec/spl_cont     ; 0.01 for synthetic spectrum
;  high_thresh = 1.0+high_rej*sig
;  low_thresh = 1.0-low_rej*sig
;  nbdhi = 0
;  if high_rej gt 0 then bdhi = where(ratio gt high_thresh,nbdhi)
;  if nbdhi gt 0 then tspec[bdhi] = !values.f_nan
;  nbdlo = 0
;  if low_rej gt 0 then bdlo = where(ratio lt low_thresh,nbdlo)
;  if nbdlo gt 0 then tspec[bdlo] = !values.f_nan
;
;  ; Grow
;  if grow gt 0 then begin
;    bad = x*0.
;    if nbdhi gt 0 then bad[bdhi]=1
;    if nbdlo gt 0 then bad[bdlo]=1
;    gkernel = fltarr(long(grow)+2)+1.0
;    bad2 = convol(bad,gkernel,/center,/edge_truncate)
;    bdgrow = where(bad2 gt 0.5,nbdgrow)
;    if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;    if nbdgrow gt 0 then tx[bdgrow] = !values.f_nan
;  endif
;
;  bd = where(finite(tspec) eq 0,nbd)
;  if nbd gt 0 then tx[bd] = !values.f_nan
;
;  ; Finished
;  if count ge niterate then flag=1  
;
;  ; Plot continuum normalization
;  ;!p.multi=[0,1,2]
;  ;plot,spec,tit='Iteration='+strtrim(count,2)
;  ;gd = where(finite(tspec) eq 1,ngd)
;  ;oplot,gd,spec[gd],ps=1,co=150
;  ;oplot,spl_cont,co=250,thick=2.0
;  ;;oplot,xbin,specbin,ps=1,co=250,sym=2
;  ;oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;  ;
;  ;plot,spec/spl_cont,yr=[0,2.0],ys=1
;  ;oplot,gd,ratio[gd],ps=1,co=150
;  ;oplot,minmax(x),[1,1],co=200
;  ;if high_rej gt 0 then oplot,high_thresh,co=250
;  ;if low_rej gt 0 then oplot,low_thresh,co=250
;  
;  ;wait,0.1
;  ;stop
;
;  count++
;
;ENDWHILE
;;
;
;pl = 1
;; plotting
;if keyword_set(pl) then begin
;  !p.multi=[0,1,2]
;  plot,spec,xs=1
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec[gd],ps=1,co=150
;  oplot,spl_cont,co=250,thick=2.0
;  ;oplot,xbin,specbin,ps=1,co=250,sym=2
;  oplot,xmnbin,specbin,ps=1,co=250,sym=3,thick=2.5
;
;  ;!p.multi=[1,1,2]
;  plot,[0],/nodata,xr=minmax(x),yr=[0.0,1.5],xs=1,ys=1
;  oplot,spec/spl_cont
;  oplot,gd,(spec/spl_cont)(gd),ps=1,co=150
;  ;oplot,spec/cont4,co=200
;  ;oplot,spec/cont2
;  ;oplot,spec/cont1,co=200
;  oplot,minmax(x),[1,1],co=200
;  if high_rej gt 0 then oplot,high_thresh,co=250
;  if low_rej gt 0 then oplot,low_thresh,co=250
;
;  wait,0.2
;  !p.multi=0
;endif
;
;cont = spl_cont
;spec2 = spec/cont
;
;stop
;return


;##########################
; FOURTH ATTEMPT
;##########################
; Algorithm based on the IRAF continuum spline3 method.
; break the spectrum into "order" segments and spline the rest of the
;   array from these.  Then remove outliers and iterate.
;
;; Dispersion
;dw = median(slope(wave))
;
;nspec = n_elements(spec)
;x = findgen(nspec)
;tspec = spec
;tx = x
;
;niterate = 10  ;20
;low_rej = 0.9 ;2
;high_rej = 0 ;1.1 ;0
;order = 15 ;10 ;10 ; 20.
;;grow = ceil(5.0/dw) ; 3 ;1
;;grow = grow > 1
;grow = 0
;
;count = 0
;flag = 0
;WHILE (flag eq 0) do begin
;
;  ; Bin the data
;  nbin = ceil( float(nspec)/float(order) )
;  BINDATA,tx,tspec,xbin,specbin,binsize=nbin,/nan,xmnbin=xmnbin,min=0
;
;  ; Maybe use some kind of "weighted" x values for the binned data
;
;  ; Spline
;  ;spl_cont = CSPLINE(xbin,specbin,x)
;  spl_cont = CSPLINE(xmnbin,specbin,x)
;
;  ; Calculate RMS
;  diff = tspec - spl_cont
;  ratio = tspec/spl_cont
;  ;rms = sqrt(mean(diff^2.0,/nan))
;
;  ; DO NOT use RMS to reject outliers since it is S/N dependent.
;  ; Need a S/N INDEPENDENT means of identifying outliers.
;
;  ; Find outliers
;  nbdhi = 0
;  ;if high_rej gt 0 then bdhi = where(diff gt high_rej*rms,nbdhi)
;  if high_rej gt 0 then bdhi = where(ratio gt high_rej,nbdhi)
;  if nbdhi gt 0 then tspec[bdhi] = !values.f_nan
;  nbdlo = 0
;  ;if low_rej gt 0 then bdlo = where(diff lt -low_rej*rms,nbdlo)
;  if low_rej gt 0 then bdlo = where(ratio lt low_rej,nbdlo)
;  if nbdlo gt 0 then tspec[bdlo] = !values.f_nan
;
;  ; Grow
;  bad = x*0.
;  if nbdhi gt 0 then bad[bdhi]=1
;  if nbdlo gt 0 then bad[bdlo]=1
;  gkernel = fltarr(long(grow)+2)+1.0
;  bad2 = convol(bad,gkernel,/center,/edge_truncate)
;  bdgrow = where(bad2 gt 0.5,nbdgrow)
;  if nbdgrow gt 0 then tspec[bdgrow] = !values.f_nan
;  if nbdgrow gt 0 then tx[bdgrow] = !values.f_nan
;
;  bd = where(finite(tspec) eq 0,nbd)
;  if nbd gt 0 then tx[bd] = !values.f_nan
;
;  ; Finished
;  if count ge niterate then flag=1  
;
;  ;plot,spec
;  ;gd = where(finite(tspec) eq 1,ngd)
;  ;oplot,gd,spec[gd],ps=1,co=150
;  ;oplot,spl_cont,co=250
;  ;;oplot,xbin,specbin,ps=1,co=250,sym=2
;  ;oplot,xmnbin,specbin,ps=1,co=250,sym=2
;  ;
;  ;wait,0.1
;  ;stop
;
;  count++
;
;ENDWHILE
;
;;; Get continuum by masked Gaussian smoothing
;;; want sigmsa=10A, fwhm=2.35*sigma = 23.5A
;;mask = float(finite(tspec))
;;fwhmpix = ceil(nspec/10.)
;;save_except = !EXCEPT & !EXCEPT=0
;;psf = psf_gaussian(npixel=nspec,fwhm=fwhmpix,ndim=1,/norm)
;;;dum = check_math()
;;!EXCEPT = save_except
;;cont4 = CONVOL(spec*mask,psf,invalid=0.0,/center,/edge_truncate,/norm)
;;
;; This continuum is slightly better on the very blue edge
;;; BUT worse in other regions
;
;;cont2 = spl_cont
;
;pl = 1
;; plotting
;if keyword_set(pl) then begin
;  !p.multi=[0,1,2]
;  plot,spec,xs=1
;  gd = where(finite(tspec) eq 1,ngd)
;  oplot,gd,spec[gd],ps=1,co=150
;  oplot,spl_cont,co=250
;;  ;oplot,cont4,co=200
;; ;oplot,xbin,specbin,ps=1,co=250,sym=2
;; oplot,xmnbin,specbin,ps=1,co=250,sym=2
;; ;oplot,cont1,co=200
;
;  ;!p.multi=[1,1,2]
;  plot,[0],/nodata,xr=minmax(x),yr=[0.0,1.2],xs=1,ys=1
;  oplot,spec/spl_cont
;  oplot,gd,(spec/spl_cont)(gd),ps=1,co=150
;  ;oplot,spec/cont4,co=200
;  ;oplot,spec/cont2
;  ;oplot,spec/cont1,co=200
;  oplot,minmax(x),[1,1],linestyle=5,co=200
;  !p.multi=0
;endif
;
;cont = spl_cont
;spec2 = spec/cont
;
;stop


;##########################
; THIRD ATTEMPT
;##########################
;
;; Actually call IRAF continuum
;irafdir = '/net/home/dln5q/iraf/'
;
;; Write the spectrum to temporary fits file
;tempin = MKTEMP('temp')
;tempin2 = tempin+'.fits'
;tempout = MKTEMP('temp')
;tempout2 = tempout+'.fits'
;FITS_WRITE,tempin2,spec
;
;; Make the scripts
;cd,current=curdir
;undefine,iraflines
;push,iraflines,'cd '+curdir
;push,iraflines,'noao'
;push,iraflines,'onedspec'
;push,iraflines,'continuum(input="'+tempin2+'",output="'+tempout2+'",ask=no,lines="*",bands="1",type="fit",'+$
;               'replace=no,wavescale=no,logscale=no,override=no,listonly=no,logfiles="",interactive=no,'+$
;               'sample="*",naverage=1,function="spline3",order=3,low_reject=2.0,high_reject=0.0,niterate=10,'+$
;               'grow=1.0,markrej=no)'
;push,iraflines,'logout'
;tempscript = MKTEmP('cont')
;contscript = curdir+'/'+tempscript
;WRITELINE,tempscript,iraflines
;IRAF_RUN,tempscript,irafdir,out=out,/silent,error=iraferror
;
;
; continuum
;        input = "test"          Input images
;       output = "none"          Output images
;          ask = "yes"                                                
;       (lines = "*")            Image lines to be fit
;       (bands = "1")            Image bands to be fit
;        (type = "ratio")        Type of output
;     (replace = no)             Replace rejected points by fit?
;   (wavescale = yes)            Scale the X axis with wavelength?
;    (logscale = no)             Take the log (base 10) of both axes?
;    (override = no)             Override previously fit lines?
;    (listonly = no)             List fit but don't modify any images?
;    (logfiles = "logfile")      List of log files
; (interactive = yes)            Set fitting parameters interactively?
;      (sample = "*")            Sample points to use in fit
;    (naverage = 1)              Number of points in sample averaging
;    (function = "spline3")      Fitting function
;       (order = 1)              Order of fitting function
;  (low_reject = 2.)             Low rejection in sigma of fit
; (high_reject = 0.)             High rejection in sigma of fit
;    (niterate = 10)             Number of rejection iterations
;        (grow = 1.)             Rejection growing radius in pixels
;     (markrej = yes)            Mark rejected points?
;    (graphics = "stdgraph")     Graphics output device
;      (cursor = "")             Graphics cursor input
;        (mode = "ql")    
;
;; Load the fit
;FITS_READ,tempout2,cont1,/no_abort
;
;; Remove temporary files
;FILE_DELETE,[tempin,tempin2,tempout,tempout2,tempscript],/allow,/quiet

;SKIPTOHERE:

;##########################
; SECOND ATTEMPT
;##########################
;
;; Dispersion
;dw = median(slope(wave))
;
;; Smooth
;smspec0 = SAVGOLSM(spec,[5,5,4])
;smspec1 = SAVGOLSM(spec,[10,10,4])
;std = mad(spec-smspec1)
;
;; IMPORTANT CHANGES NEEDING TO BE MADE!!!
;; 1.) ALWAYS MASK OUT PROMINENT LINES, I.E. Halpha, Hbeta, NaD, MgH, etc.
;; 2.) USE A LARGER GAUSSIAN
;; MAKE SURE THAT THE MAGNESIUM BAND DOES NOT GET REMOVED!!!!
;;print,'NEED TO MAKE IMPORTANT CHANGES!!!!!'
;;stop
;
;; Lines to mask out:
;; Halpha, Hbeta, Hgamma, Hdelta, NaD, MgH+MgB,
;; 
;
;;; Normalize with rebinned maximum
;;sbin = round(40./dw)
;;nbin = ceil(float(nwave)/sbin)
;;maxcont = REBINMAX(smspec0,sbin)*0.95
;;maxx = dindgen(nbin)*sbin+0.5*sbin
;;; Add end points
;;sbinhalf = ceil(sbin*0.5)
;;maxcont = [max(smspec0[0:sbinhalf-1])*0.95, maxcont, max(smspec0[nwave-sbinhalf:nwave-1])*0.95]
;;maxx = [0,maxx,nwave-1]
;;; Polynomial fit
;;nord = (nbin+2-1) < 10
;;coef = poly_fit(maxx,double(maxcont),nord)
;;xx = findgen(nwave)
;;nmaxcont = poly(xx,coef)
;;nsmspec1 = smspec0/nmaxcont
;;
;;std = mad(spec/nmaxcont - savgolsm(spec/nmaxcont,[10,10,4]) )
;;
;;; Now mask out the lines
;;mask1 = where(nsmspec1 gt 0.8 and nsmspec1 lt 1.2,nmask1)
;;bad1 = where(nsmspec1 lt 0.8 or nsmspec1 gt 1.2,nbad1)
;;tempspec = smspec0
;;;tempspec[bad1] = nmaxcont[bad1]
;;tempspec[bad1] = !values.f_nan
;
;; First continuum through vigorous smoothing
;nsm = ceil(50./dw)
;nsm = (nsm > 100.) < nwave
;cont1 = SMOOTH(spec,nsm,/edge_truncate)
;
;std = MAD(spec-cont1)
;
;; reject points far from the continuum
;bad1 = where(abs(spec-cont1) gt 3.0*std,nbad1)
;tempspec = spec
;if nbad1 gt 0 then tempspec[bad1] = !values.f_nan
;
;; smooth again
;cont2 = SMOOTH(tempspec,nsm,/edge_truncate,/nan)
;missed = where(finite(cont2) eq 0,nmissed)
;if nmissed gt 0 then cont2[missed] = cont1[missed]
;
;nspec1 = spec/cont2
;
;; Use masked Gaussian smoothing
;;------------------------------
;
;; First mask
;mask = float(nspec1 ge 0.8 and nspec1 le 1.2)
;
;; Grow the "bad" regions
;ngrow = ceil(5.0/dw)
;gkernel = fltarr(ngrow)+1.0
;invmask = 1.0-mask
;invmask2 = convol(invmask,gkernel,/center,/edge_truncate)
;invmask2 = invmask2/(invmask2 > 1.0)
;mask2 = 1.0-invmask2
;invmask3 = float(mask2 eq 0.0 and (nspec1 lt 0.95 or nspec1 gt 1.05))
;mask3 = 1.0-invmask3
;
;; Get continuum by masked Gaussian smoothing
;; want sigma=10A, fwhm=2.35*sigma = 23.5A
;fwhmpix = 75/dw ;23.5/dw
;save_except = !EXCEPT & !EXCEPT=0
;psf = psf_gaussian(npixel=nwave,fwhm=fwhmpix,ndim=1,/norm)
;dum = check_math()
;!EXCEPT = save_except
;cont3 = CONVOL(spec*mask3,psf,invalid=0.0,/center,/edge_truncate,/norm)
;nspec2 = spec/cont3
;
;; Second try, masked Gaussian smoothing
;;---------------------------------------
;
;; First mask
;mask = float(nspec2 ge 0.9 and nspec2 le 1.1)
;
;; Grow the "bad" regions
;ngrow = ceil(5.0/dw)
;gkernel = fltarr(ngrow)+1.0
;invmask = 1.0-mask
;invmask2 = convol(invmask,gkernel,/center,/edge_truncate)
;invmask2 = invmask2/(invmask2 > 1.0)
;mask2 = 1.0-invmask2
;invmask3 = float(mask2 eq 0.0 and (nspec2 lt 0.95 or nspec2 gt 1.05))
;mask3 = 1.0-invmask3
;
;; Get continuum by masked Gaussian smoothing
;; want sigmsa=10A, fwhm=2.35*sigma = 23.5A
;fwhmpix = 75/dw ;23.5/dw
;save_except = !EXCEPT & !EXCEPT=0
;psf = psf_gaussian(npixel=nwave,fwhm=fwhmpix,ndim=1,/norm)
;dum = check_math()
;!EXCEPT = save_except
;cont4 = CONVOL(spec*mask3,psf,invalid=0.0,/center,/edge_truncate,/norm)
;nspec3 = spec/cont4
;
;; Final normalized spectrum
;cont = cont4
;spec2 = nspec3
;
;
;
;##########################
; FIRST ATTEMPT
;##########################

;; Normalize with rebinned maximum
;sbin2 = sbin ; round(sbin/2.)
;sbin2half = ceil(sbin2*0.5)
;nbin2 = ceil(float(nwave)/sbin2)
;maxcont2 = REBINMAX(tempspec,sbin2,/nan)*0.95
;maxx2 = dindgen(nbin2)*sbin2+0.5*sbin2
;; Get good points
;gd = where(finite(maxcont2) eq 1,ngd)
;maxcont2 = maxcont2[gd]
;maxx2 = maxx2[gd]
;; Add end points
;;maxcont2 = [max(tempspec[0:sbinhalf-1])*0.95, maxcont2, max(tempspec[nwave-sbinhalf:nwave-1])*0.95]
;maxcont2 = [maxcont2[0], maxcont2, maxcont2[ngd-1]]  ; repeat first and last
;maxx2 = [0,maxx2,nwave-1]
;; Polynomial fit
;nord2 = (nbin2+2-1) < 10
;coef2 = poly_fit(maxx2,double(maxcont2),nord2)
;xx2 = findgen(nwave)
;nmaxcont2 = poly(xx,coef2)
;nsmspec2 = smspec0/nmaxcont2
;
;; Final normalized spectrum
;cont = nmaxcont2
;spec2 = spec/nmaxcont2
;
;
; OLD STUFF
;
;; Normalize with CONFT
;sbin = round(nwave*dw/7.)
;CONTF,smspec0,cont1,sbin=sbin,nord=4
;nsmspec1 = smspec0/cont1
;
;; Again, but with higher order
;;sbin = round(nwave*dw/10.)
;CONTF,nsmspec1,ncont2,sbin=sbin,nord=5
;cont2 = cont1/ncont2
;nsmspec2 = smspec0/cont2
;
;; Again, but now exclude lines below 0.8 and above 1.2
;mask1 = where(nsmspec1 gt 0.8 and nsmspec1 lt 1.2,nmask1)
;CONTF,nsmspec1,ncont2,sbin=70,nord=5,mask=mask1
;cont2 = cont1/ncont2
;nsmspec2 = smspec0/cont2
;
;
;std = mad(smspec2[mask1])
;nsmspec = smspec2
;oldcont = cont1
;
;; Iterate until convergence
;flag = 0
;count = 0
;WHILE (flag eq 0) do begin
;
;  hi = (1.1 > (1.0+2.0*std)) < 1.2
;  lo = (0.9 < (1.0-2.0*std)) > 0.8
;  mask = where( nsmspec gt lo and nsmspec lt hi,nmask)
;  CONTF,smspec0,ncont,sbin=70,nord=6,mask=mask
;
;  maxdiff = max(abs(ncont-oldcont)/abs(ncont))
;
;  oldcont = ncont
;  std = mad(smspec0[mask]/ncont[mask])
;
;  if maxdiff lt 0.01 then flag=1
;
;  count++
;
;  stop
;
;ENDWHILE
;
;;; Again, but now exclude lines below 0.8
;;mask1 = where(smspec1 gt 0.8,nmask1)
;;CONTF,smspec0,cont2,sbin=70,nord=5,mask=mask1
;;smspec2 = smspec0/cont2
;;
;;; Again, but now exclude lines below 0.9
;;mask2 = where(smspec2 gt 0.9,nmask2)
;;CONTF,smspec0,cont3,sbin=70,nord=5,mask=mask2
;;smspec3 = smspec0/cont3


end
