pro specfit_modelinterp,synstr,temp0,logg0,metal0,alpha0,wave,spec,error=error,$
                        silent=silent,stp=stp,head=head,nearest=nearest

;+
;
; SPECFIT_MODELINTERP
;
; This interpolates the synthetic spectra between grid points
;
; **NOTE** The stellar parameters will be rounded to the nearest 0.001.
;
; INPUTS:
;  synstr    The synthetic grid structure.  It should have TEFF, LOGG, METAL, ALPHA and FILE
;              tags.  TEFF, LOGG, METAL and ALPHA should be 1D arrays of the UNIQUE stellar
;              parameters of the synthetic grid.  FILE should be a 4D string array (or 3D if
;              ALPHA only has one element) of dimensions [Nteff,Nlogg,Nmetal,Nalpha] with
;              the absolute path of the synthetic spectrum FITS file.  It should be an empty
;              string if the spectrum does not exist.
;  temp      The Teff value for the interpolated spectrum.
;  logg      The logg value for the interpolated spectrum.
;  metal     The [M/H] value for the interpolated spectrum.
;  alpha     The [alpha/Fe] value for the interpolated spectrum.
;  /nearest  Use closest model instead of interpolating.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  wave      The wavelength array for the interpolated synthetic spectrum.
;  spec      The flux array for the interpolated synthetic spectrum.
;  =head     The synthetic spectrum header
;  =error    The error message if there was one.
;
; USAGE:
;  IDL>specfit_modelinterp,synstr,5050.,2.55,-0.64,0.1,wave,spec
;
; By D.Nidever  Nov. 2008
;-


nsynstr = n_elements(synstr)
ntemp = n_elements(temp0)
nlogg = n_elements(logg0)
nmetal = n_elements(metal0)
nalpha = n_elements(alpha0)

undefine,wave,spec,head,error

; Not enough inputs
if (nsynstr eq 0 or ntemp eq 0 or nlogg eq 0 or nmetal eq 0 or nalpha eq 0) then begin
  print,'Syntax - specfit_modelinterp,synstr,temp,logg,metal,alpha,wave,spec,nearest=nearest,head=head,error=error,silent=silent'
  return
endif

; Checking the stellar parameters
if (VALID_NUM(temp0,ftemp) eq 0) then begin
  error = temp+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(logg0,flogg) eq 0) then begin
  error = logg+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(metal0,fmetal) eq 0) then begin
  error = metal+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif
if (VALID_NUM(alpha0,falpha) eq 0) then begin
  error = alpha+' IS NOT A NUMBER'
  if not keyword_set(silent) then print,error
  return
endif

; Local variables, round to nearest 0.001
temp = round(double(temp0)*1000.0d0)/1000.0d0
logg = round(double(logg0)*1000.0d0)/1000.0d0
metal = round(double(metal0)*1000.0d0)/1000.0d0
alpha = round(double(alpha0)*1000.0d0)/1000.0d0


; Synthetic grid structure
; SYNSTR must have five tags: file, TEFF, LOGG, METAL, ALPHA
; TEFF, LOGG, METAL and ALPHA should be the arrays of UNIQUE stellar values
; FILE should be a NteffxNloggxNmetalxNalpha string array with the model filenames
;  If the model doesn't exist then file should be an empty string.
type = size(synstr,/type)
; Wrong type
if type ne 8 then begin
  error = 'SYNSTR must be the synthetic grid structure'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'FILE') eq 0 then begin
  error = 'SYNSTR does not have FILE tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'TEFF') eq 0 then begin
  error = 'SYNSTR does not have TEFF tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'LOGG') eq 0 then begin
  error = 'SYNSTR does not have LOGG tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'METAL') eq 0 then begin
  error = 'SYNSTR does not have METAL tag'
  if not keyword_set(silent) then print,error
  return
endif
if tag_exist(synstr,'ALPHA') eq 0 then begin
  error = 'SYNSTR does not have ALPHA tag'
  if not keyword_set(silent) then print,error
  return
endif


; Get the unique stellar parameter arrays
; They should already be unique
teffarr = synstr.teff
loggarr = synstr.logg
metalarr = synstr.metal
alphaarr = synstr.alpha
;file = synstr.file
nteffarr = n_elements(teffarr)
nloggarr = n_elements(loggarr)
nmetalarr = n_elements(metalarr)
nalphaarr = n_elements(alphaarr)

ndim=0
if nteffarr gt 1 then ndim++
if nloggarr gt 1 then ndim++
if nmetalarr gt 1 then ndim++
if nalphaarr gt 1 then ndim++


; **NOTE** This code takes too long to execute.  It's because of
; synstr.file.  The "size" command takes a long time on synstr.file
; and if you copy it to a local variable the copying takes a while.
; It's partly because it's a string array.
;
; Check that FILE has the proper dimensions
; 1-element dimensions on the end are left off, i.e. [10,50,15,1] == [10,50,15]
;;sz = size(synstr.file)
;sz = size(file)
;if sz[0] ne ndim then begin
;  error = 'SYNSTR.FILE does not have the right dimensions: [TEFF, LOGG, METAL, ALPHA]'
;  ;error = 'SYNSTR.FILE must have 4 dimensions: [TEFF, LOGG, METAL, ALPHA]'
;  if not keyword_set(silent) then print,error
;  return
;endif
;if sz[1] ne nteffarr then begin
;  error = 'SYNSTR.FILE[*,0,0,0] must have same number of elements as SYNSTR.TEFF'
;  if not keyword_set(silent) then print,error
;  return
;endif
;if ndim gt 1 then begin
;  if sz[2] ne nloggarr then begin
;    error = 'SYNSTR.FILE[0,*,0,0] must have same number of elements as SYNSTR.LOGG'
;    if not keyword_set(silent) then print,error
;    return
;  endif
;  if ndim gt 2 then begin
;    if sz[3] ne nmetalarr then begin
;      error = 'SYNSTR.FILE[0,0,*,0] must have same number of elements as SYNSTR.METAL'
;      if not keyword_set(silent) then print,error
;      return
;    endif
;    if ndim gt 3 then begin
;      if sz[4] ne nalphaarr then begin
;        error = 'SYNSTR.FILE[0,0,0,*] must have same number of elements as SYNSTR.ALPHA'
;        if not keyword_set(silent) then print,error
;        return
;      endif
;    endif ; ndim>3
;  endif ; ndim>2
;endif ; ndim>1


; Are the parameters outside the grid?
if (temp lt min(teffarr) or temp gt max(teffarr)) then begin
  error = 'TEMP out of bounds.  '+strtrim(min(teffarr),2)+' <= Teff <= '+strtrim(max(teffarr),2)
  if not keyword_set(silent) then print,error
  return
endif
if (logg lt min(loggarr) or logg gt max(loggarr)) then begin
  error = 'LOGG out of bounds.  '+strtrim(min(loggarr),2)+' <= LOGG <= '+strtrim(max(loggarr),2)
  if not keyword_set(silent) then print,error
  return
endif
if (metal lt min(metalarr) or metal gt max(metalarr)) then begin
  error = '[M/H] out of bounds. '+strtrim(min(metalarr),2)+' <= [M/H] <= '+strtrim(max(metalarr),2)
  if not keyword_set(silent) then print,error
  return
endif
if (alpha lt min(alphaarr) or alpha gt max(alphaarr)) then begin
  error = '[alpha/Fe] out of bounds. '+strtrim(min(alphaarr),2)+' <= [alpha/Fe] <= '+strtrim(max(alphaarr),2)
  if not keyword_set(silent) then print,error
  return
endif


;----------------------------------------------------------------
; Step I: Get grid points/lines surrounding the requested value
;----------------------------------------------------------------

mch = machar(/double)
eps = mch.eps

; TEMPERATURE
tongrid = where(abs(teffarr-temp) lt 10.*eps,ntongrid)
if ntongrid gt 0 then begin
  teffbox = teffarr[tongrid[0]]
  teffboxind = tongrid[0]
endif else begin
  ; Get temp1
  gdlo = where(teffarr lt temp,ngdlo)
  indlo = first_el(maxloc(teffarr[gdlo]))
  teff1 = teffarr[gdlo[indlo]]
  ; Get temp2
  gdhi = where(teffarr gt temp,ngdhi)
  indhi = first_el(minloc(teffarr[gdhi]))
  teff2 = teffarr[gdhi[indhi]]
  teffbox = [teff1,teff2]
  teffboxind = [gdlo[indlo],gdhi[indhi]]
endelse
nteff = n_elements(teffbox)

; LOGG
gongrid = where(abs(loggarr-logg) lt 10.*eps,ngongrid)
if ngongrid gt 0 then begin
  loggbox = loggarr[gongrid[0]]
  loggboxind = gongrid[0]
endif else begin
  ; Get logg1
  gdlo = where(loggarr lt logg,ngdlo)
  indlo = first_el(maxloc(loggarr[gdlo]))
  logg1 = loggarr[gdlo[indlo]]
  ; Get logg2
  gdhi = where(loggarr gt logg,ngdhi)
  indhi = first_el(minloc(loggarr[gdhi]))
  logg2 = loggarr[gdhi[indhi]]
  loggbox = [logg1,logg2]
  loggboxind = [gdlo[indlo],gdhi[indhi]]
endelse
nlogg = n_elements(loggbox)

; METAL
mongrid = where(abs(metalarr-metal) lt 10.*eps,nmongrid)
if nmongrid gt 0 then begin
  metalbox = metalarr[mongrid[0]]
  metalboxind = mongrid[0]
endif else begin
  ; Get metal1
  gdlo = where(metalarr lt metal,ngdlo)
  indlo = first_el(maxloc(metalarr[gdlo]))
  metal1 = metalarr[gdlo[indlo]]
  ; Get metal2
  gdhi = where(metalarr gt metal,ngdhi)
  indhi = first_el(minloc(metalarr[gdhi]))
  metal2 = metalarr[gdhi[indhi]]
  metalbox = [metal1,metal2]
  metalboxind = [gdlo[indlo],gdhi[indhi]]
endelse
nmetal = n_elements(metalbox)

; ALPHA
aongrid = where(abs(alphaarr-alpha) lt 10.*eps,naongrid)
if naongrid gt 0 then begin
  alphabox = alphaarr[aongrid[0]]
  alphaboxind = aongrid[0]
endif else begin
  ; Get alpha1
  gdlo = where(alphaarr lt alpha,ngdlo)
  indlo = first_el(maxloc(alphaarr[gdlo]))
  alpha1 = alphaarr[gdlo[indlo]]
  ; Get alpha2
  gdhi = where(alphaarr gt alpha,ngdhi)
  indhi = first_el(minloc(alphaarr[gdhi]))
  alpha2 = alphaarr[gdhi[indhi]]
  alphabox = [alpha1,alpha2]
  alphaboxind = [gdlo[indlo],gdhi[indhi]]
endelse
nalpha = n_elements(alphabox)


; Number of models
nmodels = nteff*nlogg*nmetal*nalpha

; The model already exists
if nmodels eq 1 then begin
  file = synstr.file[teffboxind,loggboxind,metalboxind,alphaboxind]
  ; We have a filename
  if strtrim(file,2) ne '' then begin
    test = FILE_TEST(file)
    ; File exists
    if test eq 1 then begin
      if not keyword_set(silent) then print,'Synthetic spectrum already exists.  Restoring '+file
      RDFITS_SPEC,file,wave,spec,head=head
      return
    endif else begin
      if not keyword_set(silent) then print,'Spectrum '+file+' NOT FOUND'
    endelse
  endif else begin
    if not keyword_set(silent) then print,'NO spectrum at this grid point'
  endelse
endif


;----------------------------------------------------------------
; Step II: Get the names and check if they exist
;----------------------------------------------------------------
filenames = strarr(nmodels)
teffvals = fltarr(nmodels)
loggvals = fltarr(nmodels)
metalvals = fltarr(nmodels)
alphavals = fltarr(nmodels)

count = 0
; Loop through Teff values
for i=0,nteff-1 do begin

  ; Loop through logg values
  for j=0,nlogg-1 do begin

    ; Loop through metal values
    for k=0,nmetal-1 do begin

      ; Loop through alpha values
      for l=0,nalpha-1 do begin

        ; Is there a filename
        file = synstr.file[teffboxind[i],loggboxind[j],metalboxind[k],alphaboxind[l]]
        file = strtrim(file,2)

        if file eq '' then begin
          error = 'Teff='+strtrim(temp,2)+' logg='+strtrim(logg,2)+' [M/H]='+strtrim(metal,2)+' [alpha/Fe]='+strtrim(alpha,2)+$
                  ' is OUTSIDE of the grid'
          if not keyword_set(silent) then print,error
          return
        endif

        ; Does the file exist        
        test = FILE_TEST(file)
        if test eq 0 then begin
          error = file+' NOT FOUND'
          if not keyword_set(silent) then print,error
          return
        endif

        filenames[count] = file
        teffvals[count] = teffbox[i]
        loggvals[count] = loggbox[j]
        metalvals[count] = metalbox[k]
        alphavals[count] = alphabox[l]

        count++

      end ; alpha loop
    end ; metal loopg
  end  ; logg loop
end  ; teff loop


;----------------------------------------------------------------
; Step III: Restore the spectra
;----------------------------------------------------------------
if not keyword_set(nearest) then begin

  for i=0,nmodels-1 do begin
    undefine,wave,spec,head
    RDFITS_SPEC,filenames[i],wave,spec,head=head

    if i eq 0 then begin
      n = n_elements(wave)
      wave0 = wave
      specarr = fltarr(nmodels,n)
      specarr[0,*] = spec
    endif else begin

      ; Not the right size
      if n_elements(spec) ne n then begin
        error = filenames[i]+' does NOT have the right size'
        if not keyword_set(silent) then print,error
        return
      endif

      specarr[i,*] = spec
    endelse

  end

end ; not /nearest


;----------------------------------------------------------------
; Step IV: Interpolate based in distance
;----------------------------------------------------------------

; Calculate distances
dteff = range(teffbox)
if dteff eq 0.0 then dteff=1.0   ; in case there's only one point
dlogg = range(loggbox)
if dlogg eq 0.0 then dlogg=1.0
dmetal = range(metalbox)
if dmetal eq 0.0 then dmetal=1.0
dalpha = range(alphabox)
if dalpha eq 0.0 then dalpha=1.0

dist = sqrt( ((teffvals-temp)/dteff)^2.0 + ((loggvals-logg)/dlogg)^2.0 + $
             ((metalvals-metal)/dmetal)^2.0 + ((alphavals-alpha)/dalpha)^2.0 )

; Interpolate
if not keyword_set(nearest) then begin

  ; Take weighted average
  ; weight by 1/distance
  wt = 1.0/dist
  wt = wt/total(wt)  ; normalize
  spec = reform( wt#specarr )

  ; Update header
  SXADDPAR,head,'TEFF',temp
  SXADDPAR,head,'LOGG',logg
  SXADDPAR,head,'FEH',metal
  SXADDPAR,head,'ALPHA_FE',alpha

; Nearest
endif else begin
  best = first_el(minloc(dist))
  undefine,wave,spec,head
  RDFITS_SPEC,filenames[best[0]],wave,spec,head=head
endelse

if keyword_set(stp) then stop

end
