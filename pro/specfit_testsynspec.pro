function specfit_testsynspec,teff,logg,metal,alpha,dir=dir,error=error,stp=stp

;+
;
; SPECFIT_TESTSYNSPEC
;
; This program tests if a synthetic spectrum of given stellar
; parameters exists.
;
; INPUTS:
;  teff    The effective temperature in Kelvin.
;  logg    The surface gravity.
;  metal   The metallicity, [M/H].
;  alpha   The alpha abundance, [alpha/Fe].
;  =dir    The directory where the synthetic spectra live.
;  /stp    Stop at the end of the program
;
; OUTPUTS:
;  test    1 if the synthetic spectrum exists, and 0 if not.
;  =error  The error message if there was one.
;
; USAGE:
;  IDL>test = specfit_testsynspec(5050.,2.50,-1.5,0.0,dir=specdir)
;
; By D.Nidever  Oct 2008
;-

undefine,error,test
test = 0

nteff = n_elements(teff)
nlogg = n_elements(logg)
nmetal = n_elements(metal)
nalpha = n_elements(alpha)
ndir = n_elements(dir)

; Not enough inputs
if nteff eq 0 or nlogg eq 0 or nmetal eq 0 or nalpha eq 0 or ndir eq 0 then begin
  print,'Syntax - test = specfit_testsynspec(teff,logg,metal,alpha,dir=dir,error=error,stp=stp)'
  error = 'Not enough inputs'
  return,-1
endif

; Get the synthetic spectrum name
synthname = mksynthname(teff,logg,metal,alpha,error=namerror)
if n_elements(namerror) gt 0 then begin
  print,'ERROR - ',namerror
  error = namerror
  return,-1
endif

; Does the spectrum exist
filename = dir+'/'+synthname+'.fits'
test = FILE_TEST(filename)

if keyword_set(stp) then stop

return,test

end