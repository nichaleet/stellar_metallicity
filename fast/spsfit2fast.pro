pro spsfit2fast,scidir,nbegin=nbegin,nfinal=nfinal
;info: convert sps_fit science structure's spectra to fast spectra catalog
;input: 
  ;scidir - the directory name that contains sps_fit.fits.gz to be printed out
  ;nbegin - optional, the first spec in the structure to be printed out
  ;nfinal - optional, the last "
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  sciall = mrdfits('/scr2/nichal/workspace2/sps_fit/data/'+scidir+'/sps_fit.fits.gz',1,hdr)

  if n_Elements(nbegin) eq 0 then nbegin = 0
  if n_Elements(nfinal) eq 0 then nfinal = n_Elements(sciall)-1
  outfile = scidir+'_'+string(strtrim(fix(nbegin),2))+'_'+string(strtrim(fix(nfinal),2))+'.spec'

  sciall = sciall[nbegin:nfinal]
  nspecs = nfinal-nbegin+1

  nlambda  = n_Elements((sciall.lambda)[*,0])
  minlambda = max((sciall.lambda)[0,*])
  maxlambda = min((sciall.lambda)[nlambda-1,*])
  dlam = median(sciall.dlam)

  wlmain = sciall[0].lambda
  goodwl = where(wlmain gt minlambda and wlmain lt maxlambda,cgoodwl)
  wlmain = wlmain(goodwl)
  
  fmt = '(I5,1X,F9.2,1X,I1'
  header = '#    bin   wl   tr'
  array = fltarr(nspecs*2+3,cgoodwl)  ;column, row
  array[1,*] = wlmain 
  array[2,*] = 1
  for i=0,nspecs-1 do begin
     wlnow = sciall[i].lambda
     specnow = sciall[i].spec*100.  ;for SDSS change from E-17 erg to E-19 erg/s/cm2/A
     ivarnow = sciall[i].ivar
     specout = interpol(specnow,wlnow,wlmain)
     ivarout = interpol(ivarnow,wlnow,wlmain)
     errout = 1./sqrt(ivarout)*100.
     badwl = where(~finite(errout),cbadwl)
     if cbadwl gt 0 then errout(badwl) = -99.
     array[i*2+3,*]=specout
     array[i*2+4,*]=errout
     fmt = fmt+',1X,F12.4,1X,F12.4'
     header = header+'   F'+strtrim(string(fix(i+nbegin)),2)+'   E'+strtrim(string(fix(i+nbegin)),2)
  endfor
  fmt = fmt+')'

  cmlt = 0 ;cumulative lambda diff
  bin  = 0 
  for i=0,cgoodwl-2 do begin
     cmlt = cmlt+(wlmain(i+1)-wlmain(i))
     if cmlt gt dlam then begin
        bin = bin+1
        cmlt= 0
     endif
     array[0,i]=bin
  endfor
  array[0,cgoodwl-1] = array[0,cgoodwl-2]

  ;writearray  
  openw,1,outfile,width=400
  printf,1,header
  for i=0,cgoodwl-1 do printf,1,array[*,i],format=fmt
  close,1
end
