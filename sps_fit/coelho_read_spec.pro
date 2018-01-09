FUNCTION read_cspec,file
  flux = readfits(file) ; should be an array of 52001 pixels
  npix = n_elements(flux)
  lambda = findgen(npix)*0.2+3000. ;angstrom
  
  basename = file_basename(file)
  basename_splt = strsplit(basename,'.',/extract)
  agegyr = float(strmid(basename_splt[4],0,strpos(basename_splt[4],'gyr')))
  model = basename_splt[1]
  
  if strmatch(model,'m05p00') then parr = [-0.5,0.,-0.5]
  if strmatch(model,'m05p04') then parr = [-0.5,0.4,-0.2]
  if strmatch(model,'p00p00') then parr = [0.,0.,0.]
  if strmatch(model,'p00p04') then parr = [0.,0.4,0.3]
  if strmatch(model,'p02p00') then parr = [0.2,0.,0.2]
  if strmatch(model,'p02p04') then parr = [0.2,0.4,0.5]
  
  str = {feh:parr[0],afe:parr[1],zh:parr[2],agegyr:agegyr,spec:flux,lambda:lambda}
  return,str
END

FUNCTION coelho_read_spec,file

  ct = n_elements(file)
  ff = file_search(file[0],count=ctf)
  IF ctf EQ 0 then begin
     print,'READ_SPEC ERROR: file not found: ',file
     return,0
  ENDIF
  
  astr = read_cspec(file[0])
  IF ct GT 1 THEN BEGIN 
     str = replicate(astr[0],ct)
     str[0] = astr
     FOR i=1,ct-1 do str[i] = read_cspec(file[i])
  ENDIF ELSE str = astr

  RETURN, str
END
