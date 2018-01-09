pro makedeimosreg, directory

  files = file_search(directory, 'spec1d*.{fits,fits.gz}', count=c)
  w=where(strmatch(files,'*serendip*') eq 0, cw)
  files = files[w]
  c = n_elements(files)
  objnames = strarr(c)
  ras  = fltarr(c)
  decs = fltarr(c)
  for i=0,c-1 do begin
     hdr = headfits(files(i),ext=1,/silent)
     extensions = strsplit(file_basename(files[i]),'.',/extract)
     objname = extensions[3]
     mask    = extensions[1]
     slit    = extensions[2]
     objnames[i] = slit+', '+objname
     ra  = sxpar(hdr,'RA_OBJ')
     dec = sxpar(hdr, 'DEC_OBJ')
     print, ra,'  ',dec
     ra  = float(strsplit(ra,':',/extract))
     ras[i] = ra[0]*15.+ra[1]/4.+ra[2]*15./3600.
     dec = float(strsplit(dec,':',/extract))  
     decs[i] = dec[0]+dec[1]/60.+dec[2]/3600.
  endfor
  write_ds9_regionfile,ras,decs,comment=objnames,filename=mask+'.reg',color='green'
end
