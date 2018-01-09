pro makehst_region
  files= file_search('./hst*.fits')
  filename=strmid(files,9,2)
  nfiles = n_Elements(files)
  ra = fltarr(nfiles)
  dec= fltarr(nfiles)
  for i=0,nfiles-1 do begin
     hdr = headfits(files(i),ext=2)
     ra(i) = sxpar(hdr,'CRVAL1')
     dec(i) = sxpar(hdr,'CRVAL2')
  endfor
write_ds9_regionfile,ra,dec,comment=filename,filename='ms0451_hst_layout.reg'

end
