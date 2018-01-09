pro addtag
  old = mrdfits('sps_fit.fits.gz',1)
  new = []
  nstr = n_Elements(old)
  for i=0,nstr-1 do begin
     sci = old[i]
     sci = create_Struct('alphafeupper',-999d,'alphafelower',-999d,sci)
     new = [new,sci]
  endfor
  mwrfits,new,'sps_fit.fits',/create,/silent
  spawn, 'gzip -f sps_fit.fits'

end
