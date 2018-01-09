pro addtag

  ;old = mrdfits('sps_fit_polycontn4to7_nomg.fits.gz',1)
  old = mrdfits('sps_fit.fits.gz',1)
  new = []
  npix = 1973
  for i=0,9 do begin
     sci = old[i]
     ;sci = create_struct('spsspecmg',dblarr(npix),'spsspecfullmg',dblarr(npix),'spscontfullmg',dblarr(npix),'spsspecfe',dblarr(npix),'spsspecfullfe',dblarr(npix),'spscontfullfe',dblarr(npix),'fitMgmask',bytarr(npix),'fitFemask',bytarr(npix),'mg',-999d,'mgerr',-999d,'fe',-999d,'feerr',-999d,sci)
     ;sci = create_struct('npoly',0,sci)
     sci = create_Struct('chisqmg',-999d,'chisqfe',-999d,sci)
     new = [new,sci]
  endfor
  mwrfits,new,'sps_fit.fits',/create,/silent
  spawn, 'gzip -f sps_fit.fits'

end
