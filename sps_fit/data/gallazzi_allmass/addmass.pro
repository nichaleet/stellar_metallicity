pro addmass
	str = mrdfits('sps_fit.fits.gz',1)
	catfile = '/scr2/nichal/workspace2/SDSSdata/gallazzi1/z_mstar_all.dat'
	readcol,catfile,plate,mjd,fiber,m16,m50,m84,z16,z50,z84,age16,age50,age84,comment='#'
	for i=0,n_elements(str)-1 do begin
	   loc = where(plate eq str(i).plate and mjd eq str(i).mjd and fiber eq str(i).fiber,cloc)
	   if cloc ne 0 then str(i).logmstar = m50(loc) else stop
	endfor
    mwrfits, str, 'sps_fit.fits', /create, /silent
    spawn, 'gzip -f sps_fit.fits'

end
