pro lickana
;read in lick measurement
lickfile = '/scr2/nichal/workspace2/lickindices/output/results/'+maskname+'_lickmeasurement.fits'
lickstr = mrdfits(lickfile,1,/silent)
;read in sps_measurement
spsfile = '/scr2/nichal/workspace2/sps_fit/data/'+maskname+'/sps_fit.fits.gz'
spsstr = mrdfits(spsfile,1,/silent)

goodobj = where(spsstr.good eq 1 and spsstr.goodfit eq 1)

end
