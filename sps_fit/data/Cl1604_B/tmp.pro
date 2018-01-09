pro tmp 

tmp = mrdfits('sps_fit.fits.gz',1)
real = mrdfits('real_sps_fit.fits.gz',1,hdr)
struct_add_field,real,'zmag',tmp.zmag
struct_Add_field,real,'gooddeimos',real.goodsky
struct_add_field,real,'goodmosfire',real.good
mwrfits,real,'sps_fit.fits.gz',hdr,/create

stop
end
