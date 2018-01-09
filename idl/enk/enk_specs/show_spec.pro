PRO show_spec, id
	COMMON zz, ss

	if (size(zz))[0] eq 0 then begin
		zz = mrdfits('/net/marc/raid2/MTR/DEEP2/zcat.egs.uniq.fits.gz', 1)
		ss = mrdfits('/usr/local/www/kai/DEEP/EGS/EGS_CATMATCH/D2ZALL/d2zall_i8_m24_m70_k_d1_ch_gx_acs_.fits',1)
	endif

	r = where(zz.objno eq id)

	if count ne 1 then begin
		print, 'FOUND NO MATCHES'
		return
	endif
	n = obj_new('npk_specs', zz[r], secondary=ss, /mtr)
END
