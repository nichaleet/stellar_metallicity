PRO mtr_specs 
	fullzz = mrdfits('/net/marc/raid2/MTR/DEEP2/zcat.egs.uniq.fits.gz', 1)
	zz = fullzz

	m1 = 0.656  & b1 = -88.248
	m2 = -0.487 & b2 = 156.865

	plot,zz.ra,zz.dec,psym=3
	rg = where(zz.ra ge 214.447, complement=rlt)
	bot = zz.ra
	bot[rg] = zz[rg].ra * m1 + b1
	bot[rlt] = zz[rlt].ra * m2 + b2

	oplot, zz.ra, bot, psym=3,color=1

	rg = where(zz.ra le 214.494, complement=rlt)
	top = zz.ra
	top[rg] = zz[rg].ra * m1 + b1 + 0.089
	top[rlt] = zz[rlt].ra * m2 + b2 + 0.136

	oplot, zz.ra, top,color=1,psym=3

	roi = where(zz.dec le top and zz.dec ge bot and zz.zquality ge 3, count)

	print, 'There are ', count, 'objects found in the MTR.'
	
	zz = zz[roi]

	ss = mrdfits('/usr/local/www/kai/DEEP/EGS/EGS_CATMATCH/D2ZALL/d2zall_i8_m24_m70_k_d1_ch_gx_acs_.fits',1)
	;r = where(ss.i8_irac1 ne -999 and ss.m24_mips24 ne -999)
	;ss = ss[r]
	;match,zz.objno, ss.d2z_objno, a, b

	;n = obj_new('npk_specs', zz, /mtr)
	n = obj_new('npk_specs', zz[*], secondary=ss, /mtr)
END
