PRO keck_proposal_specs
	fullzz = mrdfits('/net/marc/raid2/MTR/DEEP2/zcat.egs.uniq.fits.gz', 1)
	zz = fullzz
	objs = [12003798,12007729,12007879,12007892,12007903,12007955,12008087,12008164,12011707,12008288,12012181,12012244]
	match, zz.objno, objs, a, b
	zz = zz[a]

	if (size(ss))[0] eq 0 then ss = mrdfits('/net/lambda/a/data/d2zall_i8_m24_m70_k_d1_ch_gx_.fits',1)

	n = obj_new('npk_specs', zz, secondary=ss,/mtr)
END
