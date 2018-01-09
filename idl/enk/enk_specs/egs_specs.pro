PRO egs_specs
	zz = read_zcat()
	ss = mrdfits('/usr/local/www/kai/DEEP/EGS/EGS_CATMATCH/D2ZALL/d2zall_i8_m24_m70_k_d1_ch_gx_acs_.fits',1)
	n = obj_new('npk_specs', zz[*], secondary=ss, /mtr)
END
