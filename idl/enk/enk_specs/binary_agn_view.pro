PRO binary_agn_view
	fullzz =  read_zcat()
	zz = fullzz
	r = where(zz.objno eq 13025437)

	n = obj_new('npk_specs', zz[r], secondary=ss, /mtr)
END
