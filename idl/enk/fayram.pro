pro go
    readcol, '/net/schechter/a/ekirby/DEEP/ytoradec/serradec.txt', format='L,X,X,I,X,X', objno, yoffset
    zc = read_zcat()
    match, objno, zc.objno, w1, w2
    zc = zc[w2]	
    ss = mrdfits('/usr/local/www/kai/DEEP/EGS/EGS_CATMATCH/D2ZALL/d2zall_i8_m24_m70_k_d1_ch_gx_acs_.fits', 1)
    match, zc.objno, ss.d2z_objno, w1, w2
    ss = ss[w2[sort(w1)]]
    zc = zc[w1[sort(w1)]]
    if not(array_equal(ss.d2z_objno, zc.objno)) then message, 'objno mismatch'
    print, 'total objects: ', n_elements(zc)
    n = obj_new('npk_specs', zc, secondary=ss, /mtr)
    ds9 = obj_new('ds9')
    enk_stamp, objno, ds9
end
