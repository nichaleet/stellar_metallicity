;.r ucsc_zcat_gen
;.r full_zcat_gen
;.r all_zcat_gen
;.r addkcorr
; save, /routines, 'ucsc_zcat_gen', 'full_zcat_gen', 'all_zcat_gen', 'kcorr1', 'kcorr2', filename='catgen_routines.sav'

pro catgen
    restore, 'catgen_routines.sav'
    ucsc_zcat_gen
    full_zcat_gen
    all_zcat_gen
    kcorr1
    spawn, 'cp zcat.in ../kcorr/.'
    cd, '../kcorr/.'
    print 'Running kc ...  Use H0=70, AB mag, filename=zcat.in'
    spawn, './kc'
    spawn, 'rm zcat.in'
    spawn, 'mv kcorr.out ../zcat/zcat.out'
    cd, '../zcat/.'
    kcorr2
end
