function phot_iso, science, sdss=sdss, dm=dm
    n = n_elements(science)
    if ~keyword_set(dm) then dm = dblarr(n)
    ebv = dblarr(n)
    wyy = replicate(34, n)
    dmerr = replicate(0.0, n)

    ;http://physwww.mcmaster.ca/~harris/mwgc.dat
    w = where(strmid(science.mask, 0, 4) eq 'N147', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.43
        dmerr[w] = 0.04
        ebv[w] = 0.03
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq 'N185', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.23
        dmerr[w] = 0.03
        ebv[w] = 0.03
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq 'n288' or strmid(science.mask, 0, 3) eq '288', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 14.83
        ebv[w] = 0.03
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'vv124', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 25.67
        dmerr[w] = 0.04
        ebv[w] = 0.015
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 6) eq 'ng1904' or strmid(science.mask, 0, 4) eq '1904', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.59
        ebv[w] = 0.01
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n4590' or strmid(science.mask, 0, 4) eq '4590', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.21
        ebv[w] = 0.05
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 6) eq 'ng5024', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 16.32
        ebv[w] = 0.02
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 6) eq 'ng5053', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 16.19
        ebv[w] = 0.04
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n5634', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 17.16
        ebv[w] = 0.05
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 6) eq 'ng5904', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 14.46
        ebv[w] = 0.03
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n1904', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.59
        ebv[w] = 0.01
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n2419', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 19.83
        ebv[w] = 0.08
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n5272', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.07
        ebv[w] = 0.01
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq '5897', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.76
        ebv[w] = 0.09
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n6205', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 14.48
        ebv[w] = 0.02
        wyy[w] = 34
    endif
    w = where(science.mask eq '6229a', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 17.45
        ebv[w] = 0.01
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n6341' or strmid(science.mask, 0, 4) eq '6341', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 14.65
        ebv[w] = 0.02
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n6656b', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 13.60
        ebv[w] = 0.34
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq '6779', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.62    ;Hatzidimitriou et al. (2004)
        ebv[w] = 0.32
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n6838', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 13.79
        ebv[w] = 0.25
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq '6864', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 17.09
        ebv[w] = 0.16
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n7006' or science.mask eq '7006a', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 18.23
        ebv[w] = 0.05
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n7078' or science.mask eq 'n7078d' or science.mask eq 'n7078e' or strmid(science.mask, 0, 4) eq '7078', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.37
        ebv[w] = 0.10
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n7089' or strmid(science.mask, 0, 4) eq '7089', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 15.50
        ebv[w] = 0.06
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n7099' or strmid(science.mask, 0, 4) eq '7099', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 14.64
        ebv[w] = 0.03
        wyy[w] = 34
    endif
    w = where(science.mask eq 'n7492', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 17.06
        ebv[w] = 0.00
        wyy[w] = 34
    endif
    w = where(science.mask eq 'pal2', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 21.05
        ebv[w] = 1.24
        wyy[w] = 34
    endif
    w = where(science.mask eq 'pal13', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 17.21
        ebv[w] = 0.05
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 6) eq 'Pal13-', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 16.928
        dmerr[w] = 0.0
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'pal14', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 19.54
        ebv[w] = 0.04
        wyy[w] = 34
    endif
    w = where(strmid(strlowcase(science.mask), 0, 3) eq 'eri', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 19.77
        dmerr[w] = 0.20
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 4) eq 'LeoI' or strmid(science.mask, 0, 3) eq 'LIN', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 22.02
        dmerr[w] = 0.13
        ;ebv[w] = 0.037
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 2) eq 'L2', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 21.70
        dmerr[w] = 0.10
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'UMa2', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 17.5258
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'LeoT', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 23.1007
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'UMaI', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.1265
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'leo4', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.94
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 2) eq 'CB', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 18.2173
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'CVnII', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.8949
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'CVn1', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 21.62
        dmerr[w] = 0.06
        ebv[w] = 0.0
        wyy[w] = 28
    endif            
    w = where(strmid(science.mask, 0, 4) eq 'Herc', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.6994
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'BooI', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.11
        dmerr[w] = 0.08
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'Boo2' or strmid(science.mask, 0, 5) eq 'booii', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 18.12
        dmerr[w] = 0.20
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strlowcase(strmid(science.mask, 0, 6)) eq 'segue1' or strlowcase(strmid(science.mask, 0, 7)) eq 'segtide', c)
    if c gt 0 then begin
        alpha = 0
        ;dm[w] = 17.049666   ;25.7 +/- 2.5 kpc
        dm[w] = 16.81
        dmerr[w] = 0.2
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strlowcase(strmid(science.mask, 0, 5)) eq 'seg2_' or strlowcase(strmid(science.mask, 0, 3)) eq 'ari', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 17.7       ;Belokurov et al. (2009)
        dmerr[w] = 0.1
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'scl', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.67   ;Pietrzynski et al. (2008)
        dmerr[w] = 0.12
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'for', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.72   ;Rizzi et al. (2007)
        dmerr[w] = 0.04
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'sex', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.90
        dmerr[w] = 0.06
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'umi' or strmid(science.mask, 0, 3) eq 'uss', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.18
        dmerr[w] = 0.12
        ebv[w] = 0.03
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'dra', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.84
        dmerr[w] = 0.14
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'n6822', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 23.40          ;Feast et al. (2012)
        dmerr[w] = 0.05
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'i1613', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.34           ;Tammann et al. (2011)
        dmerr[w] = 0.03
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'peg', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.87           ;Tammann et al. (2011)
        dmerr[w] = 0.06
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'leoa', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.59           ;Tammann et al. (2011)
        dmerr[w] = 0.03
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'aqr', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 25.15           ;McConnachie et al. (2006)
        dmerr[w] = 0.08
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'cet', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.39           ;McConnachie et al. (2005)
        dmerr[w] = 0.07
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'sagdi', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 25.14           ;Momany et al. (2002)
        dmerr[w] = 0.18
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'ter5', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 22.52
        ebv[w] = 2.15
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 7) eq '208BoSB', c)
    if c gt 0 then begin
        alpha = 1
        dm[w] = 19.11            ;Hughes, Wallerstein, & bossi (2008)
        ebv[w] = 0.02
        wyy[w] = 34
    endif
    ;from Kevin Schlaufman, 10/23/2009:
    ;ECHOS1 = 10.5 ± 5 kpc
    ;ECHOS2 = 11 ± 2 kpc 
    w = where(strmid(science.mask, 0, 26) eq 'hectospec_candidates_list1', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 15.1059
        dmerr[w] = 1.03403
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 26) eq 'hectospec_candidates_list2', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 15.2070
        dmerr[w] = 0.394813
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 26) eq 'hectospec_candidates_list3', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 15.49
        dmerr[w] = 0.0
        ebv[w] = 0.06
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 26) eq 'hectospec_candidates_list7', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 10.2
        dmerr[w] = 0.5
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 5) eq 'Long1', c)
    if c gt 0 then begin
        w = where(strmid(science.mask, 0, 5) eq 'Long1' and strtrim(science.objname, 2) eq 'Callisto', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 0.0
            ebv[w] = 0.0
            wyy[w] = 28
        endif
        w = where(strmid(science.mask, 0, 5) eq 'Long1' and strtrim(science.objname, 2) eq 'Arcturus', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 0.257692
            ebv[w] = 0.0
            wyy[w] = 28
        endif
    endif
    w = where(strmid(science.mask, 0, 8) eq 'LVMslits', c)
    if c gt 0 then begin
        w = where(strmid(science.mask, 0, 8) eq 'LVMslits' and strmid(strtrim(science.objname, 2), 0, 7) eq 'v368Her', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 16.139
            ebv[w] = 0.0
            wyy[w] = 28
        endif
    endif
    w = where(strmid(science.mask, 0, 8) eq 'LVMslitB', c)
    if c gt 0 then begin
        w = where(strmid(science.mask, 0, 8) eq 'LVMslitB' and strmid(strtrim(science.objname, 2), 0, 6) eq 'AqrNew', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 18.8
            ebv[w] = 0.0
            wyy[w] = 28
        endif
        w = where(strmid(science.mask, 0, 8) eq 'LVMslitB' and strmid(strtrim(science.objname, 2), 0, 3) eq 'Her', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 18.625
            ebv[w] = 0.0
            wyy[w] = 28
        endif
        w = where(strmid(science.mask, 0, 8) eq 'LVMslitB' and strmid(strtrim(science.objname, 2), 0, 8) eq 'HerClump', c)
        if c gt 0 then begin
            alpha = 0
            dm[w] = 19.205
            ebv[w] = 0.0
            wyy[w] = 28
        endif
    endif
    w = where(strmid(science.mask, 0, 8) eq 'LVMslits', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 0.0
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(science.mask eq 'HerClump', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.205
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(science.mask eq 'Her3' or science.mask eq 'Her4', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 18.625
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(science.mask eq 'Aqr', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 18.8
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 2) eq 'W1', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 17.898918   ;38 +/- 7 kpc
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'HyaII', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 20.64      ;Martin et al. (2015)
        dmerr[w] = 0.16
        ebv[w] = 0.053
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 5) eq 'PscII', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 21.3       ;Belokurov et al. (2010)
        dmerr[w] = 0.3
        ebv[w] = 0.0
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 4) eq 'CrtI', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 21.81      ;Laevens et al. (2014)
        dmerr[w] = 0.5
        ebv[w] = 0.026
        wyy[w] = 28
    endif
    w = where(strmid(science.mask, 0, 3) eq 'd1_', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.36
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 3) eq 'd2_', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 24.17
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 3) eq 'VSS', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 16.087
        dmerr[w] = 0.105
        ebv[w] = 0.0
        wyy[w] = 34
    endif
    w = where(strmid(science.mask, 0, 3) eq 'Cnc', c)
    if c gt 0 then begin
        alpha = 0
        dm[w] = 19.7
        dmerr[w] = 0.2
        ebv[w] = 0.0
        wyy[w] = 28
    endif

    if keyword_set(sdss) then begin
        science.phot_color_sdss = '  '
        science.feh_padova_sdss = -999d
        science.teff_padova_sdss = -999d
        science.logg_padova_sdss = -999d
        science.feherr_padova_sdss = 0d
        science.tefferr_padova_sdss = 0d
        science.loggerr_padova_sdss = 0d
    endif

    if ~keyword_set(sdss) then begin
        science.phot_color = '  '
        science.feh_van = -999d
        science.teff_van = -999d
        science.logg_van = -999d
        science.feherr_van = 0d
        science.tefferr_van = 0d
        science.loggerr_van = 0d

        science.feh_yy = -999d
        science.teff_yy = -999d
        science.logg_yy = -999d
        science.feherr_yy = 0d
        science.tefferr_yy = 0d
        science.loggerr_yy = 0d

        science.feh_padova = -999d
        science.teff_padova = -999d
        science.logg_padova = -999d
        science.feherr_padova = 0d
        science.tefferr_padova = 0d
        science.loggerr_padova = 0d

        science.teff_ram = -999d
        science.tefferr_ram = 0d
        teffgoodram = 1
        science.teff_alonso = -999d
        science.tefferr_alonso = 0d
        teffgoodalonso = 1
    endif

    nscience = n_elements(science)
    if keyword_set(sdss) then begin
        for i=0,nscience-1 do begin
            case 1 of
                science[i].gsdss gt 0 and science[i].isdss gt 0 and science[i].gsdss lt 50 and science[i].isdss lt 50: science[i].phot_color_sdss = 'gi'
                science[i].gsdss gt 0 and science[i].rsdss gt 0 and science[i].gsdss lt 50 and science[i].rsdss lt 50: science[i].phot_color_sdss = 'gr'
                else:
            endcase
        endfor

        w = where(science.phot_color_sdss eq 'gi', c)
        if c gt 0 then begin
            wnoerror = where(science.isdsserr le 0 or science.isdsserr ge 10 or science.gsdsserr le 0 or science.gsdsserr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)
            
            if cnoerror gt 0 then begin
                fehphot = feh(science[wnoerror].isdss, science[wnoerror].gsdss-science[wnoerror].isdss, /gisdss, /padova, dm=dm[wnoerror], ddm=dmerr[wnoerror])
                science[wnoerror].feh_padova_sdss = fehphot.feh
                science[wnoerror].teff_padova_sdss = fehphot.teff
                science[wnoerror].logg_padova_sdss = fehphot.logg
            endif
            if cwitherror gt 0 then begin
                fehphot = feh(science[wwitherror].isdss, science[wwitherror].gsdss-science[wwitherror].isdss, science[wwitherror].isdsserr, sqrt((science[wwitherror].gsdsserr)^2.+(science[wwitherror].isdsserr)^2.), /gisdss, /padova, dm=dm[wwitherror], ddm=dmerr[wwitherror])
                science[wwitherror].feh_padova_sdss = fehphot.feh
                science[wwitherror].teff_padova_sdss = fehphot.teff
                science[wwitherror].logg_padova_sdss = fehphot.logg
                science[wwitherror].feherr_padova_sdss = fehphot.err_feh
                science[wwitherror].tefferr_padova_sdss = fehphot.err_teff
                science[wwitherror].loggerr_padova_sdss = fehphot.err_logg
            endif
        endif

        w = where(science.phot_color_sdss eq 'gr', c)
        if c gt 0 then begin
            wnoerror = where(science.rsdsserr le 0 or science.rsdsserr ge 10 or science.gsdsserr le 0 or science.gsdsserr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)
            
            if cnoerror gt 0 then begin
                fehphot = feh(science[wnoerror].rsdss, science[wnoerror].gsdss-science[wnoerror].rsdss, /grsdss, /padova, dm=dm[wnoerror], ddm=dmerr[wnoerror])
                science[wnoerror].feh_padova_sdss = fehphot.feh
                science[wnoerror].teff_padova_sdss = fehphot.teff
                science[wnoerror].logg_padova_sdss = fehphot.logg
            endif
            if cwitherror gt 0 then begin
                fehphot = feh(science[wwitherror].rsdss, science[wwitherror].gsdss-science[wwitherror].rsdss, science[wwitherror].rsdsserr, sqrt((science[wwitherror].gsdsserr)^2.+(science[wwitherror].rsdsserr)^2.), /grsdss, /padova, dm=dm[wwitherror], ddm=dmerr[wwitherror]) 
                science[wwitherror].feh_padova_sdss = fehphot.feh
                science[wwitherror].teff_padova_sdss = fehphot.teff
                science[wwitherror].logg_padova_sdss = fehphot.logg
                science[wwitherror].feherr_padova_sdss = fehphot.err_feh
                science[wwitherror].tefferr_padova_sdss = fehphot.err_teff
                science[wwitherror].loggerr_padova_sdss = fehphot.err_logg
            endif
        endif

        ;color = science.u-science.b
        ;colorerr = sqrt(science.uerr^2. + science.berr^2.)
        ;theta = 0.5883 + 0.2008*color - 5.931e-3*color^2.0 + 5.319e-3*color*science.feh_padova_sdss[69] - 1.000e-1*science.feh_padova_sdss[69] - 1.542e-2*(science.feh_padova_sdss[69])^2.0
        ;thetaerrclr = (0.2008 - 2.*5.931e-3*color * 5.319*science.feh_padova_sdss[69])*colorerr
        ;thetaerrfeh = (5.319e-3*color - 1.000e-1 - 2.*1.542e-2*science.feh_padova_sdss[69])*science.feherr_padova_sdss[69]
        ;science.teff_alonso = 5040.0 / theta
        ;science.tefferr_alonso = sqrt(((5040./theta^2.)^2.)*(thetaerrclr^2. + thetaerrfeh^2.) + 110.0^2.)
        return, science
    endif

    teffgoodram = bytarr(nscience)
    teffgoodalonso = bytarr(nscience)
    syserrram = dblarr(nscience)
    syserralonso = dblarr(nscience)
    color = dblarr(nscience)
    colorerr = dblarr(nscience)
    aramdwarf = dblarr(nscience, 6)
    aramgiant = dblarr(nscience, 6)
    aalonso = dblarr(nscience, 6)
    clrrangeramdwarf = dblarr(nscience, 2)
    clrrangeramgiant = dblarr(nscience, 2)
    clrrangealonso = dblarr(nscience, 2)
    pramdwarf = dblarr(nscience, 7)
    pramgiant = dblarr(nscience, 7)
    tags = tag_names(science)
    if contains(tags, 'MERR') and contains(tags, 'TERR') then begin
        for i=0,nscience-1 do begin
            if science[i].m gt 0 and science[i].t gt 0 and science[i].m lt 50 and science[i].t lt 50 then science[i].phot_color = 'MT'
        endfor
    endif else begin
        for i=0,nscience-1 do begin
            case 1 of
                science[i].v gt 0 and science[i].k gt 0 and science[i].v lt 50 and science[i].k lt 50: science[i].phot_color = 'VK'
                science[i].v gt 0 and science[i].h gt 0 and science[i].v lt 50 and science[i].h lt 50: science[i].phot_color = 'VH'
                science[i].v gt 0 and science[i].j gt 0 and science[i].v lt 50 and science[i].j lt 50: science[i].phot_color = 'VJ'
                science[i].v gt 0 and science[i].i gt 0 and science[i].v lt 50 and science[i].i lt 50: science[i].phot_color = 'VI'
                science[i].b gt 0 and science[i].r gt 0 and science[i].b lt 50 and science[i].r lt 50: science[i].phot_color = 'BR'
                science[i].v gt 0 and science[i].r gt 0 and science[i].v lt 50 and science[i].r lt 50: science[i].phot_color = 'VR'
                science[i].b gt 0 and science[i].v gt 0 and science[i].b lt 50 and science[i].v lt 50: science[i].phot_color = 'BV'
                science[i].r gt 0 and science[i].i gt 0 and science[i].r lt 50 and science[i].i lt 50: science[i].phot_color = 'RI'
                science[i].j gt 0 and science[i].k gt 0 and science[i].j lt 50 and science[i].k lt 50: science[i].phot_color = 'JK'
                else:
            endcase
        endfor
    endelse

    w = where(science.phot_color eq 'VK', c)
    if c gt 0 then begin
        wnoerror = where(science[w].verr le 0 or science[w].verr ge 10 or science[w].kerr le 0 or science[w].kerr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].k, science[w[wnoerror]].v-science[w[wnoerror]].k, /vk, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].k, science[w[wnoerror]].v-science[w[wnoerror]].k, /vk, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].k, science[w[wwitherror]].v-science[w[wwitherror]].k, science[w[wwitherror]].kerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].kerr)^2.), /vk, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].k, science[w[wwitherror]].v-science[w[wwitherror]].k, science[w[wwitherror]].kerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].kerr)^2.), /vk, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        science.teff_alonso = -999.0
        science.tefferr_alonso = -999.0

        color[w] = science[w].v - science[w].k
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.4942, 0.2809, -0.0180, -0.0294, 0.0444, -0.0008]
            aramgiant[w[i],*] = [0.4405, 0.3272, -0.0252, -0.0016, -0.0053, -0.0040]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-1425.36, 3218.36, -2566.54, 859.644, -102.554, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.896, 3.360]
                    pramgiant[w[i],*] = [-72.6664, 36.5361, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.244, 3.286]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [2.35133, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [1.060, 2.665]
                    pramgiant[w[i],*] = [86.0358, -65.4928, 10.8901, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.366, 4.474]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [-1849.46, 4577.00, -4284.02, 1770.38, -268.589, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [1.101, 2.670]
                    pramgiant[w[i],*] = [-6.96153, 14.3298, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.334, 3.549]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [215.721, -796.519, 714.423, -175.678, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [1.126, 2.596]
                    pramgiant[w[i],*] = [-943.925, 1497.64, -795.867, 138.965, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.258, 2.768]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            if color[w[i]] lt 2.25 then begin
                aalonso[w[i],*] = [0.5558, 0.2105, 1.981e-3, -9.965e-3, 1.325e-2, -2.726e-3]
                syserralonso[w[i],*] = 40d
                case 1 of
                    science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.20, 2.5]
                    science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [1.00, 2.5]
                    science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [1.20, 2.5]
                    science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [1.70, 2.5]
                    else: clrrangealonso[w[i],*] = [-900, -900]
                endcase
            endif else begin
                aalonso[w[i],*] = [0.3770, 0.3660, -3.170e-2, -3.074e-3, -2.765e-3, -2.973e-3]
                syserralonso[w[i],*] = 25d
                case 1 of
                    science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [2.00, 4.90]
                    science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [2.00, 4.60]
                    science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [2.00, 3.40]
                    science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [2.00, 2.80]
                    else: clrrangealonso[w[i],*] = [-900, -900]
                endcase
            endelse
        endfor
        syserrram[w] = 28d

        teffgoodram[w] = 1
        teffgoodalonso[w] = 1
    endif

    w = where(science.phot_color eq 'VH', c)
    if c gt 0 then begin
        wnoerror = where(science[w].verr le 0 or science[w].verr ge 10 or science[w].herr le 0 or science[w].herr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].h, science[w[wnoerror]].v-science[w[wnoerror]].h, /vh, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].h, science[w[wnoerror]].v-science[w[wnoerror]].h, /vh, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].h, science[w[wwitherror]].v-science[w[wwitherror]].h, science[w[wwitherror]].herr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].herr)^2.), /vh, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].h, science[w[wwitherror]].v-science[w[wwitherror]].h, science[w[wwitherror]].herr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].herr)^2.), /vh, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].v - science[w].h
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.4931, 0.3056, -0.0241, -0.0396, 0.0678, 0.0020]
            aramgiant[w[i],*] = [0.4354, 0.3405, -0.0263, -0.0012, -0.0049, -0.0027]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-53.5574, 36.0990, 15.6878, -8.84468]
                    clrrangeramdwarf[w[i],*] = [0.839, 3.215]
                    pramgiant[w[i],*] = [-377.022, 334.733, -69.8093, 0.0000, 0.0000]
                    clrrangeramgiant[w[i],*] = [1.194, 3.059]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [1.60629]
                    clrrangeramdwarf[w[i],*] = [1.032, 2.532]
                    pramgiant[w[i],*] = [71.7949, -55.5383, 9.61821, 0.0000, 0.0000]
                    clrrangeramgiant[w[i],*] = [1.293, 4.263]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [506.559, -1277.52, 939.519, -208.624]
                    clrrangeramdwarf[w[i],*] = [1.070, 2.535]
                    pramgiant[w[i],*] = [-27.4190, 20.7082, 0.0000, 0.0000, 0.0000]
                    clrrangeramgiant[w[i],*] = [1.273, 3.416]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [-471.588, 643.972, -199.639]
                    clrrangeramdwarf[w[i],*] = [1.093, 2.388]
                    pramgiant[w[i],*] = [-46.2946, 20.1061, 0.0000, 0.0000, 0.0000]
                    clrrangeramgiant[w[i],*] = [1.232, 2.625]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase
        endfor
        syserrram[w] = 32d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 0
    endif

    w = where(science.phot_color eq 'VJ', c)
    if c gt 0 then begin
        wnoerror = where(science[w].verr le 0 or science[w].verr ge 10 or science[w].jerr le 0 or science[w].jerr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].j, science[w[wnoerror]].v-science[w[wnoerror]].j, /vj, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].j, science[w[wnoerror]].v-science[w[wnoerror]].j, /vj, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].j, science[w[wwitherror]].v-science[w[wwitherror]].j, science[w[wwitherror]].jerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].jerr)^2.), /vj, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].j, science[w[wwitherror]].v-science[w[wwitherror]].j, science[w[wwitherror]].jerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].jerr)^2.), /vj, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].v - science[w].j
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.4050, 0.4792, -0.0617, -0.0392, 0.0401, -0.0023]
            aramgiant[w[i],*] = [0.2943, 0.5604, -0.0677, 0.0179, -0.0532, -0.0088]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [422.406, -910.603, 621.335, -132.566, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.815, 2.608]
                    pramgiant[w[i],*] = [-122.595, 76.4847, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.259, 2.400]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [-466.616, 658.349, -220.454, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.860, 2.087]
                    pramgiant[w[i],*] = [-10.3848, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.030, 3.418]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [-862.072, 1236.84, -423.729, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.927, 1.983]
                    pramgiant[w[i],*] = [4.18695, 13.8937, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [1.033, 2.679]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [-1046.10, 1652.06, -597.340, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.891, 1.932]
                    pramgiant[w[i],*] = [-67.7716, 28.9202, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.977, 2.048]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase
        endfor
        syserrram[w] = 38d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 0
    endif

    w = where(science.phot_color eq 'MT', c)
    if c gt 0 then begin
        wnoerror = where(science[w].merr le 0 or science[w].merr ge 10 or science[w].terr le 0 or science[w].terr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        i_mt = science.t
        vi_mt = 0.8*(science.m-science.t) - 0.006
        ierr_mt = science.terr
        vierr_mt = sqrt((0.8*science.merr)^2. + (0.8*science.terr)^2. + (0.014)^2.)

        if cnoerror gt 0 then begin
            fehphot = feh(i_mt[w[wnoerror]], vi_mt[w[wnoerror]], /vi, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(i_mt[w[wnoerror]], vi_mt[w[wnoerror]], /vi, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(i_mt[w[wnoerror]], vi_mt[w[wnoerror]], /vi, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(i_mt[w[wwitherror]], vi_mt[w[wwitherror]], ierr_mt[w[wwitherror]], vierr_mt[w[wwitherror]], /vi, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(i_mt[w[wwitherror]], vi_mt[w[wwitherror]], ierr_mt[w[wwitherror]], vierr_mt[w[wwitherror]], /vi, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(i_mt[w[wwitherror]], vi_mt[w[wwitherror]], ierr_mt[w[wwitherror]], vierr_mt[w[wwitherror]], /vi, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = vi_mt[w]
        colorerr[w] = vierr_mt[w]
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.3295, 0.9516, -0.2290, -0.0316, 0.0003, -0.0081]
            aramgiant[w[i],*] = [0.3575, 0.9069, -0.2025, 0.0395, -0.0551, -0.0061]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-2757.79, 9961.33, -10546.6, -1746.05, 10512.3, -6653.57, 1301.21]
                    clrrangeramdwarf[w[i],*] = [0.491, 1.721]
                    pramgiant[w[i],*] = [0.42933, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.573, 2.000]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [-22.9008, 40.2078, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.597, 1.052]
                    pramgiant[w[i],*] = [-0.14180, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.795, 1.524]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [-667.732, 1709.88, -1069.62, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.547, 1.026]
                    pramgiant[w[i],*] = [9.31011, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.870, 1.303]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.000, 0.000]
                    pramgiant[w[i],*] = [-23.0514, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.812, 1.095]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            aalonso[w[i],*] = [0.5379, 0.3981, 4.432e-2, -2.693e-2, 0.0000, 0.0000]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.20, 2.90]
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.80, 2.00]
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.85, 2.20]
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [1.00, 1.70]
                else: clrrangealonso[w[i],*] = [-900, -900]
            endcase
        endfor
        syserrram[w] = 40d
        syserralonso[w] = 125d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 2
    endif

    w = where(science.phot_color eq 'VI', c)
    if c gt 0 then begin
        wnoerror = where(science[w].verr le 0 or science[w].verr ge 10 or science[w].ierr le 0 or science[w].ierr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)            

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].v-science[w[wnoerror]].i, /vi, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].v-science[w[wnoerror]].i, /vi, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].v-science[w[wnoerror]].i, /vi, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].v-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].ierr)^2.), /vi, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].v-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].ierr)^2.), /vi, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].v-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].ierr)^2.), /vi, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].v - science[w].i
        colorerr[w] = sqrt((science[w].verr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.3295, 0.9516, -0.2290, -0.0316, 0.0003, -0.0081]
            aramgiant[w[i],*] = [0.3575, 0.9069, -0.2025, 0.0395, -0.0551, -0.0061]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-2757.79, 9961.33, -10546.6, -1746.05, 10512.3, -6653.57, 1301.21]
                    clrrangeramdwarf[w[i],*] = [0.491, 1.721]
                    pramgiant[w[i],*] = [0.42933, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.573, 2.000]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [-22.9008, 40.2078, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.597, 1.052]
                    pramgiant[w[i],*] = [-0.14180, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.795, 1.524]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [-667.732, 1709.88, -1069.62, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.547, 1.026]
                    pramgiant[w[i],*] = [9.31011, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.870, 1.303]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.000, 0.000]
                    pramgiant[w[i],*] = [-23.0514, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.812, 1.095]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            aalonso[w[i],*] = [0.5379, 0.3981, 4.432e-2, -2.693e-2, 0.0000, 0.0000]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.20, 2.90]
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.80, 2.00]
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.85, 2.20]
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [1.00, 1.70]
                else: clrrangealonso[w[i],*] = [-900, -900]
            endcase
        endfor
        syserrram[w] = 40d
        syserralonso[w] = 125d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 2
    endif

    w = where(science.phot_color eq 'VR', c)
    if c gt 0 then begin
        wnoerror = where(science[w].verr le 0 or science[w].verr ge 10 or science[w].rerr le 0 or science[w].rerr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].v-science[w[wnoerror]].r, /vr, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].v-science[w[wnoerror]].r, /vr, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].v-science[w[wnoerror]].r, /vr, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].v-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].rerr)^2.), /vr, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].v-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].rerr)^2.), /vr, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].v-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].verr)^2.+(science[w[wwitherror]].rerr)^2.), /vr, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].v - science[w].r
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.4333, 1.4399, -0.5419, -0.0481, -0.0239, -0.0125]
            aramgiant[w[i],*] = [0.3849, 1.6205, -0.6395, 0.1060, -0.0875, -0.0089]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-2666.55, 27264.5, -103923, 174663, -104940, -23249.4, 32644.9]
                    clrrangeramdwarf[w[i],*] = [0.204, 0.880]
                    pramgiant[w[i],*] = [-8.51797, 15.6675, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.299, 1.106]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [4.20153, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.284, 0.546]
                    pramgiant[w[i],*] = [-10.7764, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.387, 0.752]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [123.940, -342.217, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.264, 0.532]
                    pramgiant[w[i],*] = [61.9821, -78.7382, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.429, 0.598]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [8.55498, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.240, 0.336]
                    pramgiant[w[i],*] = [27.9886, -100.149, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.394, 0.550]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            aalonso[w[i],*] = [0.4972, 0.8841, -0.1904, -1.197e-2, -1.025e-2, -5.500e-3]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.15, 1.70]
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.45, 1.50]
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.50, 1.00]
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [0.55, 0.85]
                else: clrrangealonso[w[i],*] = [-900, -900]
            endcase
        endfor
        syserrram[w] = 41d
        syserralonso[w] = 150d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 1
    endif

    w = where(science.phot_color eq 'BV', c)
    if c gt 0 then begin
        wnoerror = where(science[w].berr le 0 or science[w].berr ge 10 or science[w].verr le 0 or science[w].verr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].v, science[w[wnoerror]].b-science[w[wnoerror]].v, /bv, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(science[w[wnoerror]].v, science[w[wnoerror]].b-science[w[wnoerror]].v, /bv, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].v, science[w[wnoerror]].b-science[w[wnoerror]].v, /bv, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].v, science[w[wwitherror]].b-science[w[wwitherror]].v, science[w[wwitherror]].verr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].verr)^2.), /bv, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].v, science[w[wwitherror]].b-science[w[wwitherror]].v, science[w[wwitherror]].verr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].verr)^2.), /bv, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].v, science[w[wwitherror]].b-science[w[wwitherror]].v, science[w[wwitherror]].verr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].verr)^2.), /bv, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].b - science[w].v
        colorerr[w] = sqrt((science[w].berr)^2. + (science[w].verr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.5002, 0.6440, -0.0690, -0.0230, -0.0566, -0.0170]
            aramgiant[w[i],*] = [0.5737, 0.4882, -0.0149, 0.0563, -0.1160, -0.0114]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-261.548, 684.977, -470.049, 79.8977, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.310, 1.507]
                    pramgiant[w[i],*] = [112.116, -372.622, 67.1254, 395.333, -203.471, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.144, 1.668]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [-324.033, 1516.44, -2107.37, 852.150, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.307, 1.202]
                    pramgiant[w[i],*] = [-12.9762, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.664, 1.558]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [30.5985, -46.7882, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.335, 1.030]
                    pramgiant[w[i],*] = [606.032, -1248.79, 627.453, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.605, 1.352]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [139.965, -292.329, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.343, 0.976]
                    pramgiant[w[i],*] = [-9.26209, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.680, 1.110]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            if color[w[i]] lt 0.75 then begin
                aalonso[w[i],*] = [0.5716, 0.5404, -6.126e-2, -4.862e-2, -1.777e-2, -7.969e-3]
                syserralonso[w[i],*] = 167d
                case 1 of
                    science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.20, 0.80]
                    science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.35, 0.80]
                    science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.35, 0.80]
                    science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [0.50, 0.80]
                    else: clrrangealonso[w[i],*] = [-900, -900]
                endcase
            endif else begin
                aalonso[w[i],*] = [0.6177, 0.4354, -4.025e-3, 5.204e-2, -0.1127, -1.385e-2]
                syserralonso[w[i],*] = 96d
                case 1 of
                    science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.70, 1.90]
                    science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.70, 1.80]
                    science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.70, 1.35]
                    science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [0.70, 1.00]
                    else: clrrangealonso[w[i],*] = [-900, -900]
                endcase
            endelse
        endfor
        syserrram[w] = 51d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 1
    endif

    w = where(science.phot_color eq 'BR', c)
    if c gt 0 then begin
        wnoerror = where(science[w].berr le 0 or science[w].berr ge 10 or science[w].rerr le 0 or science[w].rerr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].b-science[w[wnoerror]].r, /br, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].b-science[w[wnoerror]].r, /br, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].r, science[w[wnoerror]].b-science[w[wnoerror]].r, /br, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].b-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].rerr)^2.), /br, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].b-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].rerr)^2.), /br, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].r, science[w[wwitherror]].b-science[w[wwitherror]].r, science[w[wwitherror]].rerr, sqrt((science[w[wwitherror]].berr)^2.+(science[w[wwitherror]].rerr)^2.), /br, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        teffgoodram[w] = 0
        teffgoodalonso[w] = 0
    endif

    w = where(science.phot_color eq 'RI', c)
    if c gt 0 then begin
        wnoerror = where(science[w].rerr le 0 or science[w].rerr ge 10 or science[w].ierr le 0 or science[w].ierr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].r-science[w[wnoerror]].i, /ri, /van, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_van = fehphot.feh
            science[w[wnoerror]].teff_van = fehphot.teff
            science[w[wnoerror]].logg_van = fehphot.logg
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].r-science[w[wnoerror]].i, /ri, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].i, science[w[wnoerror]].r-science[w[wnoerror]].i, /ri, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].r-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].rerr)^2.+(science[w[wwitherror]].ierr)^2.), /ri, /van, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_van = fehphot.feh
            science[w[wwitherror]].teff_van = fehphot.teff
            science[w[wwitherror]].logg_van = fehphot.logg
            science[w[wwitherror]].feherr_van = fehphot.err_feh
            science[w[wwitherror]].tefferr_van = fehphot.err_teff
            science[w[wwitherror]].loggerr_van = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].r-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].rerr)^2.+(science[w[wwitherror]].ierr)^2.), /ri, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].i, science[w[wwitherror]].r-science[w[wwitherror]].i, science[w[wwitherror]].ierr, sqrt((science[w[wwitherror]].rerr)^2.+(science[w[wwitherror]].ierr)^2.), /ri, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].r - science[w].i
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aramdwarf[w[i],*] = [0.2919, 2.1141, -1.0723, -0.0756, 0.0267, -0.0041]
            aramgiant[w[i],*] = [0.4351, 1.6549, -0.7215, -0.0610, 0.0332, -0.0023]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] lt 1.0 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: begin
                    pramdwarf[w[i],*] = [-3326.97, 26263.8, -75355.8, 94246.5, -43334.8, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.242, 0.838]
                    pramgiant[w[i],*] = [61.3557, -116.711, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.413, 0.793]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: begin
                    pramdwarf[w[i],*] = [12.4740, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.300, 0.718]
                    pramgiant[w[i],*] = [-16.8645, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.383, 0.771]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: begin
                    pramdwarf[w[i],*] = [-5837.31, 41439.2, -94729.8, 69584.8, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.283, 0.551]
                    pramgiant[w[i],*] = [32.0870, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.434, 0.725]
                end
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5 and science[w[i]].feh_yy[wyy[w[i]]] ge -4.0: begin
                    pramdwarf[w[i],*] = [32.1826, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
                    clrrangeramdwarf[w[i],*] = [0.290, 0.364]
                    pramgiant[w[i],*] = [-15.6318, 0.0000, 0.0000, 0.0000, 0.0000, 0.000, 0.000]
                    clrrangeramgiant[w[i],*] = [0.364, 0.545]
                end
                else: begin
                    clrrangeramdwarf[w[i],*] = [-900, -900]
                    clrrangeramgiant[w[i],*] = [-900, -900]
                end
            endcase

            aalonso[w[i],*] = [0.4974, 1.345, -0.5008, -8.134e-2, 3.705e-2, -6.184e-3]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.15, 1.40]
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.25, 0.80]
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.35, 0.70]
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [0.40, 0.65]
                else: clrrangealonso[w[i],*] = [-900, -900]
            endcase
        endfor
        syserrram[w] = 62d
        syserralonso[w] = 150d
        teffgoodram[w] = 1
        teffgoodalonso[w] = 1
    endif

    w = where(science.phot_color eq 'JK', c)
    if c gt 0 then begin
        wnoerror = where(science[w].jerr le 0 or science[w].jerr ge 10 or science[w].kerr le 0 or science[w].kerr ge 10, cnoerror, complement=wwitherror, ncomplement=cwitherror)

        if cnoerror gt 0 then begin
            fehphot = feh(science[w[wnoerror]].k, science[w[wnoerror]].j-science[w[wnoerror]].k, /jk, /padova, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_padova = fehphot.feh
            science[w[wnoerror]].teff_padova = fehphot.teff
            science[w[wnoerror]].logg_padova = fehphot.logg
            fehphot = feh(science[w[wnoerror]].k, science[w[wnoerror]].j-science[w[wnoerror]].k, /jk, /yy, alpha=alpha, dm=dm[w[wnoerror]] - 3.1*ebv[w[wnoerror]], ddm=dmerr[w[wnoerror]])
            science[w[wnoerror]].feh_yy = fehphot.feh
            science[w[wnoerror]].teff_yy = fehphot.teff
            science[w[wnoerror]].logg_yy = fehphot.logg
        endif

        if cwitherror gt 0 then begin
            fehphot = feh(science[w[wwitherror]].k, science[w[wwitherror]].j-science[w[wwitherror]].k, science[w[wwitherror]].kerr, sqrt((science[w[wwitherror]].jerr)^2.+(science[w[wwitherror]].kerr)^2.), /jk, /padova, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_padova = fehphot.feh
            science[w[wwitherror]].teff_padova = fehphot.teff
            science[w[wwitherror]].logg_padova = fehphot.logg
            science[w[wwitherror]].feherr_padova = fehphot.err_feh
            science[w[wwitherror]].tefferr_padova = fehphot.err_teff
            science[w[wwitherror]].loggerr_padova = fehphot.err_logg
            fehphot = feh(science[w[wwitherror]].k, science[w[wwitherror]].j-science[w[wwitherror]].k, science[w[wwitherror]].kerr, sqrt((science[w[wwitherror]].jerr)^2.+(science[w[wwitherror]].kerr)^2.), /jk, /yy, alpha=alpha, dm=dm[w[wwitherror]] - 3.1*ebv[w[wwitherror]], ddm=dmerr[w[wwitherror]])
            science[w[wwitherror]].feh_yy = fehphot.feh
            science[w[wwitherror]].teff_yy = fehphot.teff
            science[w[wwitherror]].logg_yy = fehphot.logg
            science[w[wwitherror]].feherr_yy = fehphot.err_feh
            science[w[wwitherror]].tefferr_yy = fehphot.err_teff
            science[w[wwitherror]].loggerr_yy = fehphot.err_logg
        endif

        color[w] = science[w].r - science[w].i
        colorerr[w] = sqrt((science[w].rerr)^2. + (science[w].ierr)^2.)
        for i=0,c-1 do begin
            aalonso[w[i],*] = [0.5816, 0.9134, -0.1443, 0.0000, 0.0000, 0.0000]
            case 1 of
                science[w[i]].feh_yy[wyy[w[i]]] le 0.2 and science[w[i]].feh_yy[wyy[w[i]]] gt -0.5: clrrangealonso[w[i],*] = [0.00, 1.10]
                science[w[i]].feh_yy[wyy[w[i]]] le -0.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -1.5: clrrangealonso[w[i],*] = [0.20, 1.00]
                science[w[i]].feh_yy[wyy[w[i]]] le -1.5 and science[w[i]].feh_yy[wyy[w[i]]] gt -2.5: clrrangealonso[w[i],*] = [0.30, 0.90]
                science[w[i]].feh_yy[wyy[w[i]]] le -2.5: clrrangealonso[w[i],*] = [0.40, 0.80]
                else: clrrangealonso[w[i],*] = [-900, -900]
            endcase
        endfor
        syserralonso[w] = 125d
        teffgoodram[w] = 0
        teffgoodalonso[w] = 1
    endif

    w = where(teffgoodram eq 1, c)
    if c gt 0 then begin
        for i=0,c-1 do begin
            case 1 of
                science[w[i]].logg_yy[wyy[w[i]]] gt 0.0 and science[w[i]].logg_yy[wyy[w[i]]] lt 3.55: begin
                    aram = reform(aramgiant[w[i],*])
                    p = reform(pramgiant[w[i],*])
                    clrrangeram = reform(clrrangeramgiant[w[i],*])
                    goodlogg = 1
                end
                science[w[i]].logg_yy[wyy[w[i]]] ge 3.55 and science[w[i]].logg_yy[wyy[w[i]]] lt 5.0: begin
                    aram = reform(aramdwarf[w[i],*])
                    p = reform(pramdwarf[w[i],*])
                    clrrangeram = reform(clrrangeramdwarf[w[i],*])
                    goodlogg = 2
                end
                else: begin
                    goodlogg = 0
                end
            endcase

            if goodlogg eq 0 then continue
            if color[w[i]] lt clrrangeram[0] or color[w[i]] gt clrrangeram[1] then continue
            
            theta = aram[0] + aram[1]*color[w[i]] + aram[2]*(color[w[i]])^2. + aram[3]*color[w[i]]*science[w[i]].feh_yy[wyy[w[i]]] + aram[4]*science[w[i]].feh_yy[wyy[w[i]]] + aram[5]*(science[w[i]].feh_yy[wyy[w[i]]])^2.
            thetaerrclr = (aram[1] + 2.*aram[2]*color[w[i]] + aram[3]*science[w[i]].feh_yy[wyy[w[i]]])*colorerr[w[i]]
            thetaerrfeh = (aram[3]*color[w[i]] + aram[4] + 2.*aram[5]*science[w[i]].feh_yy[wyy[w[i]]])*science[w[i]].feherr_yy[wyy[w[i]]]
            teffpoly = poly(color[w[i]], p[*])
            teffpolyerr = poly(color[w[i]], (dindgen(6)+1.0)*p[1:6])*colorerr[w[i]]
            science[w[i]].teff_ram = (5040./theta) + teffpoly
            science[w[i]].tefferr_ram = sqrt(((5040./theta^2.)^2.)*(thetaerrclr^2. + thetaerrfeh^2.) + teffpolyerr^2. + syserrram[w[i]]^2.)
        endfor
    endif

    w = where(teffgoodalonso eq 1, c)
    if c gt 0 then begin
        for i=0,c-1 do begin
            ;if color[w[i]] lt clrrangealonso[w[i],0] or color[w[i]] gt clrrangealonso[w[i],1] then continue
            theta = aalonso[w[i],0] + aalonso[w[i],1]*color[w[i]] + aalonso[w[i],2]*(color[w[i]])^2. + aalonso[w[i],3]*color[w[i]]*science[w[i]].feh_yy[wyy[w[i]]] + aalonso[w[i],4]*science[w[i]].feh_yy[wyy[w[i]]] + aalonso[w[i],5]*(science[w[i]].feh_yy[wyy[w[i]]])^2.
            thetaerrclr = (aalonso[w[i],1] + 2.*aalonso[w[i],2]*color[w[i]] + aalonso[w[i],3]*science[w[i]].feh_yy[wyy[w[i]]])*colorerr[w[i]]
            thetaerrfeh = (aalonso[w[i],3]*color[w[i]] + aalonso[w[i],4] + 2.*aalonso[w[i],5]*science[w[i]].feh_yy[wyy[w[i]]])*science[w[i]].feherr_yy[wyy[w[i]]]
            science[w[i]].teff_alonso = 5040./theta
            science[w[i]].tefferr_alonso = sqrt(((5040./theta^2.)^2.)*(thetaerrclr^2. + thetaerrfeh^2.) + syserralonso[w[i]]^2.)
        endfor
    endif

    w = where(teffgoodalonso eq 2, c)
    if c gt 0 then begin
        for i=0,c-1 do begin
            ;if color[w[i]] lt clrrangealonso[w[i],0] or color[w[i]] gt clrrangealonso[w[i],1] then continue
            theta = aalonso[w[i],0] + aalonso[w[i],1]*color[w[i]] + aalonso[w[i],2]*(color[w[i]])^2. + aalonso[w[i],3]*(color[w[i]])^3.
            thetaerrclr = (aalonso[w[i],1] + 2.*aalonso[w[i],2]*color[w[i]] + 3.*aalonso[w[i],3]*(color[w[i]])^2.)*colorerr[w[i]]
            science[w[i]].teff_alonso = 5040./theta
            science[w[i]].tefferr_alonso = sqrt(((5040./theta^2.)^2.)*(thetaerrclr^2.) + syserralonso[w[i]]^2.)
        endfor
    endif
    return, science
end
