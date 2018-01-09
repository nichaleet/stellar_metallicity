function read_moog_bin_ch, teff, logg, feh, alphafe, cfe, status=status, fullres=fullres, mds=mds
    path = getenv('chome')+'gridch'+(keyword_set(mds) ? '_mds/' : '/')+(keyword_set(fullres) ? 'synths/' : 'bin/')
    n = keyword_set(fullres) ? 20001L : 2858L
    ;nall = n * 24
    acfe = double([-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2,  0.0,  0.2, 0.4,  0.6,  0.8,  1.0,  1.4, 1.8,  2.2,  2.6,  3.0,  3.5])
    icfe = where(acfe eq cfe, c)
    if c ne 1 then message, '[C/Fe] out of range.'

    teff = round(teff/100.)*100L
    logg = round(logg*10.)
    feh = round(feh*10.)
    alphafe = round(alphafe*10.)
    cfe = round(cfe*10.)

    gsign = logg ge 0 ? '_' : '-'
    fsign = feh ge 0 ? '_' : '-'
    asign = alphafe ge 0 ? '_' : '-'
    ;csign = cfe ge 0 ? '_' : '-'
    ;bingzfile = path+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), csign, abs(cfe), format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,"c",A1,I02,".bin.gz")')
    bingzfile = path+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,".bin.gz")')
    ;print, bingzfile

    status = -1
    if ~file_test(bingzfile) then begin
        print, bingzfile
        message, 'MOOG binary spectrum not found.';, /info
        status = 0
        return, [-1]
    endif

    moogspec_ch = fltarr(n)
    openr, 1, bingzfile, /compress
    if icfe gt 0 then begin
        junk = fltarr(n*icfe)
        readu, 1, junk
    endif
    readu, 1, moogspec_ch
    close, 1

    ;common moogch_baseline, baseline
    ;if n_elements(baseline) le 0 then baseline = read_moog_baseline_ch(fullres=fullres, mds=mds)
    moogspec_ch = 1.0 - moogspec_ch
    ;moogspec_ch /= baseline

    status = 1
    return, moogspec_ch
end
