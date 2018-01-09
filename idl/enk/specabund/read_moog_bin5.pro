function read_moog_bin5, teff, logg, feh, alphafe, vt, n=n, status=status, fullres=fullres
    path = '/net/puccini/c/ekirby/grid5/'+(keyword_set(fullres) ? 'synths/' : 'bin/')

    teff = round(teff/100.)*100L
    logg = round(logg*10.)
    feh = round(feh*10.)
    alphafe = round(alphafe*10.)
    vt = round(vt*100.)

    gsign = logg ge 0 ? '_' : '-'
    fsign = feh ge 0 ? '_' : '-'
    asign = alphafe ge 0 ? '_' : '-'
    bingzfile = path+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), vt, format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,"v",I03,".bin.gz")')
    ;print, bingzfile

    status = -1
    if ~file_test(bingzfile) then begin
        print, bingzfile
        message, 'MOOG binary spectrum not found.';, /info
        status = 0
        return, [-1]
    endif

    moogspec = fltarr(n)
    openr, 1, bingzfile, /compress
    readu, 1, moogspec
    close, 1

    status = 1
    return, 1.0 - moogspec
end
