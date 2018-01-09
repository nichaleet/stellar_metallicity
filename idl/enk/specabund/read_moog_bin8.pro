function read_moog_bin8, teff, logg, feh, alphafe, mgsicati, n=n, status=status, fullres=fullres
    path7 = getenv('chome')+'grid7/'+(keyword_set(fullres) ? 'synths/' : 'bin/')
    path8 = getenv('chome')+'grid8/'+(keyword_set(fullres) ? 'synths/' : 'bin/')

    if (size(mgsicati))[1] eq 0 then mgsicati = alphafe
    teff = round(teff/100.)*100L
    logg = round(logg*10.)
    feh = round(feh*10.)
    alphafe = round(alphafe*10.)
    mgsicati = round(mgsicati*10.)

    gsign = logg ge 0 ? '_' : '-'
    fsign = feh ge 0 ? '_' : '-'
    asign = alphafe ge 0 ? '_' : '-'
    msign = mgsicati ge 0 ? '_' : '-'
    bingzfile = path7+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,".bin.gz")')
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

    if alphafe ne mgsicati then begin
        bingzfile = path8+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), msign, abs(mgsicati), format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,"a",A1,I02,".bin.gz")')
        if ~file_test(bingzfile) then begin
            print, bingzfile
            message, 'MOOG binary alpha-modified spectrum not found.';, /info
            status = 0
            return, [-1]
        endif

        moogspecalpha = fltarr(keyword_set(fullres) ? 36442 : 5206)
        openr, 1, bingzfile, /compress
        readu, 1, moogspecalpha
        close, 1

        restore, getenv('M31')+'specabund/moogify/specregion/alphamask_bin.sav'
        moogspec[wbin] = moogspecalpha
    endif

    status = 1
    return, 1.0 - moogspec
end
