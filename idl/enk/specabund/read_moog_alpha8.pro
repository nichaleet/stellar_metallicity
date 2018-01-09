function read_moog_alpha8, teff, logg, feh, alphafe, mgsicati, n=n, status=status, fullres=fullres
    path = getenv('chome')+'grid8/'+(keyword_set(fullres) ? 'synths/' : 'bin/')

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

    n = keyword_set(fullres) ? 36442 : 5206
    moogspecalpha = fltarr(n)

    if alphafe ne mgsicati then begin
        bingzfile = path+string(teff, gsign, abs(logg), teff, gsign, abs(logg), fsign, abs(feh), asign, abs(alphafe), msign, abs(mgsicati), format='("t",I04,"/g",A1,I02,"/t",I04,"g",A1,I02,"f",A1,I02,"a",A1,I02,"a",A1,I02,".bin.gz")')
        if ~file_test(bingzfile) then begin
            status = 0
            print, bingzfile
            message, 'MOOG binary alpha-modified spectrum not found.';, /info
            return, [-1]
        endif

        openr, 1, bingzfile, /compress
        readu, 1, moogspecalpha
        close, 1
    endif else message, 'There are no grid8 files where alphafe = mgsicati.'

    status = 1
    return, 1.0 - moogspecalpha
end
