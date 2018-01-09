function read_moog_baseline_ch, status=status, fullres=fullres, mds=mds
    path = getenv('chome')+'gridch'+(keyword_set(mds) ? '_mds/' : '/')+(keyword_set(fullres) ? 'synths/' : 'bin/')
    n = keyword_set(fullres) ? 20001 : 2858
    bingzfile = path+'baseline.bin.gz'

    status = -1
    if ~file_test(bingzfile) then begin
        print, bingzfile
        message, 'MOOG binary spectrum not found.';, /info
        status = 0
        return, [-1]
    endif

    moogspec_ch = fltarr(n)
    openr, 1, bingzfile, /compress
    readu, 1, moogspec_ch
    close, 1

    status = 1
    return, 1.0 - moogspec_ch
end
