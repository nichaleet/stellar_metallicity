function moog_li, teff, logg, feh, alphafe, vt, ali=ali
    retry:
    semname = string(floor((systime(1) mod 10) * 1e8), format='(I09)')
    status = sem_create(semname)
    semlocked = sem_lock(semname)
    if ~semlocked then begin
        sem_delete, semname
        goto, retry
    endif
    basename = getenv('CALTECH')+'lithium/moog/'+semname
    parfile = basename+'.par'
    if file_test(parfile) then begin
        sem_release, semname
        sem_delete, semname
        goto, retry
    endif
    linefile = getenv('CALTECH')+'lithium/LiI6707.moog'
    atmfile = basename+'.atm'
    out1file = basename+'.out1'
    out2file = basename+'.out2'

    interp_atm, teff, logg, vt, feh, alphafe<1.1999, outfile=atmfile, tweakel=3, tweakabund=ali
    make_par, parfile=parfile, linefile=linefile, atmfile=atmfile, outfile=out2file, driver='synth', minlambda=6697.0, maxlambda=6717.0, extralines=['isotopes     2    1', '    3.006    29.000', '    3.007     1.036']
    spawn, 'MOOGSILENT '+parfile
    moog = read_moog_spec(out2file, /newmoog)
    file_delete, atmfile
    file_delete, parfile
    file_delete, out1file
    file_delete, out2file

    sem_release, semname
    sem_delete, semname

    return, moog
end

function get_moog_li, x, p
    common resolution, dlambda
    teff = p[0]
    logg = p[1]
    feh = p[2]
    alphafe = p[3]
    vt = p[4]
    ali = p[5]

    moogspec = moog_li(teff, logg, feh, alphafe, vt, ali=ali)

    w = where(x ge 6697 and x le 6717)
    moogspecsmooth = smooth_gauss_wrapper(moogspec.lambda, moogspec.spec, x[w], dlambda[w])

    return, moogspecsmooth
end
