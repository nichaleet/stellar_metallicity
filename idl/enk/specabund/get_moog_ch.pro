function get_moog_ch, xin, a, mds=mds
    common moog_ch, datalambda, dlam, moogmask, contmask, contdiv, contdivivar, smoothfactor, phot, moogcontrevise, pixmask
    common mooglambdacom_ch, mooglambda_ch
    common chi2common, xcont, ycont, dycont

    teff = a[0]
    logg = a[1]
    feh = a[2]
    alpha = a[3]
    cfe = a[4]
    resscale = a[5]

    if phot eq 1 then x = xin[2:n_elements(xin)-1] else x = xin

    moogspec_ch = interp_moog_ch(teff, logg, feh, alpha, cfe, mds=mds)
    moogspec_ch = smooth_gauss_wrapper(mooglambda_ch, moogspec_ch, xcont, dlam*resscale)

    match, xcont, x, wmask, wjunk
    moogspec_ch = moogspec_ch[wmask]
    return, phot eq 1 ? [teff, logg, moogspec_ch] : moogspec_ch
end
