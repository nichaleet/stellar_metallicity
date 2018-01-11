pro test_coadd_improved
    moogify = mrdfits(getenv('M31')+'data/gc/n7078/moogify.fits', 1, /silent)
    pars = mrdfits(getenv('M31')+'data/gc/n7078/moogify7_flexteff.fits', 1, /silent)
    match, strtrim(moogify.objname, 2), strtrim(pars.objname, 2), w1, w2
    moogify = moogify[w1]
    pars = pars[w2]

    w = where(pars.teff gt 4400 and pars.teff lt 4600)
    moogify = moogify[w]
    coadd = coadd_improved(moogify[w], /inorm)

    splot, coadd.lambda, coadd.spec
end
