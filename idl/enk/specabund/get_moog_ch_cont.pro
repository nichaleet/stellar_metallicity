function get_moog_ch_cont, xin, a, mds=mds
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

    won = where(moogspec_ch gt 0 and finite(dycont) and dycont gt 0.0, complement=woff, con)
    if con lt 100 then speccont = dblarr(n_elements(xcont))+1.0 else begin 
        bkpt = slatec_splinefit(xcont[won], ycont[won]/moogspec_ch[won], coeff, invvar=dycont[won], bkspace=100, upper=3, lower=3, /silent)
        speccont = slatec_bvalu(xcont, bkpt, coeff)
    endelse

    if 0 then begin
        plot, xcont, ycont
        oplot, xcont, moogspec_ch*speccont, color=fsc_color('red')
        wait, 1
    endif
    match, xcont, x, wmask, wjunk
    moogspec_ch = moogspec_ch[wmask]*speccont[wmask]
    return, phot eq 1 ? [teff, logg, moogspec_ch] : moogspec_ch
end
