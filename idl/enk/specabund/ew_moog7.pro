function ew_moog7, x, a
    common moog, datalambda, dlam, moogmask, contmask, contdiv, contdivivar, smoothfactor, phot, moogcontrevise, pixmask, hectospec
    common ew, wstart, wend, dyewmp

    teff = a[0]
    logg = a[1]

    moogspec = get_moog7(datalambda, a)

    cgoodwindow = n_elements(wstart)
    ewmoog = dblarr(cgoodwindow)
    ;for i=0,cgoodwindow-1 do ewmoog[i] = total((1.0 - moogspec[wstart[i]:wend[i]])*(dyewmp[wstart[i]:wend[i]])^(-2.0))/total((dyewmp[wstart[i]:wend[i]])^(-2.0))*(datalambda[wend[i]]-datalambda[wstart[i]])
    for i=0,cgoodwindow-1 do ewmoog[i] = int_tabulated(datalambda[wstart[i]:wend[i]], (1.0 - moogspec[wstart[i]:wend[i]])*0.5*(datalambda[wstart[i]+1:wend[i]+1]-datalambda[wstart[i]-1:wend[i]-1])*(dyewmp[wstart[i]:wend[i]])^(-2.0)) / int_tabulated(datalambda[wstart[i]:wend[i]], (dyewmp[wstart[i]:wend[i]])^(-2.0)) * (datalambda[wend[i]]-datalambda[wstart[i]])
    
    return, phot eq 1 ? [teff, logg, ewmoog] : ewmoog
end
