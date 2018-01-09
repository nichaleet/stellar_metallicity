pro cannon_makesps_struct,wlmin=wlmin,wlmax=wlmax,fileout=fileout
    if n_elements(wlmin) eq 0 then wlmin=3500.
    if n_elements(wlmax) eq 0 then wlmax=7500.
    if n_Elements(fileout) eq 0 then fileout='cannon_sps_str.fits'

    f = file_search('/scr2/nichal/workspace2/fsps-3.0/SSP/SSP_Padova_MILES_Kroupa_*.spec')
    spso = sps_read_spec(f) ;[94,22] array in age and metallicity
    dimen = size(spso,/dimension)
    ;make new structure with interested wl and with continuum normalized spectra
    lambda = spso[0].lambda
    goodwl = where(lambda gt 3500. and lambda lt 7500.,cwl)
    lambda = lambda(goodwl)
    sps = {zmet:0.,agegyr:0.,logmass:0.,loglbol:0.,logsfr:0.,lambda:fltarr(cwl),$
           spec:fltarr(cwl),cont:fltarr(cwl),contdiv:fltarr(cwl)}
    sps = replicate(sps,dimen(0),dimen(1))
    
    sps.zmet = spso.zmet
    sps.agegyr = spso.agegyr
    sps.logmass = spso.logmass
    sps.loglbol = spso.loglbol
    sps.logsfr = spso.logsfr

    ;continuum mask
    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    contmask = bytarr(cwl)+1
    for i=0,n_elements(linestart)-1 do begin
       w = where(lambda ge linestart[i] and lambda le lineend[i], c)
       if c gt 0 then contmask[w] = 0
    endfor
    
    for i=0,dimen(0)-1 do for j=0,dimen(1)-1 do begin
        spsnow = sps[i,j]
        spsonow = spso[i,j]
        spsnow.lambda = spsonow.lambda(goodwl)
        if total(spsnow.lambda-lambda) ne 0 then stop,'wavelength is wrong'
        spsnow.spec = spsonow.spec(goodwl)
        spsspec = spsnow.spec
        n = cwl
        won = where(contmask eq 1, complement=woff, con) 
        ;;fit continuum to rest wavelength
        bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
        cont = slatec_bvalu(lambda, bkpt, coeff)
        spsspec /= cont
        spsnow.cont = cont
        spsnow.contdiv = spsspec
        sps[i,j] = spsnow
    endfor
    mwrfits,sps,fileout,/create,/silent
end

