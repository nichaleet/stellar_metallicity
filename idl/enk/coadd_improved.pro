
; NAME:
;  coadd_spectra
;
; PURPOSE:
;  sums 1d spectra with ivar weighting to get composite having
;  improved SNR
;
; CALLING SEQUENCE:
;  result = coadd_spectra(list,normalize=normalize,zero=zero,inorm=inorm,ivartrim=ivartrim,weight=weight)
; 
; INPUTS:
;  list -- list of zresult structures tabulating files to coadd
;  
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;  normalize -- if set, spectrum is normalized to unity prior to
;               coaddition
;  zeroweight -- if set, Aband and Bband regions are given 0 ivar
;                (makes little difference)
;  ivartrim   -- if set, throw out points with worst ivar, as used to
;                do (can be problematic - DEPRECATED)
;  inorm   -- normalize all ivar to have mean 1, and include only
;             variations in ivar due to sky intensity, not
;             obj. intensity (~equally weight
;                each spectrum)
;  dloglam -- change sample spacing (default: 2E-5)
;  weight  -- a weight to apply to each element of list in building
;             the coadded spectra.  If weights are used, NORMALIZE and
;             INORM are also set.
;
; OUTPUTS:
;   result -- a structure of the summed spectra with tags for spectra,
;             lambda, ivar
; 
; COMMENTS:  
;   this routine is useful for improving SNR by coadding many files.
;
; REVISION HISTORY:
;  md, rs 12dec02
;  md  29jan03 ivar=0. in Aband, Bband regions
; jan, ry 2004, 2005 major changes
; enk 2009may28 change skyspec->ivar mapping into
; skyspec,coaddedspec->ivar mapping
;----------------------------------------------------------------------

function coadd_improved, moogify, z, normalize=normalize, zeroweight=zeroweight, $
                        ivartrim=ivartrim,inorm=inorm,dloglam=dloglam,weight=weight

    nspec = n_elements(moogify)
    rlamlimit = [6300., 9100.]
    llimit = alog10(rlamlimit)

    if n_elements(dloglam) eq 0 then dloglam = 2.e-5
    npoints = long((llimit[1]-llimit[0])/dloglam)
    specsum = dblarr(npoints)
    norm=fltarr(npoints)
    llambda = findgen(npoints)*dloglam +llimit[0]
    lambda = 10.^llambda
    ivar = lambda*0.d0
    if n_elements(ivartrim) eq 0 then ivartrim = 0
    if n_elements(inorm) eq 0 then inorm = 0
    if n_elements(normalize) eq 0 then normalize=0
    useweight=n_elements(weight) gt 0
    if n_elements(weight) lt n_elements(moogify) then weight=fltarr(n_elements(moogify))+1.
    weight=weight/total(weight)
    totweight=specsum*1.d0
    totn=totweight
    inorm=inorm OR useweight
    normalize = normalize OR useweight

    rescale=fltarr(nspec)

    iter = 0
    while iter lt 10 do begin
        for i=0, nspec-1 do begin ;begin loop over input list

            if n_params() eq 2 then begin
                z_in = z[i]
                if file_test(moogify[i]) then begin
                    if inorm eq 0 then begin 
                        ss1d=fill_gap(moogify[i],/tweak,/telluric,/silent,header=header) 
                        wave=ss1d.lambda
                        airtovac,wave
                        ss1d.lambda=wave
                    endif else ss1d=fill_gap(moogify[i],/tweak,/silent,header=header)
                endif else begin
                    print, strtrim(moogify[i], 2)+' not found'
                    continue
                endelse
                spec_in=ss1d.spec
                skyspec_in=ss1d.skyspec
                lambda_in=ss1d.lambda
                ivar_in=ss1d.ivar

                badpix = ss1d.nbadpix gt 2 OR (ss1d.infomask AND 21b) gt 0 $
                         OR (ss1d.ormask AND 30b) gt 0 

                badpix=dilate(badpix,intarr(5)+1)

                ivar_in=ivar_in*(badpix eq 0)     
                ss1d.ivar=ivar_in
            endif else begin
                spec_in = moogify[i].contdiv
                lambda_in = moogify[i].lambda
                ivar_in = moogify[i].contdivivar
                skyspec_in = moogify[i].skyspec
                z_in = moogify[i].zrest

                badpix = ~finite(ivar_in) or ~finite(spec_in) or ~finite(skyspec_in) or ivar_in le 0.0

                winf = where(~finite(ivar_in) or ~finite(spec_in) or ~finite(skyspec_in), cinf)
                if cinf gt 0 then begin
                    ivar_in[winf] = 0.0
                    spec_in[winf] = 0.0
                    skyspec_in[winf] = 0.0
                endif
            endelse



; corrected for DEIMOS response in fill_gap or inorm section

; sometimes want to give spectra ~equal weight, but deweight bad
; pixels - set IVAR to have median 1 in that case, and have it scale
;          ONLY with sky spectrum
                if inorm OR useweight then begin
; do two sides separately



; fit for IVAR vs. sky counts
                    if iter eq 0 then begin
                        wh=where(ivar_in NE 0 AND skyspec_in gt 2*median(skyspec_in) $
                                 AND badpix eq 0 AND $
                                 lindgen(n_elements(skyspec_in)) $
                                 lt n_elements(skyspec_in)/2 ,ct)
                        if ct gt 10 then begin
                            fit=svdfit(skyspec_in[wh],1./ivar_in[wh],2)
                            wh1=where(lindgen(n_elements(skyspec_in)) $
                                      lt n_elements(skyspec_in)/2 ,ct)
; adjust ivar
                            if ct gt 0 then ivar_in[wh1]= $
                                1./(fit[0]+fit[1]*(skyspec_in[wh1]>0.))*(ivar_in[wh1] ne 0)
                        endif    



; now do it on the other side
                        wh=where(ivar_in NE 0 AND skyspec_in gt 2*median(skyspec_in) $
                                 AND badpix eq 0 AND $
                                 lindgen(n_elements(skyspec_in)) $
                                 ge n_elements(skyspec_in)/2 ,ct)
                        if ct gt 10 then begin
                            fit=svdfit(skyspec_in[wh],1./ivar_in[wh],2)
                            wh2=where(lindgen(n_elements(skyspec_in)) $
                                      ge n_elements(skyspec_in)/2 ,ct)
                            if ct gt 0 then ivar_in[wh2]= $
                                1./(fit[0]+fit[1]*(skyspec_in[wh2]>0.))*(ivar_in[wh2] ne 0)
                        endif    
                    endif else begin
                        specsuminterpol = interpol(specsum,lambda*(1.+z_in),lambda_in)
                        wh=where(ivar_in NE 0 $
                                 AND badpix eq 0 AND $
                                 lindgen(n_elements(skyspec_in)) $
                                 lt n_elements(skyspec_in)/2 and specsuminterpol gt 0,ct)
                        if ct gt 10 then begin
                            fit=regress(transpose([[skyspec_in[wh]], [specsuminterpol[wh]]]),1./ivar_in[wh],/double, const=const)
                            wh1=where(lindgen(n_elements(skyspec_in)) $
                                      lt n_elements(skyspec_in)/2 ,ct)
; adjust ivar
                            if ct gt 0 then begin
                                ivar_in[wh1] = 1./(const+fit[0]*(skyspec_in[wh1]>0.)+fit[1]*(specsuminterpol[wh1]))*(ivar_in[wh1] ne 0)
                                ;print, 'iteration ', iter, ', blue side fit: ', const, fit, format='(A10,I2,A17,3(G10.3))'
                            endif
                        endif    



; now do it on the other side
                        wh=where(ivar_in NE 0 $
                                 AND badpix eq 0 AND $
                                 lindgen(n_elements(skyspec_in)) $
                                 ge n_elements(skyspec_in)/2 and specsuminterpol gt 0,ct)
                        if ct gt 10 then begin
                            fit=regress(transpose([[skyspec_in[wh]], [specsuminterpol[wh]]]),1./ivar_in[wh],/double, const=const)
                            wh1=where(lindgen(n_elements(skyspec_in)) $
                                      lt n_elements(skyspec_in)/2 ,ct)
; adjust ivar
                            if ct gt 0 then begin
                                ivar_in[wh1] = 1./(const+fit[0]*(skyspec_in[wh1]>0.)+fit[1]*(specsuminterpol[wh1]))*(ivar_in[wh1] ne 0)
                                if i eq 0 then print, 'iteration ', iter, ', red side fit:  ', const, fit, format='(A10,I2,A17,3(G10.3))'
                            endif
                        endif    
                    endelse

                    ivar_in=ivar_in*(skyspec_in ne 0) 

                    if total(ivar_in lt 0) gt 0 then ivar_in=ivar_in*0.
                    whok=where(ivar_in ne 0,ct)
                    if ct gt 0 then rescale[i]=1./median(ivar_in[whok]) else rescale[i]=1.
                    ivar_in=ivar_in*rescale[i]

; need to do telluric corrections after adjusting ivar for best results
                    if n_params() eq 2 then begin
                        airmass = sxpar(header, 'AIRMASS')
                        ss1d.ivar=ivar_in
                        ss1dtmp=ss1d


                        remove_telluric, ss1d, airmass,silent=silent
                        fix_response,ss1d
                        
                        ivar_in=ss1d.ivar
                        wave=ss1d.lambda

                        airtovac,wave
                        ss1d.lambda=wave
                        lambda_in=ss1d.lambda
                        spec_in=ss1d.spec
                    endif
                endif



; set to unit average flux if desired
                if normalize OR useweight then begin
                    whgood=where(ivar_in gt 0,ct)
                    if ct gt 1 then norm=median(spec_in[whgood]) else norm=-1.E10
                    if norm lt 0 then ivar_in=ivar_in*0.
                    spec_in = spec_in/norm   ;normalize
                    ivar_in = ivar_in*norm^2 ;keep snr constant

                    if n_params() eq 2 then begin
                        ss1d.spec = spec_in
                        ss1d.ivar=ivar_in
                    endif
                endif





; the data has noise in excess of Poisson due to skysub problems. 
; Get rid of the worst 5% of the pixels, plus their immediate
; neighbors
; THIS IS OLD CODE AND SHOULD BE A LAST RESORT!

                if ivartrim then begin

                    awful = where(ivar_in lt 1.e-8)
                    ivar_in[awful] = 0.
                    bad = (ivar_in*0.) 
                    good = bad             ; how is that?
                    nlam = n_elements(lambda_in)
                    cut = sort(ivar_in)
                    junk = min(ivar_in[cut] gt 0, cs) ;find location of worst non-zero point
                    cut = cut[cs:cs+nlam/20]
                    bad[cut] = 1.
                    bad = smooth(bad, 3)
                    good = bad eq 0.
                    ivar_in = ivar_in*good
                endif



; account for interpolation issues
                ivar_in = ivar_in < smooth(double(ivar_in),3)



;before coadding, set regions of telluric absorption to have ivar=0
                if keyword_set(zeroweight) then begin
                    aband = where(lambda_in gt 7592. and lambda_in lt 7675., naband)
                    bband = where(lambda_in gt 6866. and lambda_in lt 6881., nbband)
                    if (naband gt 0) then ivar_in[aband] = 0.
                    if (nbband gt 0) then ivar_in[bband] = 0.
                endif





                lambdarange=minmax(lambda_in/(1.+z_in))
                defined=lambda gt (lambdarange[0]+1) AND lambda lt (lambdarange[1]-1)
                objflux=defined*interpol(spec_in,lambda_in/(1.+z_in),lambda) ;,/spline)
                objivar=defined*interpol(ivar_in,lambda_in/(1.+z_in),lambda) ;,/spline)

                if inorm then objnorm = defined*rescale[i]
                totn=totn+(objivar gt 0)

                whdefined=where(defined gt 0)

                if useweight then ivar =  ivar + $
                                          objivar*weight[i]  $
                else ivar =  ivar + $
                             objivar  ;accumulate ivar

                if useweight then specsum =  specsum + $
                                             objflux*objivar*weight[i] $
                else specsum =  specsum + $
                                objflux*objivar

                if useweight then totweight=totweight+weight[i]*(objivar ne 0)
                
                if inorm then norm=norm + objnorm
                
            endfor
;
; now normalize 
        good = where(ivar gt 0)
        
        specsum[good] = specsum[good]/ivar[good]

        if inorm then ivar[good]=ivar[good]/norm[good]*totn[good]

        if useweight then ivar[good]=ivar[good]/totweight[good]
        iter++
    endwhile

    w = where(ivar gt 0)
    specsum = specsum[w]
    lambda = lambda[w]
    ivar = ivar[w]
    totn = totn[w]
    result = { spec: specsum, lambda: lambda, ivar: ivar, totn:totn}

    return, result
end
