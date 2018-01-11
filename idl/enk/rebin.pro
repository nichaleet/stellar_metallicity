;+
; NAME: 
;        REBIN
;
; PURPOSE: 
;       - includes ivar weighting
;       - normalizes the spectrum to the median flux between 8000A and
;         9000A
;
; AUTHOR:
;   Evan Kirby
;   UCSC (ekirby@ucolick.org)
;
; CALLING SEQUENCE:
;   rebin, indirectory='/path/to/dat/files', outdirectory='/path/to/rebin/files'
;   
; INPUTS:
;   none
; 
; OUTPUTS:
;   the files specified on line 203
;
; MODIFICATION HISTORY:
;   enk 11nov05
;   enk 22may07    make work with Marla's logarithmically spaced
;                  ASCII files
;-

FUNCTION Standardize_Spectra, $
   lambda, $	; wavelength solution of each spectrum in Angstroms [nspec, *], restframe
   spec, $	; flux, arbitrary units [nspec, *]
   ivar, $ 	; inverse variance, in ivar units of spectra, [nspec, *]
   z

    d_lambdat = 0.33                               ; wavelength spacing in Angstroms of coadded spectrum
    lambdat = d_lambdat * dindgen(10607) + 6000 ; wavelength array of coadded spectrum
    rspec = dblarr(10607) ; individual blueshifted and rebinned spectra [nspec, 10607]
    rivar = dblarr(10607) ; associated inverse variances [nspec, 10607]    
    
    lambda /= 1.+z
    
    w = findex(lambda, lambdat) ; find the index mapping lambda -> lambdat
    ix = floor(w)           ; take the integer part of w
    leftover = w - ix       ; weighting factor from findex
    d_ix = long(-1*ts_diff(ix,1)) ; see where wavelength bin changes
    
    d_ix_0 = where(d_ix eq 0) ; indices where lambdat falls into the same lambda as previous lambdat
    d_ix_1 = where(d_ix eq 1, count1) ; indices where lambdat falls into the next lambda bin
    d_ix_2 = where(d_ix eq 2, count2) ; indices where lambdat falls into two lambda bins
    d_ix_3 = where(d_ix ge 3, count3) ; indices where lambdat falls into two lambda bins
    if count3 ne 0 then begin
        print, 'i = ' + string(i)
        message, 'STANDARIZE_SPECTRA -- Rebinning can not handle input spectrum!!'
    endif
    
    ; indices of the rebinned coadded spectrum for those lambdat which fall uniquely into one lambda
    k = d_ix_0[where(ix[d_ix_0] ge 0 and ix[d_ix_0] le n_elements(spec)-1)]
    j = ix[k]               ; indices of the individual spectra

    rspec[k] = spec[j] ; for lambdat which falls into only one lambda, just copy the counts
    rivar[k] = ivar[j]    ; and the inverse variance

    if count1 gt 0 then begin
        ; indices of the rebinned coadded spectrum for those lambdat which straddle more than one lambda
        k = d_ix_1[where(ix[d_ix_1] ge 0 and ix[d_ix_1] le n_elements(spec)-2)]
        j = ix[k]           ; indices of the individual spectra
        w = where((ivar[j] + ivar[j+1]) ne 0d, countw)
        if countw gt 0 then begin
            k = k[where((ivar[j] + ivar[j+1]) ne 0d)] ; make sure denominators do not equal zero
            j = j[where((ivar[j] + ivar[j+1]) ne 0d)]

            ; weight the counts by ivar and lambda spacing
            rspec[k] = ( spec[j] * ivar[j] * (lambda[j+1] - lambdat[k]) + $
                           spec[j+1] * ivar[j+1] * (lambdat[k+1] - lambda[j+1]) ) $
                         / ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             ivar[j+1] * (lambdat[k+1] - lambda[j+1]) )
            ; do the same for inverse variance (weight ivar by itself)
            rivar[k] = ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                           ivar[j+1] * (lambdat[k+1] - lambda[j+1]))^2. $
                         / (ivar[j] * (lambda[j+1] - lambdat[k])^2. + $
                            ivar[j+1] * (lambdat[k+1] - lambda[j+1])^2.)
        endif
    endif

    if count2 gt 0 then begin
        k = d_ix_2[where(ix[d_ix_2] ge 0 and ix[d_ix_2] le n_elements(spec)-3)]
        j = ix[k]           ; indices of the individual spectra
        w = where((ivar[j] + ivar[j+1] + ivar[j+2]) ne 0d, countw)
        if countw gt 0 then begin
            k = k[w]    ; make sure denominators do not equal zero
            j = j[w]
            
            ; weight the counts by ivar and lambda spacing
            rspec[k] = ( spec[j] * ivar[j] * (lambda[j+1] - lambdat[k]) + $
                           spec[j+1] * ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                           spec[j+2] * ivar[j+2] * (lambdat[k+1] - lambda[j+2]) ) $
                         / ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                             ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                             ivar[j+2] * (lambdat[k+1] - lambda[j+2]))
            ; do the same for inverse variance (weight ivar by itself)
            rivar[k] = ( ivar[j] * (lambda[j+1] - lambdat[k]) + $
                           ivar[j+1] * (lambda[j+2] - lambda[j+1]) + $
                           ivar[j+2] * (lambdat[k+1] - lambda[j+2]))^2. $
                         / (ivar[j] * (lambda[j+1] - lambdat[k])^2. + $
                            ivar[j+1] * (lambda[j+2] - lambda[j+1])^2. + $
                            ivar[j+2] * (lambdat[k+1] - lambda[j+2])^2.)
        endif
    endif
    smooth_ivar, rspec, rivar, 1.55, specfinal, ivarfinal            ; this is where you change the smoothing
    C = median(specfinal[where(lambdat ge 8000 and lambdat le 9000)])
    rspec = specfinal / C
    rivar = ivarfinal * C^2.

    ; return a structure with wavelength, flux, and inverse variance
    rebin = {lambda:lambdat, spec:rspec, ivar:rivar}
    return, rebin
end


pro smooth_ivar, array, ivar, sigma, arrayout, ivarout
     arrayout = array*0.0
     ivarout = ivar*0.0
     if n_elements(array) ne n_elements(ivar) then begin
         print, 'Wrong array dimensions in SMOOTH_IVAR.'
     endif
     asize = n_elements(array)
     n = 31
     half_n = (n-1)/2
     window = (1. / (double(sigma)) * sqrt(2.*!PI)) * exp(-(dindgen(n)-half_n)^2./(2.*double(sigma)^2.))   ;Gaussian
     for i=0,asize-1 do begin
         arrayout[i] = total(array[(i-half_n)>0:(i+half_n)<asize-1]*ivar[(i-half_n)>0:(i+half_n)<asize-1]*window[0>(half_n-i):(n-1)<(asize-i+half_n)]) / total(ivar[(i-half_n)>0:(i+half_n)<asize-1]*window[0>(half_n-i):(n-1)<(asize-i+half_n)])
         ivarout[i] = (total(ivar[(i-half_n)>0:(i+half_n)<asize-1]*window[0>(half_n-i):(n-1)<(asize-i+half_n)]))^2. / total(ivar[(i-half_n)>0:(i+half_n)<asize-1]*(window[0>(half_n-i):(n-1)<(asize-i+half_n)])^2.)
     endfor
end


; This program coadds and plots various subsets of serendip spectra.
pro rebin, indirectory=indirectory, outdirectory=outdirectory, fits=fits
    cl = 2.99792458d5
    searchstring = keyword_set(fits) ? '*spec1d.*.fits' : '*spec1d.*.dat'
    infiles = file_search(indirectory, searchstring, count=c)
    noser = intarr(c)
    objno = lonarr(c)
    for i=0,c-1 do begin
        objnotmp = strsplit(file_basename(infiles[i]), '.', /extract)
        if strmatch(objnotmp[3], '*serendip*') then continue
        objno[i] = objnotmp[3]
        noser[i] = 1
    endfor
    wnoser = where(noser, c)
    infiles = infiles[wnoser]
    objno = objno[wnoser]
    mask = strarr(c)
    for i=0,c-1 do begin
        masktmp = strsplit(file_dirname(infiles[i]), '/', /extract, count=cd)
        if keyword_set(fits) then begin
            mask[i] = strcompress(masktmp[cd-2], /rem)
            mask[i] = strmid(mask[i], 0, strpos(mask[i], 's'))
        endif else mask[i] = strcompress(masktmp[cd-1], /rem)
        mask[i] = strmid(mask[i], 0, strpos(mask[i], '_'))
    endfor
    masku = mask[uniq(mask, sort(mask))]
    nmask = n_elements(masku)
    for i=0,nmask-1 do begin
        w = where(mask eq masku[i], cm)
        outfiles = outdirectory+'/'+masku[i]+'/'+file_basename(infiles[w])+'.rebin'
        masktbl = masku[i]
        ;masktbl3 = strmid(masktbl, 0, 3)        
        ;if masktbl3 eq 'd10' then masktbl = 'd10'
        ;if masktbl3 eq 'd1_' then masktbl = 'd1'
        ;if masktbl3 eq 'd2_' then masktbl = 'd2'
        ;if masktbl3 eq 'd3_' then masktbl = 'd3'
        if masktbl eq 'N147' then masktbl = 'n147'
        if masktbl eq 'N185' then masktbl = 'n185'
        vfile = '/net/pemla/data1/data3/kgilbert/MergeTables/'+masktbl+'.tbl_m'
        if ~file_test(vfile) then continue
        file_mkdir, outdirectory+'/'+masku[i]
        readcol, vfile, id, vobs, format='X,X,L,X,X,X,X,D,X,X,X,X', /silent
        for j=0,cm-1 do begin
            wz = where(id eq objno[w[j]], cwz)
            if cwz ne 1 then begin
                print, 'Skipping '+strcompress(j+1,/rem)+' of '+strcompress(cm, /rem)+': '+infiles[w[j]]+'; no vobs found'
                continue
            endif
            if keyword_set(fits) then begin
                spec1db = mrdfits(infiles[w[j]], 1, /silent)
                spec1dr = mrdfits(infiles[w[j]], 2, /silent)
                lambda = [spec1db.lambda, spec1dr.lambda]
                spec = [spec1db.spec, spec1dr.spec]
                ivar = [spec1db.ivar, spec1dr.ivar]
            endif else begin
                readcol, infiles[w[j]], lambda, spec, ivar, format='D,D,D', /silent
            endelse
            wz = wz[0]
            print, strcompress(j+1,/rem)+' of '+strcompress(cm, /rem)+': '+infiles[w[j]]
            rebin = Standardize_Spectra(lambda, spec, ivar, vobs[wz]/cl)

            openw, 1, outfiles[j], width=50000
            for k=0,n_elements(rebin.lambda)-1 do begin
                if rebin.spec[k] eq 0 then continue
                printf, 1, rebin.lambda[k], rebin.spec[k], rebin.ivar[k]
            endfor
            close, 1
        endfor
    endfor
end
