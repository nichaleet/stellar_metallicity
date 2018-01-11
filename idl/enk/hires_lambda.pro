function hires_lambda, hdr
    npix = sxpar(hdr, 'NAXIS1')
    norders = sxpar(hdr, 'NAXIS2')
    lambda = dblarr(npix, norders)
    for i=0,norders-1 do begin
        lambdapoly = dblarr(7)
        poly_string1 = sxpar(hdr, 'WV_0_'+string(i+1, format='(I02)'))
        poly_string2 = sxpar(hdr, 'WV_4_'+string(i+1, format='(I02)'))
        lambdapoly[0] = double(strmid(poly_string1, 0, 17)) 
        lambdapoly[1] = double(strmid(poly_string1, 17, 17))
        lambdapoly[2] = double(strmid(poly_string1, 34, 17))
        lambdapoly[3] = double(strmid(poly_string1, 51, 17))
        lambdapoly[4] = double(strmid(poly_string2, 0, 17))
        lambdapoly[5] = double(strmid(poly_string2, 17, 17))
        lambdapoly[6] = double(strmid(poly_string2, 34, 17))
        lambda[*,i] = poly(dindgen(npix)+1d, lambdapoly)
    endfor
    return, lambda
end
