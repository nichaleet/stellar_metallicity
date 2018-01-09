pro continuumnormalize_deimos,str,z,hdr
    lambda = str.lambda
    skyspec = str.skyspec
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit skyline
    fit = fitskyline(lambda,skyspec)
    ;spec resolution
    specres,fit,str
    ;telluric
    airmass = sxpar(hdr,'airmass')
    str=create_struct(str,'airmass',airmass)
    dlam = str.dlam 
    telluric,str

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   lambdarest = str.lambda/(1.+z)
   spec = str.telldiv
   ivar = str.ivar
   telldivivar = str.telldivivar
 
   ;create continuum mask
   npix = n_elements(lambda)
   contmask = bytarr(npix)+1
   readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
   for i=0,n_elements(linestart)-1 do begin
      w = where(lambdarest ge linestart[i] and lambdarest le lineend[i], c)
      if c gt 0 then contmask[w] = 0
   endfor

   w = where(~finite(spec) or ~finite(lambda) or ivar eq 0, c)
   if c gt 0 then contmask[w]=0

   tellmask = bytarr(n_elements(lambda))
   tellstart = [6864., 7591., 8938.]
   tellend = [6935., 7694., 9500.]
   for i=0,n_elements(tellstart)-1 do begin
       w = where(lambda ge tellstart[i] and lambda le tellend[i], c)
       if c gt 0 then begin
           ;contmask[w] = 0
           tellmask[w] = 1
       endif
   endfor
   nhalf = fix(npix/2)
   contmask[nhalf-3:nhalf+3]=0
   contmask[0:2] = 0
   contmask[npix-3:npix-1] = 0
   str=create_struct(str,'contmask',contmask)
   str=create_struct(str,'tellmask',tellmask)
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;continuum normalize
   won = where(contmask eq 1, complement=woff, con)
   bkspace = 150.
   bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=telldivivar[won], bkspace=bkspace, upper=3, lower=3, /silent)
   if bkpt[0] eq -1 then stop,'cannot do continnum fit'
   cont = slatec_bvalu(lambda, bkpt, coeff)
   contdiv = str.telldiv/cont
   contdivivar = str.telldivivar*cont^2
      
   str=create_struct(str,'continuum',cont)
   str=create_struct(str,'contdivivar',contdivivar)
   str=create_struct(str,'contdiv',contdiv)
   str=create_struct(str,'sn',0.)
   
   ;plot,str.lambda,str.spec,yrange=cgpercentiles(str.spec,percentiles=[0.02,0.98])
   ;oplot,str.lambda,str.telldiv,color=fsc_color('pink')
   ;oplot,str.lambda,str.continuum,color=fsc_color('yellow')
   ;stop
   wcont = where(str.contmask[3:npix-4] eq 1)+3
   wcont = wcont[where(finite(str.telldiv[wcont]) and finite(str.continuum[wcont]) and str.continuum[wcont] ne 0)]
   dev = abs((str.telldiv[wcont] - str.continuum[wcont]) / str.continuum[wcont])
   avgdev = mean(dev)
   w = where(dev lt 3.0*avgdev, c)
   if c gt 0 then str.sn = 1.0/mean(dev[w])
end
