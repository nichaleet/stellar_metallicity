function join_deimos_blu_red,strblu,strred
  npixblu = n_elements(strblu.lambda)
  npixred = n_elements(strred.lambda)
  npix = npixblu+npixred
  outstr = {spec:fltarr(npix),lambda:fltarr(npix),ivar:fltarr(npix),crmask:fltarr(npix),$
            bitmask:fltarr(npix),ormask:fltarr(npix),nbadpix:fltarr(npix),infomask:intarr(npix),$
            objpos:strblu.objpos,fwhm:strblu.fwhm,nsigma:strblu.nsigma,R1:strblu.r1,$
            R2:strblu.r2,skyspec:fltarr(npix),ivarfudge:strblu.ivarfudge}
  outstr.spec    = [strblu.spec,strred.spec]
  outstr.lambda  = [strblu.lambda,strred.lambda]
  outstr.ivar    = [strblu.ivar,strred.ivar]
  outstr.crmask  = [strblu.crmask,strred.crmask]
  outstr.bitmask = [strblu.bitmask,strred.bitmask]
  outstr.ormask  = [strblu.ormask,strred.ormask]
  outstr.nbadpix = [strblu.nbadpix,strred.nbadpix]
  outstr.infomask= [strblu.infomask,strred.infomask]
  outstr.skyspec = [strblu.skyspec,strred.skyspec]

  return,outstr

end
