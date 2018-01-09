pro stack_spec

  dir = '/scr2/nichal.re/workspace2/sps_fit/data/'
  masks = 'Cl1604_'+['B','C','D1','D2']
  mstr=[]
  for i=0,n_elements(masks)-1 do begin 
     strnow  = mrdfits(dir+masks[i]+'/sps_fit.fits.gz',1)
     mstr = [mstr,strnow]
  endfor
  
;sort by mass from low to high
  order = bsort(mstr.logmstar)
  mstr = mstr(order)
  ngals = n_elements(mstr)
;stack of 10 each
  nbins = ceil(ngals/10.)
  for i=0,nbins-1 do begin
     strnow = mstr[i*10:i*10+9<ngals-1]
     spec = strnow.contdiv
     ivar = strnow.contdivivar
     lambda = strnow.lambda
     wlrange = minmax(lambda)
     
  endfor
end
