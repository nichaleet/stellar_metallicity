
; analogous to IDL interpol, but for coeffs of a tset
; this actually does a linear fit, but COULD call interpol?
function tset_coeff_interpol, tset, X, U

  cf   = tset.coeff
  cfdim = (size(cf, /dimens))                                                  ;ENK
  ncf = cfdim[0]                                                               ;ENK
  if n_elements(cfdim) eq 1 then nset = 1 else nset = cfdim[1]                 ;ENK
  ;ncf  = (size(cf, /dimens))[0]
  ;nset = (size(cf, /dimens))[1]

  if nset NE n_elements(X) then message, 'dimension mismatch!'
  
  newcoeff = dblarr(ncf, n_elements(U))
  for i=0, ncf-1 do begin 
     lf = linfit(X, cf[i, *])
     newcoeff[i, *] = lf[0] + lf[1]*U
  endfor 

  return, newcoeff
end
