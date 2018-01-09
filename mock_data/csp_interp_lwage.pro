function csp_interp_lwage, age,lambdain
;return light weighted age at lambdain of CSP with elasped time age (in Gyr)
  common csp_spec_lwage, csp_lw, cspage_lw
      if (size(csp_lw))[1] eq 0 then begin
      ;;tau model with no dust. tau = 1 Gyr
      csp_lw = sps_read_spec('/scr2/nichal/workspace2/fsps-master/OUTPUTS/CSP_try4a.out.spec')
  endif

  cspage_lw = csp_lw[uniq(csp_lw.agegyr, sort(csp_lw.agegyr))].agegyr
  if age lt min(cspage_lw) or age gt max(cspage_lw) then message, 'Age is off grid'
  ncspage_lw = n_elements(cspage_lw)

  wa = value_locate(cspage_lw, age)
  if cspage_lw[wa] eq age then begin
    ia = wa
    da = 1d
    na = 1
  endif else begin
    if wa eq ncspage_lw-1 then wa -= 1
    ia = [wa, wa+1]
    da = abs(reverse(cspage_lw[ia])-age)
    na = 2
  endelse
  cspspec = dblarr(5994)
  for i=0, na-1 do begin
     w = where(csp_lw.agegyr eq cspage_lw[ia[i]],c)
     if i eq 0 then begin
       lambda = csp_lw[w].lambda
       cspspeci = csp_lw[w].spec
     endif else cspspeci = interpol(csp_lw[w].spec,csp_lw[w].lambda,lambda)
     cspspec += da[i]*cspspeci
  endfor
  cspspec /= na eq 2 ? (cspage_lw[ia[1]]-cspage_lw[ia[0]]) : 1d
  ageout = interpol(cspspec,lambda,lambdain)
  return,ageout
end

