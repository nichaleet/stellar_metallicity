pro ks_test,x1,y1,x2,y2,D_ks,Neff_ks,pvalue,silent=silent,noplot=noplot
;e.g ks_test,[[sdss.logmstar],[sdss.feh]],[[sci(wpassivemem).logmstar],[sci(wpassivemem).feh]],D_ks,Neff_ks

  ;make sure that data points (y) are positive
  miny = min([y1,y2])
  if miny lt 0. then begin
     y1 += miny
     y2 += miny
  endif

  ;sort data
  sort1 = sort(x1)
  x1 = x1(sort1)
  y1 = y1(sort1)
  sort2 = sort(x2)
  x2 = x2(sort2)
  y2 = y2(sort2)
	
  ;cumulative sum
  y1cum = total(y1,/cumulative)/total(y1)
  y2cum = total(y2,/cumulative)/total(y2)

  ;find effective N
  n1 = float(n_Elements(x1))
  n2 = float(n_elements(x2))
  neff = (n1*n2)/(n1+n2)

  ;find D
  nsam = 2.*neff
  minx = min([x1,x2])
  maxx = max([x1,x2])
  xsam = findgen(fix(nsam))/fix(nsam)*(maxx-minx)+minx
  model1 = interpol(y1cum,x1,xsam)
  model2 = interpol(y2cum,x2,xsam)
  D_ks = max(abs(model1-model2))
  ;find p value
  prob_ks,D_ks,Neff,pvalue
  if ~keyword_set(noplot) then begin
  set_plot,'x'
  plot,x1,y1cum,xtitle='x',ytitle='y'
  oplot,x1,y1cum,color=fsc_color('blue')
  oplot,x2,y2cum,color=fsc_color('red')
  oplot,xsam,model1,color=fsc_color('cyan'),linestyle=3
  oplot,xsam,model2,color=fsc_color('salmon'),linestyle=3
  endif
  if ~keyword_set(silent) then begin
     print, 'Dvalue =', D_ks
     print, 'P Value=', pvalue
  endif
;stop
  set_plot,'ps'
  
end

end
