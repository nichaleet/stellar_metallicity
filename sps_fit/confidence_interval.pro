function confidence_interval,x,prob_x,plot=plot,xtitle=xtitle
;return peak,first q,median, third q
valx = dblarr(6)
numx = n_elements(x)

;peak
peakval = max(prob_x,pkloc)
valx[0] = x[pkloc]

;make CDF (cumulative distribution function)
CDF = dblarr(numx)
for i=1,numx-1 do CDF[i]=int_tabulated(x[0:i],prob_x[0:i]) 
if ~keyword_set(xtitle) then xtitle=''
if keyword_set(plot) then plot,x,cdf,ytitle='CDF',xtitle=xtitle
valx[1:3] = interpol(x,cdf,[0.16,0.5,0.84])

;Check CDF of the peak value
best_cdf = cdf(pkloc)
if best_cdf ge 0.34 and best_cdf le 0.66 then valx[4:5] = interpol(x,cdf,[(best_cdf-0.34)>0,(best_cdf+0.34)<1.]) else $
if best_cdf lt 0.34 then valx[4:5] = [0.,interpol(x,cdf,best_cdf+0.34)] else $
if best_cdf gt 0.66 then valx[4:5] = [interpol(x,cdf,best_cdf-0.34),0.4] 

return,valx
end
