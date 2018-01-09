pro gals_character,rematch=rematch
;note all 'late' objects here are actually early, i got confused in naming the parameters. oops.

;data retreiving
sps = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1,/silent)
cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz',1,/silent)

;match the catalog
nspec = n_elements(sps)
file_matchedcat = 'cl0024_matchedcat.sav'
if file_test(file_matchedcat) eq 0 or keyword_set(rematch) then begin
	matched_cat = lonarr(nspec)
	for i=0,nspec-1 do begin
	    spherematch, cat.ra,cat.dec,sps[i].ra,sps[i].dec,1./3600.,w1,w2
		if n_elements(w1) gt 1 then w1=w1[0]
		matched_cat[i] = w1
	endfor
	save,matched_cat,filename=file_matchedcat
endif else restore,file_matchedcat
cat = cat(matched_cat)

;categorize the fit
wgoodobj = where(sps.goodfit+sps.good eq 2 and sps.zfit ge 0.38 and sps.zfit le 4.0, cgals)
wbadobj  = where(sps.goodfit+sps.good eq 1 and sps.zfit ge 0.38 and sps.zfit le 4.0 and sps.feh ne -999., cbadgals)
wallobj = [wgoodobj,wbadobj]

goodobj = bytarr(nspec)
badobj = bytarr(nspec)
allobj = bytarr(nspec)
goodobj(wgoodobj) = 1
badobj(wbadobj) = 1
allobj(wallobj) = 1

;read in OII EW and find min,max,quartiles
restore,'cl0024_allsci_oii_ew.sav' ;oii_ew_arr,oii_ew_err_arr
;find the type with maximum number of members
cmax=0
cmaxgoodobj=0
for i=0,5 do begin
	inmor = where(cat.morph eq i,cinmor)
	if cinmor gt cmax then cmax = cinmor
	inmorgoodobj = where(cat.morph eq i and allobj eq 1,cinmorgoodobj)
	if cinmorgoodobj gt cmaxgoodobj then cmaxgoodobj=cinmorgoodobj
endfor
data = fltarr(6,cmax)-99
datagoodobj = fltarr(6,cmaxgoodobj)-99
for i=0,5 do begin
	inmor = where(cat.morph eq i,cinmor)
	data[i,0:cinmor-1] = oii_ew_arr(inmor)
	inmorgoodobj = where(cat.morph eq i and allobj eq 1,cinmorgoodobj)
	datagoodobj[i,0:cinmorgoodobj-1] = oii_ew_arr(inmorgoodobj)
endfor 

set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
!p.charsize =1
;oii ew histogram
psname='OII_ew_histogram.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	!p.charsize=1.2
	oii_ew_plot = -1.*oii_ew_arr
	oii_ew_plot(where(oii_ew_plot lt -7.5))=-7
	plothist,oii_ew_plot,bin=2.5,xrange=[-7.5,60],xtitle='OII EW',ytitle='Number of galaxies',/fill,fcolor=fsc_color('gray'),xtickformat='(A1)',yrange=[0,45]
	plothist,oii_ew_plot(where(allobj eq 1 and oii_ew_plot lt 5)),bin=2.5,/overplot,color=fsc_color('org8'),/fill,fcolor=fsc_Color('org4')
	plothist,oii_ew_plot(wallobj),bin=2.5,/overplot,color=fsc_color('org8'),/fline,fcolor=fsc_Color('org8'),forientation=45
	axis,xaxis=0,xrange=[-7.5,60],xstyle=1
	!p.charsize=1
	al_legend,['All spectra','Spectra with SN>8 per Angstrom','Selected sample of quiescent galaxies'],polycolor=['gray','org8','org4'],psym=[10,10,10],box=0,charsize=1.,symsize=2,line_orientation=[-200,45,-200],position=[10,43],polyspace=0.1
device,/close
;morphology
psname='cl0024_morphology.eps'
device, filename = psname,xsize = 15,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	plothist,cat.morph-7,bin=1,xrange=[-1,6]-7,/fill,fcolor=fsc_color('gray'),yrange=[0,60],$
		xstyle=5,xtickformat='(A1)',ytitle='Number of galaxies',charsize=1.2,$
		position=[0.15,0.4,0.95,0.95]
	
	plothist,cat(wallobj).morph-7,bin=1,/overplot,color=fsc_color('org8'),/fill,fcolor=fsc_Color('org4')
	oplot,[2.4,2.4]-7,[49,53.28],thick=7,color=fsc_color('pbg7')
	oplot,[-0.4,-0.4]-7,[49,53.28],thick=7,color=fsc_color('pbg7')
	oplot,[2.13,2.4]-7,[53,53],thick=7,color=fsc_color('pbg7')
	oplot,[-0.4,-0.1]-7,[53,53],thick=7,color=fsc_Color('pbg7')
	xyouts,-0.1-7,51.8,'(early-type galaxies)',charsize=1
	xyouts,0.3-7,55,'Our samples',charsize=1
	axis,xaxis=1,xticks=8,xtickv=indgen(8)-1-7,xtickformat='(A1)'
        cglegend,title=['all spectra','spectra with SN>8','per Angstrom'],color=['gray','org4','org4'],psym=[15,15,0],box=0,charsize=0.8,/data,location=[-3,55],length=0,symsize=2
	cgplot,[0,1],/nodata,yrange=[-70,20],xrange=[-1,6]-7,xtickformat='(A1)',ytitle='OII EW',ystyle=1,$
		xstyle=5,position=[0.15,0.1,0.95,0.4],/noerase,charsize=1.2,yticks=4,ytickv=[-60,-40,-20,0]
	width=((!X.CRange[1] - !X.Crange[0]) / (12)) * 0.75
	cgBoxPlot,data,/overplot,xlocation=(indgen(6))-7,missing_data_value=-99,boxcolor='gray',/fillbox,$
		outliercolor='white',width=width,outlinecolor='black'
	;cgboxplot,datagoodobj,/overplot,xlocation=(indgen(6)+1)*2+1,missing_data_value=-99,boxcolor='org4',/fillbox,$
        ;        outliercolor='white',width=width
	
	galtype=['compact','E','E/S0','S0','Sa+b','S','Sc+d','Irr']
	axis, xaxis=0,xticks=8,xtickv=indgen(8)-8,xtickn=galtype,xtitle='Morphology',charsize=1.2
	axis,xaxis=1,xticks=8,xtickv=indgen(8)-8,xtickformat='(A1)'

device,/close

wgood_late = where(cat(wallobj).morph le 2 and cat(wallobj).morph ge 0,cgood_late)
wgood_late = wallobj(wgood_late)
good_late = bytarr(nspec)
good_late(wgood_late) = 1 

;mass distribution function
set_plot,'x'
;double schecter function from Tomczak2014
qparam = {mstar:10.75,alpha1:-0.47,phi1:10.d^(-2.76),alpha2:-1.97,phi2:10.d^(-5.21)} ;quiescent galaxies
aparam = {mstar:10.78,alpha1:-0.98,phi1:10.d^(-2.54),alpha2:-1.90,phi2:10.d^(-4.29)} ;all galaxies
Marr = findgen(151)/50.+8.5 ;from 8.5 to 11.5 Msun
qMdiff = 10.d^(Marr-qparam.mstar)
qphi = alog(10)*exp(-1.*qMdiff)*qMdiff*(qparam.phi1*(qMdiff^qparam.alpha1)+qparam.phi2*(qmdiff^qparam.alpha2))
aMdiff = 10.d^(Marr-aparam.mstar)
aphi = alog(10)*exp(-1.*aMdiff)*aMdiff*(aparam.phi1*(aMdiff^aparam.alpha1)+aparam.phi2*(amdiff^aparam.alpha2))
;normalize for mass range
qtot = int_tabulated(marr,qphi)
atot = int_tabulated(marr,aphi)
qphi = qphi/qtot
aphi = aphi/atot
;plot histogram
binsize = 0.5

allspec = where(sps.logmstar ge 8.5 and sps.logmstar le 11.5,callspec)
specyhist = histogram(sps(allspec).logmstar,binsize=binsize,nbins=6,min=8.5,locations=specxhist)
allspecphi = specyhist/binsize/total(specyhist)

latespec = where(sps.logmstar ge 8.5 and sps.logmstar le 11.5 and cat.morph le 2 and cat.morph ge 0,clatespec)
specyhist = histogram(sps(latespec).logmstar,binsize=binsize,nbins=6,min=8.5,locations=specxhist)
latespecphi = specyhist/binsize/total(specyhist)

obj = where(sps.logmstar ge 8.5 and sps.logmstar le 11.5 and cat.morph le 2 and allobj eq 1,cobj)
;plothist,sps(obj).logmstar,specxhist,specyhist,bin=binsize,/overplot,color=fsc_color('blue')
specyhist = histogram(sps(obj).logmstar,binsize=binsize,nbins=6,min=8.5,locations=specxhist)
objphi = specyhist/binsize/total(specyhist)

lateobj = where(sps.logmstar ge 8.5 and sps.logmstar le 11.5 and cat.morph le 2 and good_late eq 1,clateobj)
;plothist,sps(lateobj).logmstar,specxhist,specyhist,bin=binsize,/overplot,color=fsc_color('green')
specyhist = histogram(sps(lateobj).logmstar,binsize=binsize,nbins=6,min=8.5,locations=specxhist)
lateobjphi = specyhist/binsize/total(specyhist)

set_plot,'ps'
psname='cl0024_stellarmassfn.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	sunsym = sunsymbol()
	plot,marr,(qphi),xtitle='log[M/M'+sunsym+']',ytitle='normalized fraction',xrange=[8.5,11.5],xstyle=1
;	wzero = where(lateobjphi eq 0,nzero)
;	if nzero gt 0 then lateobjphi(wzero) = 1.e-6
	oplot,[specxhist[0]-binsize/2.,specxhist]+binsize/2.,[lateobjphi[0],lateobjphi]/2.,psym=10
	oplot,marr,(qphi),color=fsc_color('red')
device,/close

set_plot,'ps'
psname='cl0024_magnitudefn.eps'
device, filename = psname,xsize = 12,ysize = 16, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	minmag = 18.
	maxmag = 23.
	binsize=0.5
	numgal_mag_all = histogram(cat.f814w_auto,binsize=binsize,locations=magall,min=minmag,max=maxmag)
	numgal_mag_goodsn = histogram(sps(wallobj).f814w,binsize=binsize,locations=maggoodsn,min=minmag,max=maxmag)
	plothist,cat.f814w_auto,bin=binsize,xrange=[minmag,maxmag],ytitle='Number of galaxies',xstyle=9,xtickformat='(A1)',position=[0.15,0.5,0.95,0.9],/fill,fcolor=fsc_color('gray')
	plothist,sps(wallobj).f814w,bin=binsize,/overplot,color=fsc_color('org8'),/fill,fcolor=fsc_color('org4')
	cglegend,title=['all spectra','spectra with SN>8','per Angstrom'],color=['gray','org4','org4'],psym=[15,15,0],box=0,charsize=0.8,/data,location=[18.2,38],length=0,symsize=2
	axis,xaxis=1,xtickv=magall
	if total(magall-maggoodsn) ne 0 then stop,'Two magnitude arrays are not the same'
	frac_good=float(numgal_mag_goodsn)/float(numgal_mag_all)
	vsym,4,/fill,rot=45
	plot,magall+binsize/2.,frac_good,ytitle='Fraction with SN>8',xrange=[minmag,maxmag],xstyle=1,xtickformat='(A1)',/noerase,position=[0.15,0.3,0.95,0.5],psym=8,yrange=[0,1.1],yticks=4,ytickv=[0.2,0.4,0.6,0.8,1.0],yminor=4
	vsym,24,/fill
	plot,sps(allspec).f814w,sps(allspec).logmstar,psym=8,xrange=[minmag,maxmag],xtitle='F814W',ytitle= 'log[M/M'+sunsym+']',/noerase,position=[0.15,0.10,0.95,0.3],yrange=[8,11.5],ystyle=1,yticks=3,ytickv=[8,9,10,11],yminor=5
device,/close

set_plot,'ps'
psname='cl0024_magnitudefn_passivetype.eps'
device, filename = psname,xsize = 12,ysize = 16, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
        minmag = 18.
        maxmag = 23.
        binsize=0.5
	passivetype = where(oii_ew_arr+oii_ew_err_arr gt -5.)
	goodpassive = where(oii_ew_arr+oii_ew_err_arr gt -5. and allobj eq 1)
        numgal_mag_all = histogram(cat(passivetype).f814w_auto,binsize=binsize,locations=magall,min=minmag,max=maxmag)
        numgal_mag_goodsn = histogram(sps(goodpassive).f814w,binsize=binsize,locations=maggoodsn,min=minmag,max=maxmag)
        plothist,cat.f814w_auto,bin=binsize,xrange=[minmag,maxmag],ytitle='Number of galaxies',xstyle=9,xtickformat='(A1)',position=[0.15,0.5,0.95,0.9],/fill,fcolor=fsc_color('gray'),yrange=[0,48]
	plothist,cat(passivetype).f814w_auto,bin=binsize,/overplot,color=fsc_color('org8'),/fline,fcolor=fsc_Color('org4'),forientation=45
        plothist,sps(goodpassive).f814w,bin=binsize,/overplot,color=fsc_color('org8'),/fill,fcolor=fsc_color('org4')
        al_legend,['All Spectra','Passive galaxies','Passive galaxies with SN>8 per Angstrom'],polycolor=['gray','org4','org4'],psym=10,box=0,charsize=1,position=[18.,48],line_orientation=[-200,45,-200],polyspace=0.1
        axis,xaxis=1,xtickv=magall
        if total(magall-maggoodsn) ne 0 then stop,'Two magnitude arrays are not the same'
        frac_good=float(numgal_mag_goodsn)/float(numgal_mag_all)
        vsym,4,/fill,rot=45
        plot,magall+binsize/2.,frac_good,ytitle='Fraction with SN>8',xrange=[minmag,maxmag],xstyle=1,xtickformat='(A1)',/noerase,position=[0.15,0.3,0.95,0.5],psym=8,yrange=[0,1.1],yticks=4,ytickv=[0.2,0.4,0.6,0.8,1.0],yminor=4
        vsym,24,/fill
	earlytype_goodmass = where(sps.logmstar ge 8.5 and sps.logmstar le 11.5 and cat.morph ge 0 and cat.morph le 2)
        plot,sps(earlytype_goodmass).f814w,sps(earlytype_goodmass).logmstar,psym=8,xrange=[minmag,maxmag],xtitle='F814W',ytitle= 'log[M/M'+sunsym+']',/noerase,position=[0.15,0.10,0.95,0.3],yrange=[8,11.5],ystyle=1,yticks=3,ytickv=[8,9,10,11],yminor=5
device,/close


stop
end
