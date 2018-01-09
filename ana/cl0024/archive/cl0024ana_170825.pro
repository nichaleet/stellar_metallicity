pro mzrcurve, x, a, f, pder
	f=a[0]-alog10(1.+10.^((x-a[1])*(-1.*a[2])))-8.9
	; FEH = 12+log(O/H)-8.9 = z0-log[1+(M*/M0)^-gmma]-8.9
	;pder=[[replicate(1.0, N_ELEMENTS(X))],[a[2]/(alog(10)*(a[2]*(a[1]-x)))],[(x-a[1])/(alog(10)*(a[2]*(a[1]-x)))]]
	pder=[[replicate(1.0, N_ELEMENTS(X))],[replicate(1.0, N_ELEMENTS(X))],[replicate(1.0, N_ELEMENTS(X))]]
end

pro cl0024ana
	;main program for analysis section of the paper. combine previous plot_mass*.pros
	;need the matched catalog from gals_character.pro
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;data retreiving
	sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1,/silent)
	cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz',1,/silent)
	file_matchedcat = 'cl0024_matchedcat.sav'

	restore,file_matchedcat
	cat = cat(matched_cat)
	nsci = n_elements(sci)
 	;check the match
	for i=0,nsci-1 do begin 
		if randomu(seed) lt 0.2 then begin
			gcirc,2,cat(i).ra,cat(i).dec,sci(i).ra,sci(i).dec,dis
			if dis gt 1 then stop,'check the match of ra/dec'  ;if the distance is greater than 1 arcsecond
		endif
	endfor 
	print,'Roughly checked the match between the catalog and science catalog.'
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;categorize the fit
	wgoodmem = where(sci.goodfit+sci.good eq 2 and sci.zfit ge 0.38 and sci.zfit le 4.0 and sci.logmstar gt 7., cgals)
	wbadmem  = where(sci.goodfit+sci.good eq 1 and sci.zfit ge 0.38 and sci.zfit le 4.0 and sci.feh ne -999. and sci.logmstar gt 7., cbadgals)
	wallmem = [wgoodmem,wbadmem]

	goodmem = bytarr(nsci)
	goodmem(wgoodmem) = 1
	badmem = bytarr(nsci)
	badmem(wbadmem) = 1
	allmem = bytarr(nsci)
	allmem(wallmem) = 1
	
	wearlymem = where(cat(wallmem).morph le 2 and cat(wallmem).morph ge 0,cearlymem)
	wearlymem = wallmem(wearlymem)
	earlymem = bytarr(nsci)
	earlymem(wearlymem) = 1

	wlatemem = where(cat(wallmem).morph le 5 and cat(wallmem).morph ge 3,clatemem)
	wlatemem = wallmem(wlatemem)
	latemem = bytarr(nsci)
	latemem(wlatemem) = 1

	wgoodearlymem = where(earlymem eq 1 and goodmem eq 1,cgoodearlymem)
	wbadearlymem = where(earlymem eq 1 and badmem eq 1,cbadearlymem)

	ageform = fltarr(nsci)+999.
	ageform(wearlymem) = (galage(sci(wearlymem).zfit,1000.)/1.e9-sci(wearlymem).age)>0. ;age of universe when it was formed
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	;misc
	set_plot,'ps'
	!p.multi = [0,1,1]
	!p.font = 0
	sunsym = sunsymbol()
	Delta = '!9'+string("104B)+'!x'
	alpha = '!9'+string("141B)+'!x'

	;fix the uncertainties
	get_additional_uncertainties,sci(wearlymem).feh,sci(wearlymem).sn,deltaage,deltafeh
	sci(wearlymem).feherr = deltafeh; sqrt(sci(wearlymem).feherr^2+deltafeh^2)
	ageupper = 10^deltaage*sci(wearlymem).age
	agelower = sci(wearlymem).age/10^deltaage
	sci(wearlymem).ageerr = 0.5*(ageupper-agelower)
	;hmm
	mass = sci(wearlymem).logmstar	
	feh = sci(wearlymem).feh
	feherr = sci(wearlymem).feherr
	feherr(where(feherr eq 0.)) = median(feherr)
	sci(wearlymem).feherr = feherr
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Averaging
	massbins = [8.9,9.3,9.6,10.,10.25,10.55,11.]
	nbins = n_elements(massbins)-1

	ave_mass    = fltarr(nbins)
	ave_mass_dev= fltarr(nbins)
	ave_feh     = fltarr(nbins) 
	ave_feh_dev = fltarr(nbins)
	bndry_mass = fltarr(nbins*2)
	for i=0,nbins-1 do begin
		if i eq nbins-1 then msel = where(mass gt massbins[i],cmsel) else $
			msel = where(mass gt massbins[i] and mass lt massbins[i+1],cmsel)
		meanerr,feh(msel),feherr(msel),fehmean,sigmam,sigmad,sigmas		
		ave_feh(i) = fehmean
		ave_feh_dev(i) = sigmas
		ave_mass(i) = mean(mass(msel))
		ave_mass_dev(i) = stdev(mass(msel))
		bndry_mass[i*2:i*2+1] = [massbins[i],massbins[i+1]]
	endfor
	hifeh = interpol(ave_feh+ave_feh_dev,ave_mass,bndry_mass)
	lofeh = interpol(ave_feh-ave_feh_dev,ave_mass,bndry_mass)
	wtoolow = where(lofeh lt -1.,cwtoolow)
	if cwtoolow ge 1 then lofeh(wtoolow) = -1.
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Get OII emission line equivalent width
	;oii_ew_arr = fltarr(cearlymem)
	;oii_ew_err_arr = fltarr(cearlymem)
	set_plot,'x'
	;for i=0,cearlymem-1 do begin
	;	oii_ew_arr(i)=oii_ew(sci(wearlymem[i]).contdiv,sci(wearlymem[i]).lambda,sci(wearlymem[i]).contdivivar,sci(wearlymem[i]).zspec,widtherr=widtherr)
	;	oii_ew_err_arr(i) = widtherr
	;endfor
	;save,oii_ew_arr,oii_ew_err_arr,filename='cl0024_earlymem_oii_ew.sav'
	oii_ew_arr = fltarr(nsci)
	oii_ew_err_arr = fltarr(nsci)
	for i=0,nsci-1 do begin
		if earlymem(i) eq 1 then zin = sci(i).zspec else zin = sci(i).z
		oii_ew_arr(i)=oii_ew(sci(i).contdiv,sci(i).lambda,sci(i).contdivivar,zin,widtherr=widtherr)
		oii_ew_err_arr(i) = widtherr
		;stop
	endfor
	save,oii_ew_arr,oii_ew_err_arr,filename='cl0024_allsci_oii_ew.sav'
	set_plot,'ps'
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Get SDSS data
	sdss = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass2/sps_fit.fits.gz',1,/silent)
	nofit = where(sdss.feh eq -999 ,cnofit)
	if cnofit gt 0 then remove,nofit,sdss
	sdss(where(sdss.feherr eq 0.)).feherr = median(sdss.feherr)
	sdss.feherr = sqrt(sdss.feherr^2+0.1^2)
	restore,'/scr2/nichal/workspace2/SDSSana/gallazzi1/match_gallazzi1.sav'
	sdssem = where(sdss.haveemlines eq 1,csdssem)
	sdss_em = sdss(sdssem)
	noem = where(sdss.haveemlines eq 0, cnoem) ;no emission lines
    	sdss_ageform = (galage(sdss.zfit,1000.)/1.e9-sdss.age)>0. ;age of universe when it was formed	

	;fix the uncertainties
	get_sdss_additional_uncertainties,sdss.feh,sdss.sn,deltaage,deltafeh
	sdss.feherr = deltafeh; sqrt(sci(wearlymem).feherr^2+deltafeh^2)
	ageupper = 10^deltaage*sdss.age
	agelower = sdss.age/10^deltaage
	sdss.ageerr = 0.5*(ageupper-agelower)

	;Catalog from Gallazzi
	restore,'/scr2/nichal/workspace2/SDSSana/gallazzi_allmass2/match_gallazzi_allmass.sav'
	catsdss = cat
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Average the SDSS data
	massbins_sdss = [9,9.8,10.3,10.5,10.7,10.9,11.2,11.5]
	massbins_sdss = [9,9.4,9.8,10.1,10.4,10.7,10.9,11.2,11.5]
	nbins_sdss = n_elements(massbins_sdss)-1
	ave_mass_sdss    = fltarr(nbins_sdss)
	ave_mass_dev_sdss= fltarr(nbins_sdss)
	ave_feh_sdss     = fltarr(nbins_sdss)
	ave_feh_dev_sdss = fltarr(nbins_sdss)
	
	for i=0,nbins_sdss-1 do begin
		msel = where(sdss.logmstar gt massbins_sdss[i] and sdss.logmstar lt massbins_sdss[i+1],cmsel)
		meanerr,sdss(msel).feh,sdss(msel).feherr,fehmean,sigmam,sigmad,sigmas
		ave_feh_sdss(i) = fehmean
		ave_feh_dev_sdss(i) = sigmas
		ave_mass_sdss(i) = mean(sdss(msel).logmstar)
		ave_mass_dev_sdss(i) = stdev(sdss(msel).logmstar)
		;if i eq nbins_sdss-1 then stop
	endfor
	bndry_mass_sdss = massbins_sdss
	hifeh_sdss = interpol(ave_feh_sdss+ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
	lofeh_sdss = interpol(ave_feh_sdss-ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
	mtoohi = where(bndry_mass_sdss gt 11.5,cmtoohi)
	if cmtoohi gt 0 then remove, mtoohi, bndry_mass_sdss, hifeh_sdss,lofeh_sdss

	;Average the SDSS catalog data (measurements from G05)
	ave_mass_catsdss    = fltarr(nbins_sdss)
	ave_mass_dev_catsdss= fltarr(nbins_sdss)
	ave_feh_catsdss     = fltarr(nbins_sdss)
	ave_feh_dev_catsdss = fltarr(nbins_sdss)
	
	for i=0,nbins_sdss-1 do begin
		msel = where(catsdss.logm50 gt massbins_sdss[i] and catsdss.logm50 lt massbins_sdss[i+1],cmsel)
		meanerr,catsdss(msel).z50,(catsdss(msel).z84-catsdss(msel).z16)/2.,fehmean,sigmam,sigmad,sigmas
		ave_feh_catsdss(i) = fehmean
		ave_feh_dev_catsdss(i) = sigmas
		ave_mass_catsdss(i) = mean(catsdss(msel).logm50)
		ave_mass_dev_catsdss(i) = stdev(catsdss(msel).logm50)
		;if i eq nbins_sdss-1 then stop
	endfor
	bndry_mass_catsdss = massbins_sdss
	hifeh_catsdss = interpol(ave_feh_catsdss+ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
	lofeh_catsdss = interpol(ave_feh_catsdss-ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
	mtoohi = where(bndry_mass_catsdss gt 11.5,cmtoohi)
	if cmtoohi gt 0 then remove, mtoohi, bndry_mass_catsdss, hifeh_catsdss,lofeh_catsdss
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Get the average values measured in Gallazzi05
	g05_mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first time is actually 8.91
	g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
	;below are the stars in the figure 8 of Gallazzi05
	g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
	g05_fehlo = g05_feh-g05_feherr
	g05_fehhi = g05_feh+g05_feherr

	;below are the diamonds in figure 8 of Gallazzi05
	;g05_fehlo= [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04]
	;g05_fehhi= [0.0,0.0,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28]
	toolow = where(g05_fehlo lt -1.,ctoolow)
	if ctoolow gt 0 then g05_fehlo(toolow) = -1

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Choi's Data
	z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
	z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
	z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
	z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
	z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Xiangcheng's FIRE data
	readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0pt8.txt',xma_mass08,xma_feh08
	readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0.txt',xma_mass0,xma_feh0
	xma_feh08 = xma_feh08-0.2
	xma_feh0 = xma_feh0-0.2
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;De Rossi 2017 (EAGLE)
	;read off from Figure 5
	DeRossi_z0={mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
	DeRossi_z1={mass:[9.18,9.47,9.85,10.15,10.55],feh:[-0.39,-0.28,-0.09,0.08,0.17],feherr:[0.06,0.06,0.08,0.07,0.085]}
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Sybilska et al 2017 (hELENa, IFU from Sauron data)
	aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
	Syb_z0 = {mass:aa.field1,feh:aa.field2}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;Find the best fit MZR parameters

	;1) functional form from Zahid13
	a = [8.8,9.7,0.6] ;fit parameters
	weights= [1./sci(wearlymem).feherr^2,1./sdss.feherr^2]
	bestfit = curvefit([sci(wearlymem).logmstar,sdss.logmstar],[sci(wearlymem).feh,sdss.feh],weights,a,sigma,function_name='mzrcurve',status=status,chisq=chizahid)
        print,'fitting status: ',status
	zahid_mzr = a

	bestfit_sdss=curvefit(sdss.logmstar,sdss.feh,1./sdss.feherr^2,a,sigma,function_name='mzrcurve',status=status,chisq=chizahid_sdss)
	zahid_mzr_sdss = a
	;2) 1 degree polynomial
	lin_mzr = linfit([sci(wearlymem).logmstar,sdss.logmstar],[sci(wearlymem).feh,sdss.feh],measure_errors=[sci(wearlymem).feherr,sdss.feherr],sigma=lin_mzr_err,chisq=chilin)
	lin_mzr_sdss = linfit(sdss.logmstar,sdss.feh,measure_errors=sdss.feherr,sigma=lin_mzr_sdss_err,chisq=chilin_sdss)
	lin_mzr_cl = linfit(sci(wearlymem).logmstar,sci(wearlymem).feh,measure_errors=sci(wearlymem).feherr,sigma=lin_mzr_cl_err,chisq=chilin_cl)
	print,'best linear fit',lin_mzr
	;3) 2 degree polynomial
	poly_mzr = poly_fit([sci(wearlymem).logmstar,sdss.logmstar],[sci(wearlymem).feh,sdss.feh],2,measure_errors=[sci(wearlymem).feherr,sdss.feherr],sigma=poly_mzr_err,chisq=chipoly)
	poly_mzr_sdss = poly_fit(sdss.logmstar,sdss.feh,2,measure_errors=sdss.feherr,sigma=poly_mzr_sdss_err,chisq=chipoly_sdss)
	poly_mzr_cl = poly_fit(sci(wearlymem).logmstar,sci(wearlymem).feh,2,measure_errors=sci(wearlymem).feherr,sigma=poly_mzr_cl_err,chisq=chipoly_cl)
	;4) linear fit with intrinsic scatter
	bestparam = intrinsic_scatter_linear([sci(wearlymem).logmstar,sdss.logmstar],[sci(wearlymem).feh,sdss.feh],[sci(wearlymem).feherr,sdss.feherr],[lin_mzr,0.1])
	intrinsic_scatter = bestparam[2]
	chilin_intrinsic = bestparam[3] 
	print, 'intrinsic_scatter = ',bestparam[2]
	print, 'best linear fit with intrinsic scatter', bestparam[0:1]

	bestparam_sdss = intrinsic_scatter_linear(sdss.logmstar,sdss.feh,sdss.feherr,[lin_mzr_sdss,0.1])
	bestparam_cl = intrinsic_scatter_linear(sci(wearlymem).logmstar,sci(wearlymem).feh,sci(wearlymem).feherr,[lin_mzr_cl,0.1])
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Use Akaike's information criterion
	samplesize = n_elements([sci(wearlymem).logmstar,sdss.logmstar])
	k=2
	aic_lin = chilin+2.*k+(2.*k*(k+1.))/(samplesize-k-1) 
	k=3
	aic_lin_intrinsic = chilin_intrinsic+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
	aic_poly= chipoly+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
	aic_zahid = chizahid*(samplesize-3)+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
	print, 'AIC linear = ', aic_lin
	print, 'AIC linear with intrinsic =',aic_lin_intrinsic
	print, 'AIC 2degre polynomial = ', aic_poly
	print, 'AIC Zahid =',aic_zahid
	
	print,'Cl0024:linear, polynomial'
	samplesize=n_Elements(sci(wearlymem).logmstar)
	k=2
	print,chilin_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
	k=3
	print,chipoly_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)

	print,'SDSS:linear, polynomial'
	samplesize=n_Elements(sdss)
        k=2
        print,chilin_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
        k=3
        print,chipoly_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)

	; ok, choose linear fit as the best fit according to AIC	
	mzr_mass = findgen(104)/40.+8.9
	mzr_zahid= zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9
	mzr_poly = lin_mzr[0]+lin_mzr[1]*mzr_mass
	mzr_poly_sdss = lin_mzr_sdss[0]+lin_mzr_sdss[1]*mzr_mass
	mzr_poly_cl = lin_mzr_cl[0]+lin_mzr_cl[1]*mzr_mass
;stop
	save, zahid_mzr,poly_mzr,poly_mzr_err,lin_mzr,lin_mzr_err,filename='mzr_obs_param.sav'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;/////////////////////////////////////////////////////////////////////////////////////////////////////////////////;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	stop
	;PLOTTING

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	psname='SDSS_FeH_mass.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[9,11.5]
		yrange=[-1.,0.3]
		plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
	
		;shade the average region
		rainbow_colors,n_colors=21
		polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=30
	
		;x=[bndry_mass,reverse(bndry_mass)]
		;y=[hifeh,reverse(lofeh)]
		;polyfill,x,y,color=fsc_color('org2')
	
		x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
		y=[hifeh_sdss,reverse(lofeh_sdss)]
		x(where(x eq 8.9)) = 9.
		polyfill,x,y,color=60,/line_fill,orientation=45
		polyfill,x,y,color=60,/line_fill,orientation=135
	
		;draw axis
		axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
		axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
		axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
		axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
	
		;draw data points
		;oploterror,sdss.logmstar,sdss.feh,sdss.feherr,psym=1,errcolor=40,color=40,thick=2
		cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
		;oploterror,sdss_em.logmstar,sdss_em.feh,sdss_em.feherr,psym=1,$
		;errcolor=40,color=40,thick=2
		cgplot,sdss_em.logmstar,sdss_em.feh,psym=14,/overplot,color=40,symsize=1.3
		;plot representative uncertainties
		lowermass_pop= where(sdss.logmstar lt 10.,complement=uppermass_pop)
		lowermass_em_pop = where(sdss_em.logmstar lt 10.,complement=uppermass_em_pop)
		lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
	        oploterror,[11.1,11.3],[-0.35,-0.35],[median(0.5*catsdss(lowermass).logm84-0.5*catsdss(lowermass).logm16),median(0.5*catsdss(uppermass).logm84-0.5*catsdss(uppermass).logm16)],[median([sdss(lowermass_pop).feherr,sdss_em(lowermass_em_pop).feherr]),median([sdss(uppermass_pop).feherr,sdss_em(uppermass_em_pop).feherr])],color=40,psym=3,errcolor=40,errthick=1.5
	
		;Add DeRossi17 
		oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=90,$
		linethick=2,errcolor=90,psym=1
		oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),symsize=1.5,color=90
	
		;Add Xiangcheng's data
		oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=90,symsize=1.2

                ;Add Choi's data
                oploterror,z01.mass,z01.feh,z01.feherr,color=10,linethick=2,errcolor=10
                oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=10,symsize=2
	
		;Add Sybilska2017
		oplot,syb_z0.mass,syb_z0.feh,color=10,thick=2,linestyle=2

		;Add best fitted line
		oplot,mzr_mass,mzr_poly_sdss,color=0,thick=2
		oplot,mzr_mass,zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9,color=0,linestyle=3,thick=2
		oplot,mzr_mass,poly_mzr_sdss[0]+poly_mzr_sdss[1]*mzr_mass+poly_mzr_sdss[2]*mzr_mass^2,color=0,linestyle=3,thick=2
                ;Labelling
		cglegend,title=['Early-type SDSS','measured in this work','Gallazzi et al. 2005','Choi et al. 2014','Sybilska et al. 2017','Ma et al. 2016','De Rossi et al. 2017'],psym=[14,0,15,46,0,16,24],location=[10.6,-0.5],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=[40,40,30,10,10,90,90],symsize=1.2
		oplot,[10.53,10.67],[-0.8,-0.8],linestyle=2,color=10,thick=2
		xyouts,10.1,-0.52,'observations:',charsize=0.8
		xyouts,10.1,-0.88,'simulations:',charsize=0.8
		xyouts,9.1,0.2,'z~0 galaxies',charsize=1.
		;bracket,10.45,-0.8,10.5,-0.5,/left
		;bracket,10.45,-0.95,10.5,-0.85,/left
	device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	psname='SDSS_FeH_mass_G05measure.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[9,11.5]
		yrange=[-1.,0.3]
		plot,catsdss.logm50,catsdss.z50,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
	
		;shade the average region
		polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('powderblue')
	
		x=[bndry_mass_catsdss,reverse(bndry_mass_catsdss)]
		y=[hifeh_catsdss,reverse(lofeh_catsdss)]
		x(where(x eq 8.9)) = 9.
		polyfill,x,y,color=fsc_color('dodgerblue'),/line_fill,orientation=45
		polyfill,x,y,color=fsc_color('dodgerblue'),/line_fill,orientation=135
	
		;draw axis
		axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
		axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
		axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
		axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
	
		;draw data points
		cgerrplot,catsdss.logm50,catsdss.z16,catsdss.z84,color='ygb5',thick=2
		cgplot,catsdss.logm50,catsdss.z50,psym=14,/overplot,color=fsc_color('ygb5'),symsize=1.3
		cgerrplot,catsdss(sdssem).logm50,catsdss(sdssem).z16,catsdss(sdssem).z84,color='darkgray',thick=2
		cgplot,catsdss(sdssem).logm50,catsdss(sdssem).z50,psym=14,/overplot,color=fsc_color('darkgray'),symsize=1.3
		
		;Add DeRossi17 
		oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=fsc_color('ygb7'),$
		linethick=2,errcolor=fsc_color('ygb7')
		oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),symsize=1.5,color=fsc_color('ygb7')
	
		;Add Xiangcheng's data
		oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=fsc_color('ygb7'),symsize=1.5
	
		;Add Sybilska2017
		oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('blue'),thick=2
	device,/close
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	cleanplot,/silent
	!p.font=0
	psname='Cl0024_FeH_mass.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[8.9,11.5]
		yrange=[-1.,0.3]
		plot,sci(wearlymem).logmstar,sci(wearlymem).feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5

		;shade the average region
		;colorarr=[10,50,203,230,254]
		colorarr= fsc_color(['darkorchid','royalblue','org5','maroon','red8'])
	 	x=[bndry_mass,reverse(bndry_mass)]
	        y=[hifeh,reverse(lofeh)]
        	polyfill,x,y,color=fsc_color('org2')

		x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
                y=[hifeh_sdss,reverse(lofeh_sdss)]
                x(where(x eq 9)) = 8.9
                polyfill,x,y,color=fsc_color('blu6'),/line_fill,orientation=45
                polyfill,x,y,color=fsc_color('blu6'),/line_fill,orientation=135
		
		;draw axis
		axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
		axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
		axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
		axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

		;Add Xiangcheng's data
		oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=colorarr(0)
		oplot,xma_mass08,xma_feh08,psym=cgsymcat(16),color=colorarr(4)
		
		;Add Munoz15

		;Add DeRossi17 
		oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0)
		oploterror,derossi_z1.mass,derossi_z1.feh,derossi_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)      
		oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.1
		oplot,derossi_z1.mass,derossi_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.1

		;draw data points
		oploterror,sci(wearlymem).logmstar,sci(wearlymem).feh,sci(wearlymem).feherr,psym=3,errcolor=fsc_color('org4'),color=colorarr(2)
		oplot,sci(wearlymem).logmstar,sci(wearlymem).feh,psym=cgsymcat(14),color=colorarr(2),symsize=1.3
		oplot,sci(wearlymem).logmstar,sci(wearlymem).feh,psym=cgsymcat(4),color=fsc_color('org7'),symsize=1.3

		;Add Choi's data
		vsym,5,/fill,/star
		oploterror,z01.mass,z01.feh,z01.feherr,color=colorarr(1),linethick=2,errcolor=colorarr(1)
		oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color('org7'),linethick=2,errcolor=fsc_color('org7')
		oploterror,z06.mass,z06.feh,z06.feherr,color=colorarr(3),linethick=2,errcolor=colorarr(3)
		oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=colorarr(1),symsize=2
		oplot,z03.mass,z03.feh,psym=cgsymcat(46),color=colorarr(2),symsize=2
		oplot,z06.mass,z06.feh,psym=cgsymcat(46),color=colorarr(3),symsize=2
		oplot,z01.mass,z01.feh,psym=cgsymcat(45),symsize=2
		oplot,z03.mass,z03.feh,psym=cgsymcat(45),symsize=2
		oplot,z06.mass,z06.feh,psym=cgsymcat(45),symsize=2

		;add best fitted MZR curve
		;oplot,mzr_mass,mzr_zahid,linestyle=2
		oplot,mzr_mass,mzr_poly,thick=2
		oplot,mzr_mass,mzr_poly_sdss,thick=2,linestyle=2
		oplot,mzr_mass,mzr_poly_cl,thick=2,linestyle=2
		;Labelling
		al_legend,['z=0','z=[0.1,0.2]','z=[0.3,0.4]','z=[0.6,0.7]','z=[0.8,1]'],psym=15,color=colorarr,box=0,position=[10.6,-0.6]
		al_legend,['Current work','Choi et al. 2014','Ma et al. 2016','De Rossi et al. 2017'],color='black',psym=[14,46,16,24],position=[10.6,-0.32],box=0,charsize=0.9,symsize=[1.3,1.7,1.,1.1]		
		xyouts,9.,0.2,'z~0.4 galaxies',charsize=1.
	device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Gas phase metal from Zahid2013 (12+log(O/H))
	gas_oh = [{redshift:0.08,z0:9.121,M0:8.999,gmma:0.85},{redshift:0.29,z0:9.130,M0:9.304,gmma:0.77},{redshift:0.78,z0:9.161,M0:9.661,gmma:0.65},{redshift:1.4,z0:9.06,M0:9.6,gmma:0.7},{redshift:2.26,z0:9.06,M0:9.7,gmma:0.6}]
	;the equation is 12+log(O/H) = z0-log[1+(M*/M0)^-gmma]
	sun_oh = 8.9
	;alpha_Fe = 0.12*mass-1.1 ;linear plot by eye to Choi14 Fig 8, mass is in log scale
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	psname='Formation_Redshift_cl0024.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		zcat = [1000,2,1,0.7,0.4,0.]
	;	zcat = [1000.,1.5,0.7,0]
		agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif	
		plot,sci(wearlymem).logmstar,sci(wearlymem).feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2],/nodata
		
		nage = n_Elements(agecat)-1
		rainbow_colors
		zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
		;loop over formation times
		for i=0,nage-1 do begin
			selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
			if cselsdss gt 0 then begin
				cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.7
			endif
		endfor
        	for i=0,nage-1 do begin
            		sel = where(earlymem eq 1 and ageform gt agecat(i) and ageform le agecat(i+1), csel)
            		if csel gt 0 then begin
				oploterror,sci(sel).logmstar,sci(sel).feh,sci(sel).feherr,psym=1,errcolor=zcolor(i)-10,thick=0.5
				cgplot,sci(sel).logmstar,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
				cgplot,sci(sel).logmstar,sci(sel).feh,psym=4,/overplot,color='darkgray',symsize=1.3
            		endif
        	endfor
		;add gas phase MZR
		marr = findgen(101)/40.+9. ;logmass from 9 to 11.5
		a_fe = 0.12*marr-1.1
		for i=0,n_Elements(gas_oh)-1 do begin
			o_h = gas_oh(i).z0-alog10(1.+10.^((marr-gas_oh(i).m0)*(-1.)*gas_oh(i).gmma))-sun_oh;-a_fe
			colornowi = value_locate(zcat,gas_oh(i).redshift)
			colornow = zcolor(colornowi)
			oplot,marr,o_h,color=colornow
		endfor
		;Labelling
		zarr_str = strarr(n_elements(zcat)-1)
		for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z$\tex_{form}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
		zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
		al_Legend,zarr_str,psym=15,color=zcolor,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
        	al_Legend,['SDSS subsample','Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,charsize=1,position=[10.6,-0.3],font=0
		xyouts,9.1,0.1,'(a)',charsize=1.2
	device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        psname='Formation_Redshift_cl0024_extreme.eps'
        device, filename = psname,xsize = 15,ysize = 10, $
                        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
                zcat = [1000,2,0.7,0.4,0.]
                agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif       
                plot,sci(wearlymem).logmstar,sci(wearlymem).feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2],/nodata

                nage = n_Elements(agecat)-1
                rainbow_colors
                zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
                ;loop over formation times
                for i=0,nage-1 do begin
			if i eq 0 then begin
                        	selsdss = where(sdss_ageform+sdss.ageerr gt agecat(i) and sdss_ageform+sdss.ageerr le agecat(i+1),cselsdss)
                        	if cselsdss gt 0 then begin
                			cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh+sdss(selsdss).feherr,psym=16,/overplot,color=zcolor(i),symsize=0.5
                        	endif
			endif
			if i eq 3 then begin
                                selsdss = where(sdss_ageform-sdss.ageerr gt agecat(i) and sdss_ageform-sdss.ageerr le agecat(i+1),cselsdss)
                                if cselsdss gt 0 then begin
                                        cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh-sdss(selsdss).feherr,psym=16,/overplot,color=zcolor(i),symsize=0.5
                                endif
			endif
                endfor
        	for i=0,nage-1 do begin
			if i eq 0 then begin
            			sel = where(earlymem eq 1 and ageform+sci.ageerr gt agecat(i) and ageform+sci.ageerr le agecat(i+1), csel)
				if csel gt 0 then begin
                		;oploterror,sci(sel).logmstar,sci(sel).feh,sci(sel).feherr,psym=1,errcolor=zcolor(i)
               			cgplot,sci(sel).logmstar,sci(sel).feh+sci(sel).feherr,psym=14,/overplot,color=zcolor(i),symsize=1.3
            			endif
			endif
			if i eq 2 then begin
                                sel = where(earlymem eq 1 and ageform-sci.ageerr gt agecat(i) and ageform-sci.ageerr le agecat(i+1), csel)
                                if csel gt 0 then begin
                                ;oploterror,sci(sel).logmstar,sci(sel).feh,sci(sel).feherr,psym=1,errcolor=zcolor(i)
                                cgplot,sci(sel).logmstar,sci(sel).feh-sci(sel).feherr,psym=14,/overplot,color=zcolor(i),symsize=1.3
                                endif
                        endif

        	endfor
                ;Labelling
		xyouts,9.1,0.12,'Extreme limits along age-metallicity degeneracy'
                zarr_str = strarr(n_elements(zcat)-1)
                for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z$\tex_{form,min}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
                zarr_str(0) = 'z$\tex_{form,max}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
                al_Legend,zarr_str([0,2,3]),psym=15,color=zcolor([0,2,3]),box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
        al_Legend,['SDSS subsample','Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,charsize=1,position=[10.6,-0.4],font=0
        device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	psname='FEH_deviation_obs.eps'
        mzrparam=[9.01704,9.91704,0.817038]
	fehideal = mzrparam[0]-alog10(1.+10.^((sci(wearlymem).logmstar-mzrparam[1])*(-1.*mzrparam[2])))-8.9
	sdss_fehideal = mzrparam[0]-alog10(1.+10.^((sdss.logmstar-mzrparam[1])*(-1.*mzrparam[2])))-8.9
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		plot,sdss_ageform,sdss.feh-sdss_fehideal,psym=1,xtitle='Age of Universe at Formation',ytitle=delta+'[Fe/H]',xrange=[0,14],xstyle=9,yrange=[-0.4,0.4],/nodata,position=[0.15,0.15,0.95,0.9]
		axis,xaxis=0,xstyle=1,xrange=[0,14]
		axis,xaxis=1,xticks=4,xtickv=[1.558,3.316, 5.903, 8.628,13.712],xtickn=['4','2','1','0.5','0'],xtitle='z'
		cgplot,sdss_ageform,sdss.feh-sdss_fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
		cgplot,ageform(wearlymem),sci(wearlymem).feh-fehideal,psym=14,/overplot,color='org4'
		x=[sdss_ageform,ageform(wearlymem)]
		y=[sdss.feh-sdss_fehideal,sci(wearlymem).feh-fehideal]
		dx = [sdss.ageerr,sci(wearlymem).ageerr]
		dy = [sdss.feherr,sci(wearlymem).feherr]
		outlier = where(y gt 0.02+x*0.4/8.,coutlier,complement=goodones)
		linparam=linfit(x(goodones),y(goodones))
		fitexy,x(goodones),y(goodones),A,B,x_sig=dx(goodones),y_sig=dy(goodones),sigma_a_b,chi_sq,q
		oplot,[0,14],A+B*[0,14],thick=2
		;oplot,[0,14],linparam[0]+linparam[1]*[0,14],color=cgcolor('black'),thick=2
		print,A,B
		mockparam=[-0.0677188,0.0040122]
		;oplot,[0,14],mockparam[0]+mockparam[1]*[0,14],color=cgcolor('darkgray'),linestyle=2,thick=2
		oploterror,[12.5],[-0.25],[median(dx(goodones))],[median(dy(goodones))]
		xyouts,0.5,0.32,'(c) real observation',charsize=1.2
	device,/close
	;try fitting separately
	;SDSS
	fitexy,sdss_ageform,sdss.feh-sdss_fehideal,Asdss,Bsdss,x_sig=sdss.ageerr,y_sig=sdss.feherr,sigma_a_b_sdss,chisqasdss,qsdss
	x=ageform(wearlymem)
	y=sci(wearlymem).feh-fehideal
	dx = sci(wearlymem).ageerr
	dy = sci(wearlymem).feherr
	;outlier = where(y gt 0.02+x*0.4/8.,coutlier,complement=goodones)
	outlier = where(y gt 0.3,coutlier,complement=goodones)
	fitexy,x(goodones),y(goodones),Acl,Bcl,x_sig=dx(goodones),y_sig=dy(goodones),sigma_a_b_cl,chisqcl,qcl
	stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	psname='Cl0024_age_mass.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[8.9,11.5]
		plot,sci(wearlymem).logmstar,sci(wearlymem).age,xrange=xrange,xstyle=1,psym=cgsymcat(14)
		oplot,sdss.logmstar,sdss.age,psym=cgsymcat(16),symsize=0.5
	device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	psname='mass_alpha.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[8.9,11.5]
		yrange=[-0.05,0.45]
		plot,sci(wearlymem).logmstar,sci(wearlymem).alphafe,/nodata,xrange=xrange,xstyle=1,yrange=yrange,ystyle=1,xtitle='Log(M/M'+sunsym+')',ytitle='['+alpha+'/Fe]'
		;draw data points
		cgplot,sdss.logmstar,sdss.alphafe,psym=16,color='ygb5',/overplot
		cgerrplot,sci(wearlymem).logmstar,sci(wearlymem).alphafelower,sci(wearlymem).alphafeupper,color='org4'
		cgplot,sci(wearlymem).logmstar,sci(wearlymem).alphafe,psym=14,/overplot,color='org4',symsize=1.3
	device,/close

	psname='age_alpha.eps'
	device, filename = psname,xsize = 14,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[0,13]
		massrange=[9,11.5]
		yrange=[-0.05,0.45]
		plot,sci(wearlymem).age,sci(wearlymem).alphafe,/nodata,xrange=xrange,xstyle=1,yrange=yrange,ystyle=1,xtitle='Age(GYR)',ytitle='['+alpha+'/Fe]',position=[0.15,0.15,0.9,0.8]
		;draw data points
		sdssmasscolor=bytscl(sdss.logmstar,max=max(massrange),min=min(massrange))
		masscolor=bytscl(sci(wearlymem).logmstar,max=max(massrange),min=min(massrange))
		burd_colors
		cgscatter2d,sdss.age,sdss.alphafe,color=sdssmasscolor,psym=16,/overplot,/nodisplay,fit=0
		cgerrplot,sci(wearlymem).age,sci(wearlymem).alphafelower,sci(wearlymem).alphafeupper,color='gray'
		cgscatter2d,sci(wearlymem).age,sci(wearlymem).alphafe,color=masscolor,psym=14,/overplot,/nodisplay,fit=0,symsize=1.3
		cgcolorbar,divisions=4,minor=5,format='(F4.1)',range=massrange,title='Log(M/M'+sunsym+')',tlocation='TOP',charsize=0.8

	device,/close
stop
end

