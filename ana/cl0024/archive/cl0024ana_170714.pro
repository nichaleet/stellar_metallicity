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

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;PLOTTING
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	psname='Cl0024_FeH_mass.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[8.9,11.5]
		yrange=[-1.,0.3]
		plot,sci(wearlymem).logmstar,sci(wearlymem).feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5

		;shade the average region
	 	x=[bndry_mass,reverse(bndry_mass)]
	        y=[hifeh,reverse(lofeh)]
        	polyfill,x,y,color=fsc_color('org2')

		;draw axis
		axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
		axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
		axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
		axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

		;draw data points
		oploterror,sci(wearlymem).logmstar,sci(wearlymem).feh,sci(wearlymem).feherr,psym=1,errcolor=fsc_color('org4'),color=fsc_color('org4'),thick=2
		cgplot,sci(wearlymem).logmstar,sci(wearlymem).feh,psym=14,/overplot,color=fsc_color('org4'),symsize=1.3

		;Add Choi's data
		vsym,5,/fill,/star
		oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('ygb5'),linethick=2,errcolor=fsc_color('ygb5')
		oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color('org4'),linethick=2,errcolor=fsc_color('org4')
		oploterror,z06.mass,z06.feh,z06.feherr,color=fsc_color('maroon'),linethick=2,errcolor=fsc_color('maroon')
		oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('ygb5'),symsize=2
		oplot,z03.mass,z03.feh,psym=cgsymcat(46),color=fsc_color('org4'),symsize=2
		oplot,z06.mass,z06.feh,psym=cgsymcat(46),color=fsc_color('maroon'),symsize=2

		;Add Xiangcheng's data
		oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=fsc_color('ygb3')
		oplot,xma_mass08,xma_feh08,psym=cgsymcat(16),color=fsc_color('red8')
		
		;Add Munoz15

		;Add DeRossi17 
		oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=fsc_color('ygb3'),linethick=2,errcolor=fsc_color('ygb3')
		oploterror,derossi_z1.mass,derossi_z1.feh,derossi_z1.feherr,color=fsc_color('red8'),linethick=2,errcolor=fsc_color('red8')      
		oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),color=fsc_color('ygb3')
		oplot,derossi_z1.mass,derossi_z1.feh,psym=cgsymcat(24),color=fsc_color('red8')	
		;Labelling
		al_legend,['z=0','z=[0.1,0.2]','z=[0.3,0.4]','z=[0.6,0.7]','z=[0.8,1]'],psym=15,color=['ygb3','ygb5','org4','maroon','red8'],box=0,position=[10.1,-0.6]
		al_legend,['Current work','Choi et al. 2014','Ma et al. 2016','De Rossi et al. 2017'],color='black',psym=[14,46,16,24],position=[10.7,-0.6],box=0,charsize=0.9		
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
		ploterror,sci(wearlymem).logmstar,sci(wearlymem).feh,sci(wearlymem).feherr,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2],/nodata
		
		nage = n_Elements(agecat)-1
		rainbow_colors
		zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
		;loop over formation times
		for i=0,nage-1 do begin
			selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
			if cselsdss gt 0 then begin
            	cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.5
			endif
		endfor
        	for i=0,nage-1 do begin
            		sel = where(earlymem eq 1 and ageform gt agecat(i) and ageform le agecat(i+1), csel)
            		if csel gt 0 then begin
                	oploterror,sci(sel).logmstar,sci(sel).feh,sci(sel).feherr,psym=1,errcolor=zcolor(i)
               		cgplot,sci(sel).logmstar,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
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
	psname='SDSS_FeH_mass.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		;make the outline of plots
		xrange=[9,11.5]
		yrange=[-1.,0.3]
		plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
	
		;shade the average region
		polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('powderblue')
	
		;x=[bndry_mass,reverse(bndry_mass)]
		;y=[hifeh,reverse(lofeh)]
		;polyfill,x,y,color=fsc_color('org2')
	
		x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
		y=[hifeh_sdss,reverse(lofeh_sdss)]
		x(where(x eq 8.9)) = 9.
		polyfill,x,y,color=fsc_color('dodgerblue'),/line_fill,orientation=45
		polyfill,x,y,color=fsc_color('dodgerblue'),/line_fill,orientation=135
	
		;draw axis
		axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
		axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
		axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
		axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
	
		;draw data points
		oploterror,sdss.logmstar,sdss.feh,sdss.feherr,psym=1,errcolor=fsc_color('ygb5'),color=fsc_color('ygb5'),thick=2
		cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=fsc_color('ygb5'),symsize=1.3
		oploterror,sdss_em.logmstar,sdss_em.feh,sdss_em.feherr,psym=1,$
		errcolor=fsc_color('darkgray'),color=fsc_color('darkgray'),thick=2
		cgplot,sdss_em.logmstar,sdss_em.feh,psym=14,/overplot,color=fsc_color('darkgray'),symsize=1.3
		
		;Add DeRossi17 
		oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=fsc_color('ygb7'),$
		linethick=2,errcolor=fsc_color('ygb7')
		oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),symsize=1.5,color=fsc_color('ygb7')
	
		;Add Xiangcheng's data
		oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=fsc_color('ygb7'),symsize=1.5

                ;Add Choi's data
                vsym,5,/fill,/star
                oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('ygb5'),linethick=2,errcolor=fsc_color('ygb5')
                oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('ygb5'),symsize=2
	
		;Add Sybilska2017
		oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('blue'),thick=2
	device,/close
	
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

stop
end

