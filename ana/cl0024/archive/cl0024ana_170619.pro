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
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;De Rossi 2017 (EAGLE)
	;read off from Figure 5
	DeRossi_z0={mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
	DeRossi_z1={mass:[9.18,9.47,9.85,10.15,10.55],feh:[-0.39,-0.28,-0.09,0.08,0.17],feherr:[0.06,0.06,0.08,0.07,0.085]}
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;PLOTTING
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
        ;Get SDSS data
        sdss = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi1/sps_fit.fits.gz',1,/silent)
	;remove bad measurement
  	nofit = where(sdss.feh eq -999 ,cnofit)
  	if cnofit gt 0 then remove,nofit,sdss
	sdss(where(sdss.feherr eq 0.)).feherr = median(sdss.feherr)
;	restore,'/scr2/nichal/workspace2/SDSSana/gallazzi1/match_gallazzi1.sav'
        noem = where(sdss.haveemlines eq 0, cnoem) ;no emission lines
        sdss = sdss(noem)  ;only used non star-forming galaxies

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Average the SDSS data
	massbins_sdss = [9.8,10.3,10.5,10.7,10.9,11.2,11.5,12.]
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
        endfor
	bndry_mass_sdss = massbins_sdss
        hifeh_sdss = interpol(ave_feh_sdss+ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
        lofeh_sdss = interpol(ave_feh_sdss-ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
	mtoohi = where(bndry_mass_sdss gt 11.5,cmtoohi)
	if cmtoohi gt 0 then remove, mtoohi, bndry_mass_sdss, hifeh_sdss,lofeh_sdss
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Get the average values measured in Gallazzi05
	g05_mass = [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5]
	g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
	;below are the stars in the figure 8 of Gallazzi05
	g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
	g05_fehlo = g05_feh-g05_feherr
        g05_fehhi = g05_feh+g05_feherr

	;below are the diamonds in figure 8 of Gallazzi05
;	g05_fehlo= [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04]
;	g05_fehhi= [0.0,0.0,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28]
	toolow = where(g05_fehlo lt -1.,ctoolow)
	if ctoolow gt 0 then g05_fehlo(toolow) = -1
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        psname='SDSS_FeH_mass.eps'
        device, filename = psname,xsize = 15,ysize = 10, $
                xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
                ;make the outline of plots
                xrange=[8.9,11.5]
                yrange=[-1.,0.3]
                plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5

                ;shade the average region
		polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('ygb3')

                x=[bndry_mass,reverse(bndry_mass)]
                y=[hifeh,reverse(lofeh)]
                polyfill,x,y,color=fsc_color('org2')

		x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
                y=[hifeh_sdss,reverse(lofeh_sdss)]
		polyfill,x,y,color=fsc_color('blue'),/line_fill,orientation=45
		polyfill,x,y,color=fsc_color('blue'),/line_fill,orientation=135

               ;draw axis
                axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
                axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
                axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
                axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'

                ;draw data points
                oploterror,sdss.logmstar,sdss.feh,sdss.feherr,psym=1,errcolor=fsc_color('ygb5'),color=fsc_color('ygb5'),thick=2
                cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=fsc_color('ygb5'),symsize=1.3

	device,/close
stop
end

