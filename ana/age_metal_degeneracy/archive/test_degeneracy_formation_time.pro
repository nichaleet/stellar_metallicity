pro test_degeneracy_formation_time
;this program is to test whether age-degeneracy can cause the scatter in the observed MZR
;and create a separation into formation time.
;It first assumes that galaxies follow a tight correlation of MZR from the best-fit MZR in
;cl0024ana.pro that is [z0,M0,gmma] = [9.01704,9.91704,0.817038]
;Then galaxies adopt random values within chisq=chisq+min(chisq)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	set_plot,'ps'
	!p.font =  0
;open or make the mock set of galaxies
	if file_test('mock_gal_degeneracy.fits') eq 0 then make_mock_gal_degeneracy 
	sdss = mrdfits('mock_gal_degeneracy.fits',1)
	if file_test('mock_cl0024_degeneracy.fits') eq 0 then make_mock_cl0024_degeneracy
	sci = mrdfits('mock_cl0024_degeneracy.fits',1)
	
;= {sn:0.,mass:0.,fehideal:0.,ageideal:0.,feh:0.,feherr:0.,age:0.,ageerr:0.}
 	sunsym = sunsymbol()
	Delta = '!9'+string("104B)+'!x'
	alpha = '!9'+string("141B)+'!x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Gas phase metal from Zahid2013 (12+log(O/H))
	gas_oh = [{redshift:0.08,z0:9.121,M0:8.999,gmma:0.85},{redshift:0.29,z0:9.130,M0:9.304,gmma:0.77},{redshift:0.78,z0:9.161,M0:9.661,gmma:0.65},{redshift:1.4,z0:9.06,M0:9.6,gmma:0.7},{redshift:2.26,z0:9.06,M0:9.7,gmma:0.6}]
	;the equation is 12+log(O/H) = z0-log[1+(M*/M0)^-gmma]
	sun_oh = 8.9
	;alpha_Fe = 0.12*mass-1.1 ;linear plot by eye to Choi14 Fig 8, mass is in log scale
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plotting
	nsdssgals = n_elements(sdss)
	obs_z = fltarr(nsdssgals)+0.05 ;observation redshift (check with make_mock_gal_degeneray.pro)
	sdss_ageform = (galage(0.05,1000.)/1.e9-sdss.age)>0.
	ageform = (galage(sci.z,1000.)/1.e9-sci.age)>0.

	psname='Formation_Redshift_mock.eps'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		zcat = [1000,2,1,0.7,0.4,0.]
		agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif	
		plot,sdss.mass,sdss.feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2],/nodata
		
		nage = n_Elements(agecat)-1
		rainbow_colors
		zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
		;loop over formation times
		for i=0,nage-1 do begin
			selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
			if cselsdss gt 0 then begin
            			cgplot,sdss(selsdss).mass,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.5
			endif
		endfor

		for i=0,nage-1 do begin
			sel = where(ageform gt agecat(i) and ageform le agecat(i+1), csel)
			if csel gt 0 then begin
				oploterror,sci(sel).mass,sci(sel).feh,sci(sel).feherr,psym=1,errcolor=zcolor(i)
				cgplot,sci(sel).mass,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
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
        al_Legend,['mock SDSS','mock Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,charsize=1,position=[10.6,-0.3],font=0
	device,/close

	psname='FEH_deviation_mock.eps'
	Delta = '!9'+string("104B)+'!x'
	Delta = '!9'+string("104B)+'!x'
	device, filename = psname,xsize = 15,ysize = 10, $
			xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
		plot,sdss_ageform,sdss.feh-sdss.fehideal,psym=1,xtitle='Age of Universe at Formation',ytitle=delta+'[Fe/H]',xrange=[0,14],xstyle=9,yrange=[-0.4,0.4],/nodata,position=[0.15,0.15,0.95,0.9]
		axis,xaxis=0,xstyle=1,xrange=[0,14]
		axis,xaxis=1,xticks=4,xtickv=[1.558,3.316, 5.903, 8.628,13.712],xtickn=['4','2','1','0.5','0'],xtitle='z'
		cgplot,sdss_ageform,sdss.feh-sdss.fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
		cgplot,ageform,sci.feh-sci.fehideal,psym=14,/overplot,color='org4'
		badageerr = where(sdss.ageerr eq 0.)
		sdss(badageerr).ageerr = median(sdss.ageerr)
		badageerr = where(sci.ageerr eq 0.)
		sci(badageerr).ageerr = median(sci.ageerr)
		linparam=linfit([sdss_ageform,ageform],[sdss.feh-sdss.fehideal,sci.feh-sci.fehideal])
		print,linparam
		oplot,[0,14],linparam[0]+linparam[1]*[0,14],color=cgcolor('black'),thick=2
		obsparam=[-0.13512675,0.020072897]
		oplot,[0,14],obsparam[0]+obsparam[1]*[0,14],color=cgcolor('darkgray'),linestyle=2,thick=2
	device,/close		

end 



