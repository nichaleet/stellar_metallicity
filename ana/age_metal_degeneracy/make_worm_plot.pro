pro make_worm_plot,dir

	fehplot = [-0.44,-0.20,-0.12,0.03]
	ageplot = [2.5,3.5,4.5,5.5,6.5,7.5]
	nfeh = n_Elements(fehplot)
	
	set_plot,'ps'
	psname=dir+'/chisq_worm_plot_'+dir+'.eps'
	!p.font=0
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	plot,findgen(100)/125.-0.6,alog10(findgen(126)/10.+0.5)+9.,xrange=[-0.6,0.2],yrange=alog10([2,12.])+9,/nodata,$
		xtitle='[Fe/H]',xstyle=1,ystyle=5;,position=[0.15,0.1,0.85,0.9]
	axis, yaxis=0,ytitle='log(Age/yr)',yrange=alog10([2,12.])+9,ystyle=1
	axis, yaxis=1,yticks=6,ytickv=alog10([2,3,4,5,6,8,11.])+9.,ytickn=['2','3','4','5','6','8','11'],ycharsize=0.7
	burd_colors
	fehcolor=[10,60,190,254]
	for i=0,nfeh-1 do for j=0,n_elements(ageplot)-1 do begin
		fehref = fehplot(i)
		ageref = ageplot(j)
		;read data
		fitsfile = dir+"/chisq_age"+strtrim(string(fix(ageref*10)),2)+"_feh"+strtrim(string(fix(fehref*100)),2)+".fits"	
		if file_test(fitsfile) eq 0 then stop,'no fits file found'
		chisqarr = mrdfits(fitsfile,0,header)
		;you should check below if it matches those in degeneracy.pro
		grid_feh = rebin(findgen(100)/125.-0.6,100,126) ;range from [-0.6,0.192] dex
		grid_age = transpose(rebin(findgen(126)/10.+0.5,126,100)) ;range from [0.5,13] Gyr
		feharr = findgen(100)/125.-0.6
                agearr = findgen(126)/10.+0.5

		Lgrid = -0.5*(chisqarr-min(chisqarr))
                probgrid = exp(double(Lgrid))
                volume = 0.
		deltaage = agearr(1)-agearr(0)
		deltafeh = feharr(1)-feharr(0)
                cube_dimen = size(chisqarr,/dimensions)
                for ii=0,cube_dimen(0)-1 do begin
                   for jj=0,cube_dimen(1)-1 do begin
                     volume = volume+deltafeh*deltaage*probgrid[ii,jj]
                   endfor
                endfor
                probgrid = probgrid/volume

	        probfeh = fltarr(cube_dimen(0))
	        probage = fltarr(cube_dimen(1))
       		for ii=0,cube_dimen(0)-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
       		for jj=0,cube_dimen(1)-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
       		probfeh = probfeh/int_tabulated(feharr,probfeh)
       		probage = probage/int_tabulated(agearr,probage)

       		cumprobfeh = fltarr(cube_dimen(0))
       		for jj=1,cube_dimen(0)-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
       		cumprobage = fltarr(cube_dimen(1))
       		for jj=1,cube_dimen(1)-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])
                
                midage = interpol(agearr,cumprobage,0.5)
                midfeh = interpol(feharr,cumprobfeh,0.5)
                maxprob = max(probgrid,wmaxprob)
		maxprobage =grid_age(wmaxprob)
		maxprobfeh =grid_feh(wmaxprob)
                ;shift midage and midfeh to the center of the fehref and ageref
;		shiftfeh = fehref-midfeh
		shiftfeh = fehref-maxprobfeh
;		shiftage = ageref-midage
		shiftage = ageref-maxprobage
		contour,probgrid,grid_feh+shiftfeh,alog10(grid_age+shiftage)+9.,levels=[maxprob/1.65],c_linestyle=0,/overplot,color=fehcolor(i)
		oplot,[fehref],alog10([ageref])+9.,psym=7,color=fehcolor(i)
		
	endfor
	device,/close
	stop
end
