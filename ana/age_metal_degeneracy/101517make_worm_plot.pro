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
		grid_feh = rebin(indgen(120)/120.-0.8,120,126) ;range from [-0.6,0.192] dex
		grid_age = transpose(rebin(alog10(findgen(126)/10.+0.5)+9.,126,120)) ;range from [0.5,13] Gyr
		minchisq = min(chisqarr)
		contour,chisqarr,grid_feh,grid_age,levels=[minchisq+0.02,minchisq+1],c_linestyle=0,/overplot,color=fehcolor(i)
		oplot,[fehref],alog10([ageref])+9.,psym=7,color=fehcolor(i)
		
	endfor
	device,/close
	stop
end
