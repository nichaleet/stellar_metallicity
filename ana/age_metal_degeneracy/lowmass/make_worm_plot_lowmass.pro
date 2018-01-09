pro make_worm_plot_lowmass,dir

	fehplot = [-1.09,-0.820,-0.553,-0.287]
	ageplot = [1.5,2.83,4.17,5.50,6.83]
	nfeh = n_Elements(fehplot)

	;you should check below if it matches those in degeneracy_lowmass.pro
	grid_feh = rebin(indgen(126)/75.-1.5,126,121)
	grid_age = transpose(rebin(alog10(findgen(121)/15.+0.1)+9.,121,126)) 
	
	set_plot,'ps'
	psname=dir+'/chisq_worm_plot_'+dir+'.eps'
	!p.font=0
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	plot,findgen(126)/75.-1.5,alog10(findgen(121)/15.+0.1)+9.,xrange=[-1.5,0.2],yrange=alog10([1,8.])+9,/nodata,$
		xtitle='[Fe/H]',xstyle=1,ystyle=5;,position=[0.15,0.1,0.85,0.9]
	axis, yaxis=0,ytitle='log(Age/yr)',yrange=alog10([1,8.])+9,ystyle=1
	axis, yaxis=1,yticks=6,ytickv=alog10([1,2,3,4,5,6,8])+9.,ytickn=['1','2','3','4','5','6','8'],ycharsize=0.7
	burd_colors
	fehcolor=[10,60,190,254]
	for i=0,nfeh-1 do for j=0,n_elements(ageplot)-1 do begin
		fehref = fehplot(i)
		ageref = ageplot(j)
		;read data
		fitsfile = dir+"/chisq_age"+sigfig(ageref,3)+"_feh"+sigfig(fehref,3)+".fits"	
		if file_test(fitsfile) eq 0 then stop,'no fits file found'
		chisqarr = mrdfits(fitsfile,0,header)
		minchisq = min(chisqarr)
		contour,chisqarr,grid_feh,grid_age,levels=[minchisq+0.02,minchisq+1],c_linestyle=0,/overplot,color=fehcolor(i)
		oplot,[fehref],alog10([ageref])+9.,psym=7,color=fehcolor(i)
		
	endfor
	device,/close
	stop
end
