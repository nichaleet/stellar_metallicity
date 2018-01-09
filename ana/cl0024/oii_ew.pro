function oii_ew,contdiv,lambda,contdivivar,z,widtherr=widtherr
	fixz=0
	lambdain = lambda
	while fixz eq 0 do begin
		lamrange=[3720.,3735.]
		midlam = mean(lamrange)
		lambda  =lambdain/(z+1.)
		dlambda = lambda-shift(lambda,1)
		inrange = where(lambda gt lamrange[0] and lambda le lamrange[1] and finite(contdiv) and finite(contdivivar),cinrange)
	
		width = total(dlambda(inrange)*(1.-contdiv(inrange)))
		widtherr = sqrt(total((dlambda(inrange))^2/contdiv(inrange)))
		;plotting
		rangeplot=[3700,3750]
		plot,lambda,contdiv,xrange=rangeplot,/nodata
		if width gt 0 then factor = 0 else factor = 1
		polyfill,[midlam-width/2.,midlam+width/2.,midlam+width/2.,midlam-width/2.],[0,0,1,1]+factor,color=fsc_color('yellow')
		oplot,lambda,contdiv
		oplot,[3726,3726],!y.crange,color=fsc_Color('red'),linestyle=1
		oplot,[3729,3729],!y.crange,color=fsc_Color('red'),linestyle=1
		oplot,[lamrange[0],lamrange[0]],!y.crange,color=fsc_color('green'),linestyle=2
		oplot,[lamrange[1],lamrange[1]],!y.crange,color=fsc_color('green'),linestyle=2
		xyouts, 0.2,0.3,'EW= '+sigfig(width,3)+'+/-'+sigfig(widtherr,3),/normal
		strfixz=''
		read,strfixz,prompt='fix redshift? (y/n)'
		if strfixz eq 'y' then fixz=0 else fixz=1
		if fixz eq 0 then begin
			read,newoiiwl,prompt='type in the location of oii that you see: '
			z = (newoiiwl*(z+1)/3728.)-1.
		endif
	endwhile
			
	return,width
end
