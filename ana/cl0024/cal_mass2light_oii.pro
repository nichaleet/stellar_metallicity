pro cal_mass2light_oii
;calculate mass to light ratio from SPS spectra at OII continuum wl
;read in SPS spectra with Z=-0.103000 
	arr = sps_read_spec('/scr2/nichal/workspace2/fsps-3.0/SSP/SSP_Padova_MILES_Kroupa_zmet19.out.spec')
	selage = [0.1,0.5,1,2] ;Gyr

	solarlum = 3.839d+33 ;erg per second
	clight = 3.d18 ;angstrom per second
	for i=0,n_elements(selage)-1 do begin
		mindiff = min(abs(arr.agegyr-selage[i]),locmin)
		arrnow = arr(locmin)
		spec = arrnow.spec
		lambda = arrnow.lambda
		goodlamb = where(lambda gt 3600. and lambda lt 5500.)
		spec = spec(goodlamb)
		lambda = lambda(goodlamb)

		;find continuum
		;masking
		readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
		mask = bytarr(n_elements(lambda))+1
		for j=0,n_elements(linestart)-1 do begin
			w = where(lambda ge linestart[j] and lambda le lineend[j], c)
			if c gt 0 then mask[w] = 0
		endfor
		won = where(mask eq 1, complement=woff, con)
		bkpt = slatec_splinefit(lambda[won], spec[won], coeff,bkspace=200, /silent)
		if bkpt[0] eq -1 then stop,'oh no'
		cont = slatec_bvalu(lambda, bkpt, coeff)
		;plot
		plot,lambda,spec
		oplot,lambda,cont,color=fsc_Color('green')
		print, 'Age = '+sigfig(arrnow.agegyr,3)+' Gyr'
		;oii continuum flux
		oiireg = where(lambda gt 3720 and lambda lt 3730)
		oiicont = mean(cont(oiireg))
		Loii = oiicont*solarlum  ;erg/s/hz
		;change to erg/s/angstrom
		Loii =  Loii*clight/(3727.)^2 ;erg/s/angstrom
		mass = 10^arrnow.LOGMASS ;msun
		
		mass2light = mass/Loii
		print, 'M/L = '+sigfig(mass2light,4,/sci)+' Msun s A/erg'
		
		;Calculate sSFR based on Kewley04
		alpha = 6.58d-42
		EW = 5 ;EW of Oii = 5 A
		sSFR = alpha*EW/mass2light
		print, 'sSFR ='+sigfig(sSFR,4,/sci)+' per year'
		print, ''
		wait,2
	endfor		
end
