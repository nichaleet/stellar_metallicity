pro csp_ana
  
	mockarr= mrdfits('sps_fit.fits.gz',1)
	restore,'/scr2/nichal/workspace2/mock_data/mock_es_try4/variables.sav'
	mockage =  params_save.age

	;choose age to consider
	age = 5.
	wage = where(mockage eq age,cwage)
	if cwage ne 0 then mock = mockarr(wage) else stop
	z = 0.55
	feh = params_save(wage).z
	reallambda = mock.lambda
	restlambda = mock.lambda/(1.+z)
    won = where(mock.fitmask eq 1 and finite(mock.contdiv) and finite(mock.contdivivar) and mock.contdivivar gt 0 and reallambda/(1.+z) gt 3500. and reallambda/(1.+z) lt 7400., con)
	window,0,xsize=1500,ysize=300
 	plot,restlambda,mock.contdiv,xrange=[3700,4500]
   ;get ssp for various age
	agearr = float(indgen(8)+1)
	colorarr=['purple','cyan','navy','darkgreen','green','yellow','orange','red']
	for i=0,n_Elements(agearr)-1 do begin
		
 	    spsstruct = sps_interp(feh, agearr(i))
        lambda = spsstruct.lambda
        spsspec = spsstruct.spec

    	w = where(lambda gt 0 and lambda lt 100000, c)
	    if c lt 25 then message, 'Not enough pixels.'
	    lambda = lambda[w]
   		spsspec = spsspec[w]
    	clight = 299792.458

    	spsspec = spsspec*clight/lambda^2    ;change fnu(Lsun/Hz) to flambda
    	spsspec = spsspec/median(spsspec)    ;normalize to around 1
	    ;smooth to data wavelengths
		datalam = mock.lambda
		redshift = 0.55
		nlambda = n_Elements(datalam)
		if min(mock.dlam) gt 0.2 and max(mock.dlam) lt 10.0 then begin
        	dlam = mock.dlam
   		endif else dlam = replicate(3.9/2.35, nlambda)
		vdisp = 250./2.35

    	spsspec = smooth_gauss_wrapper(lambda*(redshift+1.), spsspec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))
    	lambda  = datalam ;datalam is mock.lambda

		;continuum normalize them
		bkpt = slatec_splinefit(restlambda[won], spsspec[won], coeff, invvar=mock.contdivivar[won]*(spsspec[won])^2, bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then stop
   
        cont = slatec_bvalu(restlambda, bkpt, coeff)
        mocksps = spsspec/cont
		oplot,restlambda,mocksps,color=fsc_color(colorarr(i))
	endfor	
	;oplot,restlambda,mock.contdiv
		
stop
end
