pro make_mosfire_xtcalc_template

	strb = mrdfits('/scr2/nichal/workspace2/mock_data/mock_es_try4/spec1d.csp_try4.006.mockspec.fits',1)
	strr = mrdfits('/scr2/nichal/workspace2/mock_data/mock_es_try4/spec1d.csp_try4.006.mockspec.fits',2)
	lambdared = strr.lambda
	lambdablue = strb.lambda
	lambda=[strb.lambda,strr.lambda]
	spec=[strb.spec,strr.spec]
	;fix deimos throughput
	readcol,'thr_g900_70_gg455.asc',wltput,tput
    tputred  = interpol(tput,wltput,lambdared)
    tputblue = interpol(tput,wltput,lambdablue)
    tputred = tputred/median(tputred)
    tputblue = tputblue/median(tputblue)
	spec = [strb.spec/tputblue,strr.spec/tputred]
	lambda = [lambdablue,lambdared]/(strb.redshift+1.)


	writecol,'/scr2/nichal/XTcalc/csp_4Gyr_template.txt',lambda,spec,fmt='(F6.1,2x,F7.5)'
	stop
end
