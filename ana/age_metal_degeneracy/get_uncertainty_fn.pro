pro get_uncertainty_fn
;to get theoretical uncertainties as a function of Fe/H and SN
firstit = 1
;loop over signal to noise
;for iu=8,25 do begin
for iu=20,58,2 do begin
	dir='sdssoutputsn'+strtrim(string(fix(iu)),2)
	;dir='outputsn'+strtrim(string(fix(iu)),2)
	print, 'doing '+dir
	files = file_search(dir+'/chisq*.fits')
	;get values from each file
	for i=0,n_Elements(files)-1 do begin
		;get age and metallicities from file name
		words = strsplit(files(i),'_',/extract)
		age = float(strmid(words[1],3))/10.
		feh = float(strmid(words[2],3,strpos(words[2],'.fits')))/100.
		if feh eq 0.03 then feh = 0.04
		;read the file
		chisqarr = mrdfits(files(i),0,header)
		;you should check below if it matches those in degeneracy.pro
		grid_Feh = rebin(indgen(100)/125.-0.6,100,126)
		grid_age = transpose(rebin(alog10(findgen(126)/10.+0.5)+9.,126,100)) ;range from [0.5,13] Gyr		
		minchisq = min(chisqarr)
		goodchisq = where(chisqarr lt minchisq+1.,ngoodchisq)
		if ngoodchisq lt 2 then stop
		agerange = minmax(grid_age(goodchisq))
		fehrange = minmax(grid_feh(goodchisq))
		str = {age:age,feh:feh,agerange:agerange,fehrange:fehrange}
		if i eq 0 then fixsnstr = str else fixsnstr=[fixsnstr,str]	
	endfor
	feharr = fixsnstr(uniq(fixsnstr.feh,sort(fixsnstr.feh))).feh
	nfeh = n_Elements(feharr)
	sigmaage = fltarr(nfeh) 
	sigmafeh = fltarr(nfeh)
	;run loop over each fixed metallicity
	for i=0,nfeh-1 do begin
		;choose all str with the metallicities, multiple ages)
		selfeh = where(fixsnstr.feh eq feharr(i),cselfeh)
		delage = 0.5*(fixsnstr(selfeh).agerange[1]-fixsnstr(selfeh).agerange[0])
		delfeh = 0.5*(fixsnstr(selfeh).fehrange[1]-fixsnstr(selfeh).fehrange[0])
		sigmaage(i) = median(delage) ;in dex
		sigmafeh(i) = median(delfeh) ;in dex
	endfor
	!p.multi=[0,1,2]
	plot,feharr,sigmaage,xtitle='[Fe/H]',ytitle='sigma log(age)',title=dir
	plot,feharr,sigmafeh,xtitle='[Fe/H]',ytitle='sigma FEH',yrange=[0,0.2]
	wait,1
	outstr = {sn:iu,feh:feharr,sigmaage:sigmaage,sigmafeh:sigmafeh}
	if firstit eq 1 then finstr = outstr else finstr=[finstr,outstr]
	firstit = 0
endfor

;mwrfits,finstr,'sdss_uncertainty_fn_feh_sn.fits',/silent,/create
end
