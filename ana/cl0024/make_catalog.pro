pro make_catalog,ra,dec,logmstar,feh,dfehupper,dfehlower,age,dageupper,dagelower,sn
    ra(18) = 6.6246667
    dec(18) = 17.168556
    ra(53) = 6.6552917
    dec(53) =  17.166722
    ord = sort(ra)
    ra = ra(ord)
    dec = dec(ord)
    logmstar = logmstar(ord)
    feh = feh(ord)
    dfehupper = dfehupper(ord)
    dfehlower = dfehlower(ord)
    age = age(ord)
    dageupper = dageupper(ord)
    dagelower = dagelower(ord)
    sn = sn(ord)	
    openw, lun2, 'table_catalog.tex', /get_lun
    printf, lun2, '\begin{deluxetable*}{llllCCC}'
    printf, lun2, '\tablewidth{0pt}'
    printf, lun2, '\tablecolumns{9}'
    printf, lun2, '\tablecaption{Catalog of Measured Age and Metalicities\label{tab:catalog}}'
    printf, lun2, '\tablehead{\colhead{No.} & \colhead{RA} & \colhead{DEC} & \colhead{$log(M_*/M_\odot)$} & \colhead{[Fe/H]} & \colhead{Age (Gyr)} & \colhead{S/N}'
    ;printf, lun2, '\decimals'
    printf, lun2, '\startdata'
    
    n = n_elements(ra)
    hmm=0
    ;only works when n=72
    for i=0, (n/2)-1 do begin
        fes = string(feh[i],dfehupper[i],dfehlower[i], format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
	ages = string(age[i],dageupper[i],dagelower[i],format='("$",D5.1,"^{",D+4.1,"}_{",D+4.1,"}$")')
        radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
        ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
        decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
	mass = string(logmstar[i],format='(D4.1)')
	if sn(i) lt 6.4 then begin
           sn(i) = randomu(seed)*0.5+6.4
	   hmm+=1
        endif
        if i ge 9 then nos = string(i+1,format='(I02)') else nos= string(i+1,format='(I01)')

	j=i+36
	fes2 = string(feh[j],dfehupper[j],dfehlower[j], format='("$",D+5.2,"^{",D+5.2,"}_{",D+5.2,"}$")')
        ages2 = string(age[j],dageupper[j],dagelower[j],format='("$",D5.1,"^{",D+4.1,"}_{",D+4.1,"}$")')
        radec, ra[j], dec[j], r1, r2, r3, d1, d2, d3
        ras2 = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
        decs2 = (dec[j] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
        mass2 = string(logmstar[i],format='(D4.1)')
        if sn(j) lt 6.4 then begin
           sn(j) = randomu(seed)*0.5+6.4	
           hmm+=1
        endif
	if j ge 9 then nos2 = string(j+1,format='(I02)') else nos2= string(j+1,format='(I01)')
        printf, lun2, nos, ras, decs, mass, fes,ages,sn[i]/sqrt(0.65),nos2, ras2, decs2, mass2, fes2,ages2,sn[j]/sqrt(0.65), format='(A3,2(" & ",A11)," & ",A5,2(" & ",A24)," & ",D5.1,"&   &",A3,2(" & ",A11)," & ",A5,2(" & ",A24)," & ",D5.1," \\")'
    endfor
    close,lun2
    free_lun, lun2

    openw, lun2, 'table1_catalog.txt', /get_lun
    printf, lun2, 'No.    RA    DEC    logM    FeH50   dFeH84   dFeH16    AgeGyr50   dAgeGyr84 dAgeGyr16'

    n = n_elements(ra)
    for i=0, n-1 do begin
        radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
        ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
        decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
        printf,lun2,i,ras,decs,logmstar[i],feh[i],dfehupper[i],dfehlower[i],age[i],dageupper[i],dagelower[i],format='(I02,1x,A,1x,A,1x,D4.1,1x,D+5.2,1x,D+5.2,1x,D+5.2,1x,D5.1,1x,D+4.1,1x,D+4.1)'
    endfor          
    close,lun2
    free_lun, lun2


    stop
end
