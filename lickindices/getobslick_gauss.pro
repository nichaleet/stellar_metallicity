pro gaussian7,x,A,F,pder
  if A[0] le 0 and A[5] le 0 then begin
     expterm = EXP(-0.5*((x-A[1])/A[2])^2)
     expterm2= EXP(-0.5*((x-A[6])/A[2])^2)
     F = A[0]*expterm+A[3]+A[4]*x+A[5]*expterm2
  endif else F = fltarr(n_elements(x))-99.
  if n_params() ge 4 then $
     pder = [[expterm],[A[0]*expterm*(x-A[1])/(A[2])^2],[A[0]*expterm*(x-A[1])^2/(A[2])^3],[replicate(1.0, N_ELEMENTS(X))],[x],[expterm2],[A[5]*expterm2*(x-A[6])/(A[2])^2]]
end

pro gaussian9,x,A,F,pder
  if A[0] le 0 and A[5] le 0 and A[7] le 0 then begin
     expterm = EXP(-0.5*((x-A[1])/A[2])^2)
     expterm2= EXP(-0.5*((x-A[6])/A[2])^2)
     expterm3= EXP(-0.5*((x-A[8])/A[2])^2)
     F = A[0]*expterm+A[3]+A[4]*x+A[5]*expterm2+A[7]*expterm3
  endif else F=fltarr(n_elements(x))-99.
  if n_params() ge 4 then $
     pder = [[expterm],[A[0]*expterm*(x-A[1])/(A[2])^2],[A[0]*expterm*(x-A[1])^2/(A[2])^3],[replicate(1.0, N_ELEMENTS(X))],[x],[expterm2],[A[5]*expterm2*(x-A[6])/(A[2])^2],[expterm3],[A[7]*expterm3*(x-A[8])/(A[2])^2]]
  
end

pro gauss_emission,x,A,F,pder
  if A[0] gt 0 and A[4] lt 0 and A[4] gt -1. and A[2] lt 700./300000.*A[1] and A[5] lt 400./300000.*A[1] and a[3] lt 0.1 and a[3] gt 0. then begin
     emission = A[0]*EXP(-0.5*((x-A[1])/A[2])^2)+A[3]
     absorption= A[4]*EXP(-0.5*((x-A[1])/A[5])^2)+1.-A[3]
     F = emission+absorption
  endif else F=fltarr(n_elements(x))-99.
  
  if n_params() ge 4 then $
     pder = [[EXP(-0.5*((x-A[1])/A[2])^2)],[(emission-A[3])*(x-A[1])/(A[2])^2],[(emission-A[3])*(x-A[1])^2/(A[2])^3],[replicate(1.0, N_ELEMENTS(X))],[EXP(-0.5*((x-A[1])/A[5])^2)],[A[4]*EXP(-0.5*((x-A[1])/A[5])^2)*(x-A[1])^2/(A[5])^3]]

end

pro getobslick_gauss,maskname,name,plot=plot,briefplot=briefplot,objno=objno,indexno=indexno
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Purpose: get the lick indices of all galaxies of specified mask names
;Only use the indices that can accurately reporduce those of globular cluster in Table 2 of Thomas11b (same indices as in Onodera15)
;The indices definitions are according to Thomass,Maraston,and Johansson 2010 which are from Trager1998 and Worthey&Ottaviani for HdA,HdF,HgA,HgF
; Set 'plot' for plotting everything
; Set 'brief plot' for plotting histogram of equivalent width only
; EXP: getobslick_gauss,'all_cl0024','all_cl0024',/briefplot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;set number of iterations per each index
niter = 500.

;input directory
directory = '/scr2/nichal/workspace2/sps_fit/data/'
nmasks = n_elements(maskname)

;output directory
outdir = '/scr2/nichal/workspace2/lickindices/output/ew/'

;read index definition.
deffile= '/scr2/nichal/workspace2/lickindices/lick_thomas11b.txt'
readcol,deffile,indno,indpassi,indpassf,indcontbi,indcontbf,indcontri,indcontrf,indunit,indname,indtype,comment='#',format='L,D,D,D,D,D,D,B,A,A'
mcontb = 0.5*(indcontbi+indcontbf)
mcontr = 0.5*(indcontri+indcontrf)
mpass  = 0.5*(indpassi+indpassf)
nlines = n_elements(indno)
indgauss = intarr(nlines)
indgauss_line1 = fltarr(nlines)-99.
indgauss_line2 = fltarr(nlines)-99.
indgauss_line3 = fltarr(nlines)-99.
for i=0,nlines-1 do begin
   if strmid(indtype[i],0,1) eq 'S' then indgauss[i]=0
   if strmid(indtype[i],0,1) eq 'G' then begin
      indgauss[i]=1
      word = strsplit(indtype[i],'-',/extract)
      indgauss_line1[i] = float(word[1])
   endif
   if strmid(indtype[i],0,1) eq 'D' then begin
      indgauss[i]=2
      word = strsplit(indtype[i],'-',/extract)
      indgauss_line1[i] = float(word[1])
      indgauss_line2[i] = float(word[2])
   endif
   if strmid(indtype[i],0,1) eq 'T' then begin
      indgauss[i]=3
      word = strsplit(indtype[i],'-',/extract)
      indgauss_line1[i] = float(word[1])
      indgauss_line2[i] = float(word[2])
      indgauss_line3[i] = float(word[3])
   endif
endfor

!p.font = 0
for nm=0,nmasks-1 do begin
   scifits = directory+maskname[nm]+'/sps_fit.fits.gz'
   sciall = mrdfits(scifits,1,/silent)   
   nobs = n_elements(sciall)
   goodobj = where(sciall.good eq 1 and sciall.zfit ge 0.38 and sciall.zfit le 4.0, cgals)
   ;output structure for indices
   indstr = {fitsfile:scifits,objname:'',mask:'',slit:0L,ra:0.d,dec:0.d,jdobs:-99.,zspec:-99.,good:0B,goodfspsfit:0B,CN1:[-99.,-99.],CN2:[-99.,-99.],Ca4227:[-99.,-99.],G4300:[-99.,-99.],C24668:[-99.,-99.],Mg1:[-99.,-99.],Mg2:[-99.,-99.],Mgb:[-99.,-99.],Fe4383:[-99.,-99.],Fe4531:[-99.,-99.],Fe5270:[-99.,-99.],Fe5335:[-99.,-99.],Fe5406:[-99.,-99.],HdA:[-99.,-99.],HgA:[-99.,-99.],HgF:[-99.,-99.]}
   tagout = tag_names(indstr)
   inds = replicate(indstr,nobs)
   for io=0,nobs-1 do begin     ;io = index of object
      inds[io].objname = sciall[io].objname
      inds[io].mask = sciall[io].mask
      inds[io].slit = sciall[io].slit
      inds[io].ra = sciall[io].ra
      inds[io].dec = sciall[io].dec
      inds[io].jdobs = sciall[io].jdobs
      inds[io].zspec = sciall[io].zspec
      inds[io].good = sciall[io].good
      inds[io].goodfspsfit = sciall[io].goodfit
      ;get spectra
      zspec = sciall[io].zspec
      if sciall[io].zfit ne 0 then zspec = sciall[io].zfit
      lambda  = sciall[io].lambda/(1.+zspec)
      contdiv = sciall[io].contdiv
      contdivivar = sciall[io].contdivivar
      contdiverr = 1./sqrt(contdivivar)
      spsspec = sciall[io].spsspec
     ; contdiv = sciall[io].telldiv
     ; contdivivar = sciall[io].telldivivar
     ; contdiverr = 1./sqrt(contdivivar)
     ; spsspec = sciall[io].spsspec*sciall[io].continuum
      npix = n_elements(lambda)
 
      ;begin loop over indices
      if inds[io].good eq 1 then begin

         ;mask the pixel whose flux is outliers or with large uncertainty
         mask = bytarr(npix)
         outliers = where(abs(contdiv-spsspec)/spsspec gt 3. or abs(contdiv-spsspec)/contdiverr gt 3. or contdivivar eq 0.,cout)
         if cout gt 0 then mask(outliers) = 1

         ;mask the middle parts of Hydrogen line
         ;look for the emission line. 
         ;If there is, fit and subtract the emission line.       
         !p.multi=[0,1,2]
         ;Hda line
         hreg = where(lambda ge 4098. and lambda le 4109,creg)
         widehreg = [-1.*reverse(indgen(5)+1)+min(hreg),hreg,indgen(5)+1+max(hreg)]
         a = [0.5,4103,2.,0.05,-0.3,3]
         if max(contdiv(hreg)) gt 2. then a[0]=max(contdiv(hreg))-1
         fita = [1,1,1,1,1,1]
         hfit = curvefit(lambda(hreg),contdiv(hreg),contdivivar(hreg),a,fita=fita,function_name='gauss_emission',itmax=1e5,status=status)
         plot,lambda(widehreg),contdiv(widehreg),psym=10
         if a[0] gt 0.05 and abs(a[1]-4103.) lt 2. and status ne 1 then begin
            emission = A[0]*EXP(-0.5*((lambda(hreg)-A[1])/A[2])^2)+A[3]
            absorption = A[4]*EXP(-0.5*((lambda(hreg)-A[1])/A[5])^2)+1.-A[3]
            oplot,lambda(hreg),emission,color=fsc_color('green')            
            oplot,lambda(hreg),absorption,color=fsc_color('green')            
            contdiv(hreg) = absorption/(1.-A[3])
            oplot,lambda(hreg),contdiv(hreg),color=fsc_color('orange'),psym=10
         endif
         if status eq 1 then begin
            hfit = gaussfit(lambda(hreg),contdiv(hreg),a,nterms=4)
            if a[0] gt 0.05 and abs(a[1]-4103.) lt 1.5 then begin
               oplot,lambda(hreg),hfit,color=fsc_color('green')
               contdiv(hreg) = contdiv(hreg)-hfit+a[3]
               oplot,lambda(hreg),contdiv(hreg),color=fsc_color('orange'),psym=10
            endif 
         endif

         ;hga line
         hreg = where(lambda ge 4338. and lambda le 4346,creg)
         widehreg = [-1.*reverse(indgen(5)+1)+min(hreg),hreg,indgen(5)+1+max(hreg)]
         a = [0.5,4341,2.,0.05,-0.3,3]
         if max(contdiv(hreg)) gt 2. then a[0]=max(contdiv(hreg))-1
         fita = [1,1,1,1,1,1]
         hfit = curvefit(lambda(hreg),contdiv(hreg),contdivivar(hreg),a,fita=fita,function_name='gauss_emission',itmax=1e5,status=status)
         plot,lambda(widehreg),contdiv(widehreg),psym=10 
         if a[0] gt 0.05 and abs(a[1]-4342.) lt 2 and status ne 1 then begin
            emission = A[0]*EXP(-0.5*((lambda(hreg)-A[1])/A[2])^2)+A[3]
            absorption = A[4]*EXP(-0.5*((lambda(hreg)-A[1])/A[5])^2)+1.-A[3]
            oplot,lambda(hreg),emission,color=fsc_color('green')
            oplot,lambda(hreg),absorption,color=fsc_color('green')
            contdiv(hreg) = absorption/(1.-A[3])
            oplot,lambda(hreg),contdiv(hreg),color=fsc_color('orange'),psym=10
         endif
         if status eq 1 then begin
            hfit = gaussfit(lambda(hreg),contdiv(hreg),a,nterms=4)
            if a[0] gt 0.05 and abs(a[1]-4342.) lt 1.5 then begin
               oplot,lambda(hreg),hfit,color=fsc_color('green')
               contdiv(hreg) = contdiv(hreg)-hfit+a[3]
               oplot,lambda(hreg),contdiv(hreg),color=fsc_color('orange'),psym=10
            endif
         endif
         wait,2
         !p.multi=[0,1,1]
         for ii=0,15 do begin   ;ii= index of indices
            incontb = where(lambda ge indcontbi[ii] and lambda le indcontbf[ii] and mask eq 0,cb)
            incontr = where(lambda ge indcontri[ii] and lambda le indcontrf[ii]and mask eq 0,cr)
            inpass =  where(lambda ge indpassi[ii] and lambda le indpassf[ii]and mask eq 0,cp)   
            alllamb = where(lambda ge indcontbi[ii]-5. and lambda le indcontrf[ii]+5.)
            if cr ge 5 and cb ge 5 and cp ge 5 then begin 
               lambcontb = lambda(incontb)
               dlambcontb = -1.*ts_diff([lambcontb,lambda(incontb[cb-1]+1)],1)
               dlambcontb = dlambcontb[0:cb-1]

               lambcontr = lambda(incontr)
               dlambcontr = -1.*ts_diff([lambcontr,lambda(incontr[cr-1]+1)],1)
               dlambcontr = dlambcontr[0:cr-1]

               lambpass = lambda(inpass)
               dlambpass = -1.*ts_diff([lambpass,lambda(inpass[cp-1]+1)],1)
               dlambpass = dlambpass[0:cp-1]
               if indunit[ii] eq 0 then unit=' A' else unit = ' mag'
               EWarr = fltarr(niter)
              ;iterate niter times
               for j=0,niter-1 do begin                  
                  fluxcontb = contdiv(incontb)+randomn(seed,cb)*contdiverr(incontb)
                  fluxcontr = contdiv(incontr)+randomn(seed,cr)*contdiverr(incontr)
                  fluxpass  = contdiv(inpass)+randomn(seed,cp)*contdiverr(inpass)
                                ;find average pseudocontinuum flux
                  fluxblue =  total(fluxcontb*dlambcontb)/(max(lambcontb)-min(lambcontb))
                  fluxred =  total(fluxcontr*dlambcontr)/(max(lambcontr)-min(lambcontr))
                  
                                ;linear equation for pseudocontinuum line
                  slope = (fluxred-fluxblue)/(mcontr[ii]-mcontb[ii])
                  const = fluxblue-(slope*mcontb[ii])
                  fluxpseudo = fluxblue+(lambpass-mcontb[ii])*slope
                                ;fitting gaussian to the indices
                  if indgauss[ii] eq 0 then fluxpass_gauss = fluxpass
                  if indgauss[ii] eq 1 then begin
                     A = [-0.3,indgauss_line1[ii],2.,const,slope,0,indgauss_line1[ii]]
                     fita = [1,0,1,0,0,0,0] ;fix the peak wl, const and slope
                     fluxpass_gauss = curvefit(lambpass,fluxpass,contdiverr(inpass),A,fita=fita,function_name='gaussian7',status=status,itmax=1e5)                 
                     ;if status ne 0 then fluxpass_gauss = fluxpass
                     if status ne 0 then print,'GAUSSIAN FIT FAILED at '+'objno.'+strtrim(string(io),2)+' '+indname[ii]
                  endif
                  if indgauss[ii] eq 2 then begin
                     A = [-0.3,indgauss_line1[ii],2.,const,slope,-0.3,indgauss_line2[ii]]
                     fita = [1,0,1,0,0,1,0] ;fix the peak wl, const and slope
                     fluxpass_gauss = curvefit(lambpass,fluxpass,contdiverr(inpass),A,fita=fita,function_name='gaussian7',status=status,itmax=1e5) 
                     ;if status ne 0 then fluxpass_gauss = fluxpass
                     if status ne 0 then print,'GAUSSIAN FIT FAILED at '+'objno.'+strtrim(string(io),2)+' '+indname[ii]
                  endif
                  if indgauss[ii] eq 3 then begin
                     A = [-0.3,indgauss_line1[ii],2.,const,slope,-0.3,indgauss_line2[ii],-0.3,indgauss_line3[ii]]
                     fita = [1,0,1,0,0,1,0,1,0] ;fix the peak wl, const and slope
                     fluxpass_gauss = curvefit(lambpass,fluxpass,contdiverr(inpass),A,fita=fita,function_name='gaussian9',status=status,itmax=1e5) 
                     ;if status ne 0 then fluxpass_gauss = fluxpass
                     if status ne 0 then print,'GAUSSIAN FIT FAILED at '+'objno.'+strtrim(string(io),2)+' '+indname[ii]
                  endif
                                ;calculate EW 
                  if indunit[ii] eq 0 then EW = total(((1.-fluxpass_gauss/fluxpseudo))*dlambpass) ;angstrom
                  if indunit[ii] eq 1 then EW = -2.5*alog10(total(((fluxpass_gauss/fluxpseudo))*dlambpass)/(max(lambpass)-min(lambpass))) ;magnitude
                  EWarr[j] = EW
                  if keyword_set(plot) and j mod 10 eq 0 then begin
                     plot,lambda(alllamb),contdiv(alllamb),xtitle='lambda',yrange=[0,2]
                     oplot,lambcontb,fluxcontb,psym=10,color=fsc_color('cyan')
                     oplot,lambpass,fluxpass,psym=10,color=fsc_color('green')
                     oplot,lambpass,fluxpass_gauss,color=fsc_color('orange')
                     oplot,lambcontr,fluxcontr,psym=10,color=fsc_color('red')
                     oplot,[indgauss_line1[ii],indgauss_line1[ii]],!y.crange,linestyle=1,color=fsc_color('white')
                     oplot,[indgauss_line2[ii],indgauss_line2[ii]],!y.crange,linestyle=1,color=fsc_color('white')
                     oplot,[indgauss_line3[ii],indgauss_line3[ii]],!y.crange,linestyle=1,color=fsc_color('white')
                     oplot,lambda,lambda*slope+const,linestyle=2

                     if indunit[ii] eq 0 then oplot,[mpass[ii]-EW/2.,mpass[ii]+EW/2],[0.1,0.1],thick=4,color=fsc_color('orange')
                     xyouts,0.2,0.9,maskname+' objno.'+strtrim(string(io),2)+'/'+strtrim(string(nobs),2),alignment=0,/normal,color=fsc_color('red')
                     xyouts,0.2,0.85,'index '+strtrim(string(ii+1),2)+'/16 '+indname[ii]+' iteration:'+strtrim(string(j),2),alignment=0,/normal,color=fsc_color('red')
                     xyouts,0.2,0.8,'EW ='+strtrim(string(EW,format="(F5.2)"),2)+unit,alignment=0,/normal,color=fsc_color('red')
                     xyouts,0.2,0.75,'Gaussian type '+indtype[ii],alignment=0,/normal,color=fsc_color('red')
                     if indgauss[ii] eq 0 then wait,0.05
                     if indgauss[ii] eq 1 then wait,0.1
                     if indgauss[ii] eq 2 then wait,0.2
                     if indgauss[ii] eq 3 then wait,0.2
                     
                  endif
               endfor
               ;find the mean EW and standard deviation
               EW = mean(EWarr)
               EWerr = stdev(EWarr)
               if keyword_set(plot) or keyword_set(briefplot) then begin
                  plothist,EWarr,xtitle='Equivalent width ('+unit+')',bin=Ewerr/2.
                  oplot,[EW,EW],!y.crange,color=fsc_Color('red')
                  oplot,[EW-EWerr,EW-EWerr],!y.crange,color=fsc_Color('red'),linestyle=2
                  oplot,[EW+EWerr,EW+EWerr],!y.crange,color=fsc_Color('red'),linestyle=2
                  xyouts,0.2,0.9,sciall[io].mask+' objno.'+strtrim(string(io),2)+'/'+strtrim(string(nobs),2),alignment=0,/normal,color=fsc_color('red')
                  xyouts,0.2,0.85,'index '+strtrim(string(ii+1),2)+'/16 '+indname[ii],alignment=0,/normal,color=fsc_color('red')
                  xyouts,0.2,0.8,'EW ='+strtrim(string(EW,format="(F5.2)"),2)+'+/-'+strtrim(string(EWerr,format="(F5.2)"),2)+unit,alignment=0,/normal,color=fsc_color('red')
                  ;wait,1
               endif
               match = strmatch(tagout,indname[ii],/fold_case)
               if total(match) ne 1 then stop
               postag = where(match eq 1)
               inds[io].(postag) = [EW,EWerr]
            endif
         endfor
      endif
   endfor
   outfile = outdir+maskname[nm]+'_lickindices_gauss_nomaxflux_emisssub.fits'
   mwrfits,inds,outfile,/create,/silent
endfor

end
