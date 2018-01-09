pro getobslick_sum,maskname,name,plot=plot,briefplot=briefplot,objno=objno,indexno=indexno
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Purpose: get the lick indices of all galaxies of specified mask names
;Only use the indices that can accurately reporduce those of globular cluster in Table 2 of Thomas11b (same indices as in Onodera15)
;The indices definitions are according to Thomass,Maraston,and Johansson 2010 which are from Trager1998 and Worthey&Ottaviani for HdA,HdF,HgA,HgF
; Set 'plot' for plotting everything
; Set 'brief plot' for plotting histogram of equivalent width only
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
   for io=0,nobs-1 do begin ;io = index of object
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
      ;mask the pixel whose flux is outliers or with large uncertainty
      mask = bytarr(npix)
      outliers = where(abs(contdiv-spsspec)/spsspec gt 3. or abs(contdiv-spsspec)/contdiverr gt 3. or contdivivar eq 0.,cout)
      if cout gt 0 then mask(outliers) = 1
      ;mask the middle parts of Hydrogen line
      ;mask for 2 A
      hdAreg = where(lambda ge 4101. and lambda le 4103,creg)
      if creg gt 0 then mask(hdareg) = 1
      hgAreg = where(lambda ge 4340. and lambda le 4342,creg)
      if creg gt 0 then mask(hdareg) = 1

      ;begin loop over indices
      if inds[io].good eq 1 then begin
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
                  fluxpseudo = fluxblue+(lambpass-mcontb[ii])*(fluxred-fluxblue)/(mcontr[ii]-mcontb[ii])
                  slope = (fluxred-fluxblue)/(mcontr[ii]-mcontb[ii])
                  const = fluxblue-(slope*mcontb[ii])
                                ;calculate EW
                  if indunit[ii] eq 0 then EW = total(((1.-fluxpass/fluxpseudo))*dlambpass) ;angstrom
                  if indunit[ii] eq 1 then EW = -2.5*alog10(total(((fluxpass/fluxpseudo))*dlambpass)/(max(lambpass)-min(lambpass))) ;magnitude
                  EWarr[j] = EW
                  if keyword_set(plot) then begin
                     plot,lambda(alllamb),contdiv(alllamb),xtitle='lambda'
                     oplot,lambcontb,fluxcontb,psym=10,color=fsc_color('cyan')
                     oplot,lambpass,fluxpass,psym=10,color=fsc_color('green')
                     oplot,lambcontr,fluxcontr,psym=10,color=fsc_color('red')
                     oplot,lambda,lambda*slope+const,linestyle=2

                     if indunit[ii] eq 0 then oplot,[mpass[ii]-EW/2.,mpass[ii]+EW/2],[0.1,0.1],thick=4,color=fsc_color('orange')
                     xyouts,0.2,0.9,maskname+' objno.'+strtrim(string(io),2)+'/'+strtrim(string(nobs),2),alignment=0,/normal,color=fsc_color('red')
                     xyouts,0.2,0.85,'index '+strtrim(string(ii+1),2)+'/16 '+indname[ii]+' iteration:'+strtrim(string(j),2),alignment=0,/normal,color=fsc_color('red')
                     xyouts,0.2,0.8,'EW ='+strtrim(string(EW,format="(F5.2)"),2)+unit,alignment=0,/normal,color=fsc_color('red')

                     wait,0.05
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

                  wait,0.25
               endif
               match = strmatch(tagout,indname[ii],/fold_case)
               if total(match) ne 1 then stop
               postag = where(match eq 1)
               inds[io].(postag) = [EW,EWerr]
            endif
         endfor
      endif
   endfor
   outfile = outdir+maskname[nm]+'_lickindices_sum.fits'
   mwrfits,inds,outfile,/create,/silent
endfor

end
