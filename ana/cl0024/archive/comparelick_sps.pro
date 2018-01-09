pro comparelick_sps
outdir = '/scr2/nichal/workspace2/ana/cl0024/comparelick_sps/'
maskname='all_cl0024'
;read in lick measurement
lickfile = '/scr2/nichal/workspace2/lickindices/output/results/'+maskname+'_lickmeasurement.fits'
lickstr = mrdfits(lickfile,1,/silent)

;read in sps_measurement
spsfile = '/scr2/nichal/workspace2/sps_fit/data/'+maskname+'/sps_fit.fits.gz'
spsstr = mrdfits(spsfile,1,/silent)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(lickstr) ne n_Elements(spsfile) then print,'oops the two files do not match'
nobj = n_elements(lickstr)
set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
!p.charsize=1
position1 = [0.1,0.74,0.9,0.94]
position2 = [0.1,0.42,0.9,0.64]
position3 = [0.1,0.1,0.9,0.32]
vsym,4,/fill,rot=45

repeatflag = 1
if repeatflag ne 0 then begin
for i=0,nobj-1 do begin
   ;check objname 
   if lickstr(i).good eq 1 and lickstr(i).goodfspsfit eq 1 then begin
      if lickstr(i).objname ne spsstr(i).objname then print,'oops the two files do not match'
      mask = lickstr(i).mask
      slit = lickstr(i).slit
      objname = lickstr(i).objname
      psname=outdir+strtrim(mask,2)+'_'+strtrim(objname,2)+'_'+strtrim(string(slit),2)+'.eps'
      device, filename = psname,xsize = 10,ysize = 15, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;plot
      agethomas = lickstr(i).prob_age[*,0]
      probage  = lickstr(i).prob_age[*,1]
      zthomas = lickstr(i).prob_zh[*,0]
      probz  = lickstr(i).prob_zh[*,1]
      alphathomas = lickstr(i).prob_alpha[*,0]
      probalpha  = lickstr(i).prob_alpha[*,1]
      bestage = lickstr(i).best_age
      bestz = lickstr(i).best_zh
      bestalpha = lickstr(i).best_alpha
      

      spsage = spsstr(i).age
      spsageerr = spsstr(i).ageerr
      spsz = spsstr(i).feh
      spszerr = spsstr(i).feherr
      
      ;plot age
      plot,10^(agethomas-9.),probage,xstyle=4,ystyle=4,position=position1,/nodata
      agesps = findgen(200.)*(max(!x.crange)-min(!x.crange))/200.+min(!x.crange)
      probagesps = gaussian(agesps,[0.8*max(!y.crange),spsage,spsageerr])
      cgcolorfill,agesps,probagesps,color='mediumgray'
      plot,10^(agethomas[*,0,0]-9.),probage,/noerase,xtitle='Age(Gyr)',position=position1,title=strtrim(mask,2)+' '+strtrim(objname,2)+' '+strtrim(string(slit),2)+' '+strtrim(string(spsstr(i).snfit,format='(F4.1)'),2)
      oplot,[spsage],[0.8*max(!y.crange)],psym=8,color=fsc_color('green')
      for ip=1,3 do oplot,10^([bestage[ip],bestage[ip]]-9.),!y.crange,color=fsc_color('red'),linestyle=2
      oplot,10^([bestage[0],bestage[0]]-9.),!y.crange,color=fsc_color('blue'),linestyle=1
      ;plot z
      plot,zthomas,probz,xtitle='[z/H]',position=position2,/nodata,xstyle=4,ystyle=4,/noerase
      zsps = findgen(200.)*(max(!x.crange)-min(!x.crange))/200.+min(!x.crange)
      probzsps = gaussian(zsps,[0.8*max(!y.crange),spsz,spszerr])
      cgcolorfill,zsps,probzsps,color='mediumgray'
      plot,zthomas[*,0,0],probz,/noerase,xtitle='[Z/H]',position=position2
      oplot,[spsz],[0.8*max(!y.crange)],psym=8,color=fsc_color('green')
      for ip=1,3 do oplot,[bestz[ip],bestz[ip]],!y.crange,color=fsc_color('red'),linestyle=2
      oplot,[bestz[0],bestz[0]],!y.crange,color=fsc_color('blue'),linestyle=1
      ;plotalpha
      plot,alphathomas,probalpha,xtitle='[alpha/Fe]',position=position3,/noerase
      for ip=1,3 do oplot,[bestalpha[ip],bestalpha[ip]],!y.crange,color=fsc_color('red'),linestyle=2
      oplot,[bestalpha[0],bestalpha[0]],!y.crange,color=fsc_color('blue'),linestyle=1
      device,/close
   endif
endfor
endif 
;Find the differences in the measurements of age and metallicity
goodobjs = where(lickstr.good eq 1 and lickstr.goodfspsfit eq 1)

lickage = 10^(reform(lickstr(goodobjs).best_age[0,*])-9.)
lickz = reform(lickstr(goodobjs).best_zh[0,*])
lickalpha = reform(lickstr(goodobjs).best_alpha[0,*])

lickageerr = reform(10^(lickstr(goodobjs).best_age[3,*]-9.)>0.1-10^(lickstr(goodobjs).best_age[1,*]-9)<15)/2.
lickzerr = reform(lickstr(goodobjs).best_zh[3,*]<0.67-lickstr(goodobjs).best_zh[1,*]>(-2.25))/2.
lickalphaerr = reform(lickstr(goodobjs).best_alpha[3,*]<0.1-lickstr(goodobjs).best_alpha[1,*]>(-0.3))/2.

spsage = spsstr(goodobjs).age
spsageerr = spsstr(goodobjs).ageerr
spsz = spsstr(goodobjs).feh
spszerr = spsstr(goodobjs).feherr

diffage = lickage-spsage
diffageerr = sqrt(lickageerr^2+spsageerr^2)
diffz = lickz-spsz
diffzerr = sqrt(lickzerr^2+spszerr^2)
  
;plot the differences as a function of alpha, mass, signal to noise
!p.multi=[0,1,3]
!p.charsize=1.5

psname=outdir+'diffage.eps'
device, filename = psname,xsize = 10,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,lickalpha,diffage,lickalphaerr,diffageerr,psym=1,xtitle='[alpha/Fe]',ytitle='lick Age - SPS Age (Gyr)',xrange=[-0.3,0.5]
vsym,4,/fill,rot=45
oplot,lickalpha,diffage,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).snfit,diffage,lickalphaerr,diffageerr,psym=1,xtitle='S/N',ytitle='lick Age - SPS Age (Gyr)'
oplot,spsstr(goodobjs).snfit,diffage,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).logmstar,diffage,lickalphaerr,diffageerr,psym=1,xtitle='log M*',ytitle='lick Age - SPS Age (Gyr)',xrange=[9,11.5]
oplot,spsstr(goodobjs).logmstar,diffage,psym=8,color=fsc_color('blue')
device,/close

psname=outdir+'diffz.eps'
device, filename = psname,xsize = 10,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,lickalpha,diffz,lickalphaerr,diffzerr,psym=1,xtitle='[alpha/Fe]',ytitle='lick [Z/H] - SPS [Z/H]',xrange=[-0.3,0.5]
oplot,lickalpha,diffz,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).snfit,diffz,lickalphaerr,diffzerr,psym=1,xtitle='S/N',ytitle='lick [Z/H] - SPS [Z/H]'
oplot,spsstr(goodobjs).snfit,diffz,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).logmstar,diffz,lickalphaerr,diffzerr,psym=1,xtitle='log M*',ytitle='lick [Z/H] - SPS [Z/H]',xrange=[9,11.5]
oplot,spsstr(goodobjs).logmstar,diffz,psym=8,color=fsc_color('blue')
device,/close
stop
end
