pro comparelick_sps_fixalpha,lickfile,lickfixalphafile
;e.g; comparelick_sps_fixalpha,'all_cl0024_gauss_nomaxflux_emisssub_lickmeasurement.fits','all_cl0024_gauss_nomaxflux_emisssub_lickmeasurement_fixedalpha.fits'
outdir = '/scr2/nichal/workspace2/ana/cl0024/comparelick_sps/fixalpha/'
;read in lick measurement
lickfile = '/scr2/nichal/workspace2/lickindices/output/results/'+lickfile
lickstr = mrdfits(lickfile,1,/silent)

;read in lick measurement with fixed alpha =0
lickfixfile = '/scr2/nichal/workspace2/lickindices/output/results/'+lickfixalphafile
lickfixstr = mrdfits(lickfixfile,1,/silent)

;read in sps_measurement
spsfile = '/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz'
spsstr = mrdfits(spsfile,1,/silent)

;read index definition.
deffile= '/scr2/nichal/workspace2/lickindices/lick_thomas11b.txt'
readcol,deffile,indno,indpassi,indpassf,indcontbi,indcontbf,indcontri,indcontrf,indunit,indname,indtype,comment='#',format='L,D,D,D,D,D,D,B,A,A'

;get the thomas11b prediction grid 
res=0.025
filename = '/scr2/nichal/workspace2/lickindices/thomas_fine'+strtrim(string(res,format='(F5.3)'),2)+'.sav'
filethomas = file_search(filename,count=fcount)
restore,filethomas[0]
agearr = reform(thomasgrid[*,0,0].age)
zarr = reform(thomasgrid[0,*,0].zh)
alphaarr = reform(thomasgrid[0,0,*].alphafe)
tagthomas = tag_names(thomasgrid)
tagmatch  = [where(strmatch(tagthomas,'cn1',/fold_case) eq 1),where(strmatch(tagthomas,'cn2',/fold_case) eq 1),where(strmatch(tagthomas,'ca4227',/fold_case) eq 1),where(strmatch(tagthomas,'g4300',/fold_case) eq 1),where(strmatch(tagthomas,'c24668',/fold_case) eq 1),where(strmatch(tagthomas,'mg1',/fold_case) eq 1),where(strmatch(tagthomas,'mg2',/fold_case) eq 1),where(strmatch(tagthomas,'mgb',/fold_case) eq 1),where(strmatch(tagthomas,'fe4383',/fold_case) eq 1),where(strmatch(tagthomas,'fe4531',/fold_case) eq 1),where(strmatch(tagthomas,'fe5270',/fold_case) eq 1),where(strmatch(tagthomas,'fe5335',/fold_case) eq 1),where(strmatch(tagthomas,'fe5406',/fold_case) eq 1),where(strmatch(tagthomas,'hda',/fold_case) eq 1),where(strmatch(tagthomas,'hga',/fold_case) eq 1),where(strmatch(tagthomas,'hgf',/fold_case) eq 1)]
tagthomasout = tagthomas(tagmatch)
if n_elements(tagmatch) ne 16 then stop, 'Oh no'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if n_elements(lickstr) ne n_Elements(spsstr) then print,'oops the two files do not match'
nobj = n_elements(lickstr)
set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
!p.charsize=1
position1 = [0.05,0.64,0.47,0.94]
position2 = [0.53,0.64,0.95,0.94]
position3 = [0.05,0.26,0.35,0.58]
position4 = [0.4,0.45,0.95,0.58]
position5 = [0.4,0.26,0.95,0.39]
position6 = [0.05,0.05,0.95,0.2]

repeatflag = 1
if repeatflag ne 0 then begin
for i=0,nobj-1 do begin
   ;check objname 
   if lickstr(i).good eq 1 and lickstr(i).goodfspsfit eq 1 then begin
      if lickstr(i).objname ne spsstr(i).objname then print,'oops the two files do not match'
      mask = lickstr(i).mask
      slit = lickstr(i).slit
      objname = lickstr(i).objname
      psname=outdir+strtrim(mask,2)+'_'+strtrim(objname,2)+'_'+strtrim(string(slit),2)+'_fixalpha.eps'
      device, filename = psname,xsize = 30,ysize = 20, $
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
      
      agethomasfix = lickfixstr(i).prob_age[*,0]
      probagefix  = lickfixstr(i).prob_age[*,1]
      zthomasfix = lickfixstr(i).prob_zh[*,0]
      probzfix  = lickfixstr(i).prob_zh[*,1]
      bestagefix = lickfixstr(i).best_age
      bestzfix = lickfixstr(i).best_zh

      spsage = spsstr(i).age
      spsageerr = spsstr(i).ageerr
      spsz = spsstr(i).feh
      spszerr = spsstr(i).feherr

      ;plot age
      plot,10^(agethomas-9.),probage,xstyle=4,ystyle=4,position=position1,/nodata
      agesps = findgen(200.)*(max(!x.crange)-min(!x.crange))/200.+min(!x.crange)
      maxy = max(!y.crange)
      probagesps = gaussian(agesps,[0.8*maxy,spsage,spsageerr])
      cgcolorfill,agesps,probagesps,color=fsc_color('mediumgray')
      plot,10^(agethomas-9.),probage,/noerase,xtitle='Age(Gyr)',position=position1,ytitle='P(Age)/dex'
      vsym,4,/fill,rot=45
      oplot,[spsage],[0.8*maxy],psym=8,color=fsc_color('green')
      oplot,10^([bestage[0],bestage[0]]-9.),[0,maxy],linestyle=1,color=fsc_color('black')
      maxprob = max(probagefix)
      oplot,10^(agethomasfix[*,0]-9.),probagefix/maxprob*0.9*max(!y.crange),color=fsc_color('red')
      for ip=0,3 do oplot,10^([bestagefix[ip],bestagefix[ip]]-9.),!y.crange,linestyle=1,color=fsc_color('red')
      al_legend,['Lick','Lick@[a/Fe]=0','SPS'],psym=15,color=fsc_Color(['black','red','darkgray']),/right,/top,box=0
      ;plot z
      plot,zthomas,probz,xtitle='[z/H]',position=position2,/nodata,xstyle=4,ystyle=4,/noerase
      maxy = max(!y.crange)
      zsps = findgen(200.)*(max(!x.crange)-min(!x.crange))/200.+min(!x.crange)
      probzsps = gaussian(zsps,[0.8*maxy,spsz,spszerr])
      cgcolorfill,zsps,probzsps,color=fsc_Color('mediumgray')
      plot,zthomas,probz,/noerase,xtitle='[Z/H]',position=position2,ytitle='P([Z/H])/dex'
      vsym,4,/fill,rot=45
      oplot,[spsz],[0.8*maxy],psym=8,color=fsc_color('green')
      
      oplot,[bestz[0],bestz[0]],[0,maxy],linestyle=1,color=fsc_color('black')
      maxprob = max(probzfix)
      oplot,zthomasfix[*,0],probzfix/maxprob*0.9*max(!y.crange),color=fsc_color('red')
      for ip=0,3 do oplot,[bestzfix[ip],bestzfix[ip]],!y.crange,linestyle=1,color=fsc_color('red')

      ;plotalpha
      plot,alphathomas,probalpha,xtitle='[alpha/Fe]',position=position3,/noerase,ytitle='P([alpha/Fe])/dex'
      for ip=1,3 do oplot,[bestalpha[ip],bestalpha[ip]],!y.crange,linestyle=2
      oplot,[bestalpha[0],bestalpha[0]],!y.crange,color=fsc_color('blue'),linestyle=1
      ;plot spec
      plot,spsstr(i).lambda/(1.+spsstr(i).zspec),spsstr(i).contdiv,xrange=[4050,4580],yrange=[0,1.5],position=position4,/noerase,xstyle=1
      oplot,spsstr(i).lambda/(1.+spsstr(i).zspec),spsstr(i).spsspec,color=fsc_color('blue'),thick=2
      for ip=0,15 do oplot,[indpassi[ip],indpassf[ip],indpassf[ip],indpassi[ip],indpassi[ip]],[0.2,0.2,1.3,1.3,0.2],color=fsc_color('cyan')
      for ip=0,15 do if indpassf[ip] lt 4600. then xyouts,indpassi[ip],1.25*(ip mod 2)+0.1,indname[ip],/data,charsize=0.5

      plot,spsstr(i).lambda/(1.+spsstr(i).zspec),spsstr(i).contdiv,xrange=[4600,5400],yrange=[0.5,1.5],position=position5,/noerase,xstyle=1,ystyle=1
      oplot,spsstr(i).lambda/(1.+spsstr(i).zspec),spsstr(i).spsspec,color=fsc_color('blue'),thick=2
      for ip=0,15 do oplot,[indpassi[ip],indpassf[ip],indpassf[ip],indpassi[ip],indpassi[ip]],[0.6,0.6,1.3,1.3,0.2],color=fsc_color('cyan')
      for ip=0,15 do if indpassi[ip] gt 4600. then xyouts,indpassi[ip],0.85*(ip mod 2)+0.55,indname[ip],/data,charsize=0.5

      ;plot chisquare contribution
         ;find where the best values are
      for ic=0,1 do begin
         if ic eq 0 then obj = lickstr(i)
         if ic eq 1 then obj = lickfixstr(i)
         minagediff = min(abs(agearr-obj.best_age[0]),ageloc)
         minzdiff = min(abs(zarr-obj.best_zh[0]),zloc)
         if ic eq 0 then minalphadiff = min(abs(alphaarr-obj.best_alpha[0]),alphaloc)
         if ic eq 1 then minalphadiff = min(abs(alphaarr-obj.fixedalpha),alphaloc)

         index = [obj.cn1[0],obj.cn2[0],obj.ca4227[0],obj.g4300[0],obj.c24668[0],obj.mg1[0],obj.mg2[0],obj.mgb[0],obj.fe4383[0],obj.fe4531[0],obj.fe5270[0],obj.fe5335[0],obj.fe5406[0],obj.hda[0],obj.hga[0],obj.hgf[0]]
         index_err = [obj.cn1[1],obj.cn2[1],obj.ca4227[1],obj.g4300[1],obj.c24668[1],obj.mg1[1],obj.mg2[1],obj.mgb[1],obj.fe4383[1],obj.fe4531[1],obj.fe5270[1],obj.fe5335[1],obj.fe5406[1],obj.hda[1],obj.hga[1],obj.hgf[1]]
         goodind = where(index ne -99 and index_err ne -99,cgoodind) ;goodind of observed spec
         chiarr = fltarr(16)-99.
         if cgoodind gt 0 then begin
            index = index(goodind)
            index_err = index_err(goodind)
            tagmatchnow = tagmatch(goodind)
            for ii=0,cgoodind-1 do begin ;loop over each index
               chiarr(goodind(ii)) = (index(ii)-thomasgrid[ageloc,zloc,alphaloc].(tagmatchnow(ii)))^2/(index_err(ii))^2
            endfor
         endif
         chiarr =chiarr/float(cgoodind)
         vsym,4,/fill,rot=45
         xval = findgen(16)
         colorplot = (ic eq 0) ? 'black':'red'
         if ic eq 0 then plot,xval,chiarr,min_value=min(chiarr(goodind))-5.,xtickv=xval,xticks=15,xtickname=tagthomasout,ytitle='reduced chi square',psym=8,xrange=[-0.5,15.5],xstyle=1,position=position6,yrange=[0,max(chiarr)+0.2],/noerase
         if ic eq 1 then oplot,xval,chiarr,psym=8,color=fsc_color('red')
         cutindex = where(obj.keepindex eq 0,ccutindex)
         if ccutindex gt 0 then oplot,xval(cutindex),chiarr(cutindex),psym=7,symsize=3,color=fsc_color(colorplot)
      endfor
      title=strtrim(mask,2)+' '+strtrim(objname,2)+' '+strtrim(string(slit),2)+' SN='+strtrim(string(spsstr(i).snfit,format='(F4.1)'),2)
      xyouts,0.5,0.96,title,/normal,alignment=0.5
      device,/close
   endif
endfor
endif 
;Find the differences in the measurements of age and metallicity
goodobjs = where(lickstr.good eq 1 and lickstr.goodfspsfit eq 1)

lickfixage = 10^(reform(lickfixstr(goodobjs).best_age[0,*])-9.)
lickfixz = reform(lickfixstr(goodobjs).best_zh[0,*])
lickalpha = reform(lickstr(goodobjs).best_alpha[0,*])

lickfixageerr = reform(10^(lickfixstr(goodobjs).best_age[3,*]-9.)>0.1-10^(lickfixstr(goodobjs).best_age[1,*]-9)<15)/2.
lickfixzerr = reform(lickfixstr(goodobjs).best_zh[3,*]<0.67-lickfixstr(goodobjs).best_zh[1,*]>(-2.25))/2.
lickalphaerr = reform(lickstr(goodobjs).best_alpha[3,*]<0.1-lickstr(goodobjs).best_alpha[1,*]>(-0.3))/2.

spsage = spsstr(goodobjs).age
spsageerr = spsstr(goodobjs).ageerr
spsz = spsstr(goodobjs).feh
spszerr = spsstr(goodobjs).feherr

diffage = lickfixage-spsage
diffageerr = sqrt(lickfixageerr^2+spsageerr^2)
diffz = lickfixz-spsz
diffzerr = sqrt(lickfixzerr^2+spszerr^2)
  
;plot the differences as a function of alpha, mass, signal to noise
!p.multi=[0,1,3]
!p.charsize=1.5

psname=outdir+'diffage_fixedalpha.eps'
device, filename = psname,xsize = 10,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,lickalpha,diffage,lickalphaerr,diffageerr,psym=1,xtitle='[alpha/Fe]',ytitle='lick Age (a=0) - SPS Age (Gyr)',xrange=[-0.3,0.5]
vsym,4,/fill,rot=45
oplot,lickalpha,diffage,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).snfit,diffage,diffageerr,psym=1,xtitle='S/N',ytitle='lick Age (a=0)- SPS Age (Gyr)'
oplot,spsstr(goodobjs).snfit,diffage,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).logmstar,diffage,diffageerr,psym=1,xtitle='log M*',ytitle='lick Age (a=0)- SPS Age (Gyr)',xrange=[9,11.5]
oplot,spsstr(goodobjs).logmstar,diffage,psym=8,color=fsc_color('blue')
device,/close

psname=outdir+'diffz_fixedalpha.eps'
device, filename = psname,xsize = 10,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ploterror,lickalpha,diffz,lickalphaerr,diffzerr,psym=1,xtitle='[alpha/Fe]',ytitle='lick [Z/H] (a=0)- SPS [Z/H]',xrange=[-0.3,0.5]
oplot,lickalpha,diffz,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).snfit,diffz,lickalphaerr,diffzerr,psym=1,xtitle='S/N',ytitle='lick [Z/H] (a=0) - SPS [Z/H]'
oplot,spsstr(goodobjs).snfit,diffz,psym=8,color=fsc_color('blue')
ploterror,spsstr(goodobjs).logmstar,diffz,lickalphaerr,diffzerr,psym=1,xtitle='log M*',ytitle='lick [Z/H] (a=0)- SPS [Z/H]',xrange=[9,11.5]
oplot,spsstr(goodobjs).logmstar,diffz,psym=8,color=fsc_color('blue')
device,/close
stop
end
