pro fitchoi,redo=redo
common sps_spec,sps, spsz, spsage

if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)

massbinval = [10.5,10.8,11.0,11.3,10.8,11.1,11.3,10.9,11.0,11.3]
fehchoi    = [-0.11,-0.05,-0.02,-0.03,-0.07,-0.04,-0.05,-0.15,-0.02,-0.05]
fehchoierr = [0.03,0.01,0.01,0.02,0.02,0.01,0.02,0.07,0.03,0.03]
agechoi    = [3.27,3.47,4.55,5.61,2.99,3.28,4.00,2.67,2.49,3.06]
agechoierr = [0.21,0.12,0.13,0.15,0.15,0.13,0.18,0.20,0.12,0.11]
CFe  = [0.19,0.16,0.18,0.19,0.18,0.18,0.14,0.16,0.24,0.26]
NFe  = [0.26,0.25,0.24,0.26,0.26,0.34,0.25,0.18,0.58,0.73]
NFeerr = [0.10,0.05,0.03,0.05,0.08,0.05,0.06,0.30,0.09,0.09]
MgFe = [0.18,0.22,0.21,0.29,0.23,0.23,0.30,0.05,0.09,0.19]
MgFeerr = [0.04,0.02,0.01,0.03,0.04,0.02,0.03,0.13,0.05,0.04]
CaFe = [0.04,0.03,0.04,0.04,0.05,-0.03,0.08,0.06,-0.03,0.00]
set_plot,'ps'

if keyword_set(redo) then begin
   dirdat = '/scr2/nichal/workspace2/Choi14_spec/'
   datfiles = file_search(dirdat+'*.dat',count=cdat)
   loadct,4
   psname = '/scr2/nichal/workspace2/sps_fit/choiout/mpfit/fitchoispec2.eps'
   device, filename = psname,xsize = 30,ysize = 30, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   !p.multi = [0,2,5]
   !p.font = 0
   for i=0,cdat-1 do begin
      fnow  = datfiles[i]
      readcol,fnow,wl,flux,fluxerr,mask
   ;cover CaHK line
   ;   badwl = where(wl lt 4050.,cbadwl)
   ;   if cbadwl ne 0 then mask(badwl) = 0.
   ;cover Mgb line
   ;   badwl = where(wl lt 5200. and wl gt 5150,cbadwl)
   ;   if cbadwl ne 0 then mask(badwl) = 0.
      
      fnow  = strmid(fnow,strpos(fnow,'ages'))
      substr= strsplit(fnow,'_',/extract)
      z     = float((strsplit(substr[1],'0',/extract))[1:2])
      zmean = 0.5*total(z)
      if z[0] eq 0.3  then zbin =0
      if z[0] eq 0.4  then zbin =1
      if z[0] eq 0.55 then zbin =2
      massbin = fix(strmid(substr[2],7,1))-2
      nlamb   = n_elements(wl)
      science={lambda:wl,contdiv:flux,contdivivar:1./fluxerr^2,fitmask:mask,zspec:zmean,spscont:dblarr(nlamb),spsspec:dblarr(nlamb),zfit:zmean,feh:-999.,feherr:-999.,age:-999.,ageerr:-999.,vdisp:-999.,vdisperr:-999,objname:fnow,pcor:dblarr(4,4),covar:dblarr(4,4),agefix:-999.,agefixerr:-999.,fehfix:-999.,fehfixerr:-999.,covarfix:dblarr(4,4),spsspecfix:dblarr(nlamb),spscontfix:dblarr(nlamb),fehgrid:dblarr(4)-999.,agegrid:dblarr(4)-999.,spsspecgrid:dblarr(nlamb),spscontgrid:dblarr(nlamb),agechoi:agechoi[i],agechoierr:agechoierr[i],fehchoi:fehchoi[i],fehchoierr:fehchoierr[i]}
      fitspschoi,science        ;FIT HERE
      fitspschoi_fix,science    ;FIT HERE FIXVELDISP AND Z
      plot,science.lambda,science.contdiv/science.spscont,xtitle='lambda',title=fnow,xrange=minmax(science.lambda),xstyle=1,yrange=[0.4,1.3],ystyle=1,charsize=2
      oplot,science.lambda,science.spsspec,color=120
      masked = 0.42*(1.-science.fitmask)
      oplot,science.lambda,masked,psym=10,thick=15.,color=fsc_color('pink')
      plotlines,0
      mwrfits,science,'/scr2/nichal/workspace2/sps_fit/choiout/'+fnow+'.fits',/create,/silent
   endfor
   device,/close
endif

;make grid search for age and metallicity
   dirsci = '/scr2/nichal/workspace2/sps_fit/choiout/'
   scifiles = file_search(dirsci+'*.fits',count=cfits)
   for i=0,cfits-1 do begin
;   for i=5,5 do begin
      science = mrdfits(scifiles[i],1)
      nstepage = 25
      nstepfe  = 25
      minage = min([science.age-3*science.ageerr,science.agechoi-3*science.agechoierr])
      maxage = max([science.age+3*science.ageerr,science.agechoi+3*science.agechoierr])
      ;minage = 2.5
      rangeage = [minage,maxage]
      minfe    = min([science.feh-3*science.feherr,science.fehchoi-3*science.fehchoierr])
      maxfe    = max([science.feh+3*science.feherr,science.fehchoi+3*science.fehchoierr])
      maxfe    = 0.07
      rangefe  =[minfe,maxfe]
      gridsearch_agefe,science,nstepage,nstepfe,rangeage,rangefe
      mwrfits,science,scifiles[i],/create,/silent
      psname = '/scr2/nichal/workspace2/sps_fit/choiout/grid/'+science.objname+'.gridspec.eps'
      cgerase
      device, filename = psname,xsize = 15,ysize = 15, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.font = 0
      loadct,4
      multiplot,[1,2]
      plot,science.lambda,science.contdiv/science.spscont,xrange=minmax(science.lambda),xstyle=1,yrange=[0.4,1.3],ystyle=1,ytitle='Normalized Flux',title=science.objname,font=0
      oplot,science.lambda,science.spsspec,color=fsc_color('dark green')
      oplot,science.lambda,science.spsspecgrid,color=fsc_color('red')
      masked = 0.42*(1.-science.fitmask)
      oplot,science.lambda,masked,psym=10,thick=15.,color=fsc_color('pink')
      plotlines,0
      al_Legend,['Grid Search','MPFIT SSP'],psym=[0,0],linestyle=[2,2],color=[fsc_color('red'),fsc_color('dark green')],box=0,thick=2,symsize=1.5,/right,/top,font=0
      multiplot
      data = science.contdiv/science.spscont
      plot,science.lambda,100.*(science.spsspec-data)/data,ytitle='residuals(%)',/nodata,xrange=minmax(science.lambda),xstyle=1,yrange=[-10,10],ystyle=1,xtitle='lambda',font=0
      oplot,science.lambda,100.*(science.spsspec-data)/data,color=fsc_color('dark green')
      oplot,science.lambda,100*(science.spsspecgrid-data)/data,color=fsc_color('red')
      oplot,!x.crange,[0,0],linestyle=1
      multiplot,/reset
      device,/close
      cgerase
   endfor

;Retrieve Science
dirsci = '/scr2/nichal/workspace2/sps_fit/choiout/'
scifiles = file_search(dirsci+'*.fits',count=cfits)
namefit = ['MPFIT','Grid_Search']
for jj=0,1 do begin
   fehnicha   = fehchoi*0.
   fehnichaerr= fehnicha
   fehnichahi= fehnicha
   fehnichalo= fehnicha

   agenicha   = fehnicha
   agenichaerr= fehnicha
   agenichahi = fehnicha
   agenichalo = fehnicha

   zarr       = fehnicha
   if jj eq 0 then for i=0,cfits-1 do begin
      science = mrdfits(scifiles[i],1)
      fnow  = strmid(scifiles[i],strpos(scifiles[i],'ages'))
      if science.objname+'.fits' ne fnow then stop
      zarr[i] = science.zspec-0.05
      fehnicha[i] = science.feh
      fehnichaerr[i] = science.feherr
      fehnichahi[i] = science.feh+science.feherr
      fehnichalo[i] = science.feh-science.feherr
      agenicha[i] = science.age
      agenichaerr[i] = science.ageerr
      agenichahi[i] = science.age+science.ageerr
      agenichalo[i] = science.age-science.ageerr
   endfor
   if jj eq 1 then for i=0,cfits-1 do begin
      science = mrdfits(scifiles[i],1)
      fnow  = strmid(scifiles[i],strpos(scifiles[i],'ages'))
      if science.objname+'.fits' ne fnow then stop
      zarr[i]       = science.zspec-0.05
      fehnicha[i]   = science.fehgrid[0]
      fehnichahi[i] = science.fehgrid[3]
      fehnichalo[i] = science.fehgrid[1]
      fehnichaerr[i] = 0.5*(science.fehgrid[3]-science.fehgrid[1])

      agenicha[i]   = science.agegrid[0]
      agenichahi[i] = science.agegrid[3]
      agenichalo[i] = science.agegrid[1]
      agenichaerr[i] = 0.5*(science.agegrid[3]-science.agegrid[1])

   endfor
   psname = namefit[jj]+'_fehageoffset2.eps'
   device, filename = psname,xsize = 20,ysize = 12, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   !p.multi = [0,2,1]
   !X.OMargin = [0, 0]
   !Y.OMargin = [0, 0]
   !p.font =0
   num = strtrim(string(indgen(10)),2)
   labels=num+':'+strtrim(string(fix(zarr),format='(F3.1)'),2)+','+strtrim(string(massbinval,format='(F4.1)'),2)
   labels = strjoin(labels,' ')
   ploterror,fehchoi,fehnicha,fehchoierr,fehnichaerr,xtitle='Feh_Jieun',ytitle='Feh_Nicha',psym=1,position=[0.1,0.2,0.45,0.9],yrange=[-0.2,0.2],ystyle=1,/nodata,color=fsc_color('black'),errcolor=fsc_color('white')
   cgerrplot,fehnicha,fehchoi-fehchoierr,fehchoi+fehchoierr,/Horizontal,width=.05
   cgerrplot,fehchoi,fehnichalo,fehnichahi,width=.05
   oplot,!x.crange,!x.crange,linestyle=2,color=120
   xyouts,fehchoi,fehnicha,num,/data
   ploterror,agechoi,agenicha,agechoierr,agenichaerr,xtitle='Age_Jieun',ytitle='Age_Nicha',psym=1,position=[0.55,0.2,0.95,0.9],yrange=[2.5,4.5],ystyle=1,/nodata,color=fsc_color('black'),errcolor=fsc_color('white')
   cgerrplot,agenicha,agechoi-agechoierr,agechoi+agechoierr,/horizontal,width=.05
   cgerrplot,agechoi,agenichalo,agenichahi,width=.05
   oplot,!x.crange,!x.crange,linestyle=2,color=120
   xyouts,agechoi,agenicha,num,/data
   xyouts,[0.01],[0.05],[labels],/normal
   cgText, 0.5, 0.95, ALIGNMENT=0.5, CHARSIZE=1.25, /NORMAL, namefit[jj]+' Results'
   device,/close

;plot differences in age and Fe as a function of [alpha/Fe], mass
   fehdiff    = fehnicha-fehchoi
   fehdiffhi  = fehdiff+sqrt(fehchoierr^2+(fehnichahi-fehnicha)^2)
   fehdifflo  = fehdiff-sqrt(fehchoierr^2+(fehnichalo-fehnicha)^2)
   agediff    = agenicha-agechoi
   agediffhi  = agediff+sqrt(agechoierr^2+(agenichahi-agenicha)^2)
   agedifflo  = agediff-sqrt(agechoierr^2+(agenichalo-agenicha)^2)
   yfehrange  = minmax([fehdiffhi+0.01,fehdifflo-0.01])
   yagerange  = minmax([agediffhi+0.1,agedifflo-0.1])
   mgferange  = minmax([mgfe-mgfeerr-0.01,mgfe+mgfeerr+0.01])
   Nferange   = minmax([Nfe-Nfeerr-0.01,Nfe+Nfeerr+0.01])
   massrange  = max(massbinval)-min(massbinval)
   masscolor  = fix(100*(massbinval-min(massbinval))/massrange)
   zsym = [15,15,15,15,16,16,16,14,14,14]
   
   cgLoadCT, 33, ncolors=100
   psname = namefit[jj]+'_offsetvsalpha.eps'
   device, filename = psname,xsize = 30,ysize = 15, $
           xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   !X.OMargin = [0, 0]
   !Y.OMargin = [6, 2]
   !p.multi = [0,3,2]
   !p.font  = 0

   cgplot,mgfe,fehdiff,ytitle='$\tex [Fe/H]_{NL}-[Fe/H]_{C14}$',xtitle='[Mg/Fe]',charsize=2,/nodata,xrange=mgferange,yrange=yfehrange
   for i=0,9 do begin
      cgerrplot,fehdiff[i],mgfe[i]+mgfeerr[i],mgfe[i]-mgfeerr[i],color=masscolor[i],/horizontal,width=.05
      cgerrplot,mgfe[i],fehdifflo[i],fehdiffhi[i],color=masscolor[i],width=.05
      cgplot,mgfe[i],fehdiff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgplot,Nfe,fehdiff,ytitle='$\tex [Fe/H]_{NL}-[Fe/H]_{C14}$',xtitle='[N/Fe]',charsize=2,/nodata,xrange=Nferange,yrange=yfehrange
   for i=0,9 do begin
      cgerrplot,fehdiff[i],Nfe[i]+Nfeerr[i],Nfe[i]-Nfeerr[i],/horizontal,color=masscolor[i],width=.05
      cgerrplot,Nfe[i],fehdifflo[i],fehdiffhi[i],color=masscolor[i],width=.05
      cgplot,Nfe[i],fehdiff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgplot,massbinval,fehdiff,ytitle='$\tex [Fe/H]_{NL}-[Fe/H]_{C14}$',xtitle='Mass',charsize=2,/nodata,yrange=yfehrange
   for i=0,9 do begin
      cgerrplot,massbinval[i],fehdifflo[i],fehdiffhi[i],color=masscolor[i],width=.05
      cgplot,massbinval[i],fehdiff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgplot,mgfe,agediff,ytitle='$\tex Age_{NL}-Age_{C14}$',xtitle='[Mg/Fe]',charsize=2,/nodata,xrange=mgferange,yrange=yagerange
   for i=0,9 do begin
      cgerrplot,agediff[i],mgfe[i]+mgfeerr[i],mgfe[i]-mgfeerr[i],/horizontal,color=masscolor[i],width=.05
      cgerrplot,mgfe[i],agedifflo[i],agediffhi[i],color=masscolor[i],width=.05
      cgplot,mgfe[i],agediff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgplot,Nfe,agediff,ytitle='$\tex Age_{NL}-Age_{C14}$',xtitle='[N/Fe]',charsize=2,/nodata,xrange=Nferange,yrange=yagerange
   for i=0,9 do begin
      cgerrplot,agediff[i],Nfe[i]+Nfeerr[i],Nfe[i]-Nfeerr[i],/horizontal,color=masscolor[i],width=.05
      cgerrplot,Nfe[i],agedifflo[i],agediffhi[i],color=masscolor[i],width=.05
      cgplot,Nfe[i],agediff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgplot,massbinval,agediff,ytitle='$\tex Age_{NL}-Age_{C14}$',xtitle='Mass',charsize=2,/nodata,yrange=yagerange
   for i=0,9 do begin
      cgerrplot,massbinval[i],agedifflo[i],agediffhi[i],color=masscolor[i],width=.05
      cgplot,massbinval[i],agediff[i],psym=zsym[i],color=masscolor[i],/overplot
   endfor
   oplot,!x.crange,[0,0],linestyle=2,color=fsc_color('grey')
   
   cgText, 0.5, 0.95, ALIGNMENT=0.5, CHARSIZE=1.25, /NORMAL, namefit[jj]+' Results'
   cgCOLORBAR, NCOLORS=100, RANGE=minmax(massbinval), DIVISIONS=3,position=[0.1, 0.08, 0.4, 0.11],charsize=2,font=0
   cgText,0.25,0.13,charsize=1.25,alignment=0.5,font=0,/normal,'Mass'
   al_Legend,['0.3<z<0.4','0.4<z<0.55','0.55<z<0.7'],psym=[15,16,14],box=0,thick=2,symsize=1.5,font=0,position=[0.6,0.12 ],/normal,/horizontal
   
   device,/close
endfor
stop
end
