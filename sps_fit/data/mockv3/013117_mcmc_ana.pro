pro mcmc_ana,filename
  restore,filename
  nloop = n_elements(str.chisq)
  chisq=str.chisq[1:nloop-1]
  param = str.param[*,1:nloop-1]
  perror = str.perror[*,1:nloop-1]
  paraname= str.paraname
  minchisq=min(chisq,pos)

  if minchisq eq 0 then begin
     chisq = chisq[0:pos-1]
     param = param[*,0:pos-1]
     perror = perror[*,0:pos-1]
     nloop = pos+1
     minchisq=min(chisq,pos)
  endif
  loopnum = indgen(nloop-1)+1

  set_plot,'ps'
  psname = filename+'.eps'
  device, filename = psname,xsize = 30,ysize = 18, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi = [0,1,1]
  !p.font=0
  !p.charsize=0

  plot,loopnum,chisq,yrange=minmax(chisq),position=[0.1,0.75,0.75,0.95]
  oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
  yrange=!y.crange
  plothist,chisq,xrange=yrange,xstyle=5,/autobin,position = [0.75,0.95,0.95,0.75]

  ploterror,loopnum,param[0,*],perror[0,*],psym=1,ytitle='[Fe/H]',position=[0.1,0.55,0.75,0.75],/noerase
  oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
  oplot,!x.crange,[-0.1,-0.1],color=fsc_color('darkgreen')
  yrange=!y.crange
  plothist,param[0,*],xrange=yrange,xstyle=5,/autobin,position = [0.75,0.75,0.95,0.55]

  ploterror,loopnum,param[1,*],perror[0,*],psym=1,ytitle='age',position=[0.1,0.35,0.75,0.55],/noerase
  oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
  oplot,!x.crange,[3,3],color=fsc_color('darkgreen')
  yrange=!y.crange
  plothist,param[1,*],xrange=yrange,xstyle=5,/autobin,position = [0.75,0.55,0.95,0.35]

  ploterror,loopnum,param[2,*],perror[0,*],psym=1,ytitle='Vel disp',position=[0.1,0.15,0.75,0.35],/noerase,xtitle='loop continuum fit'
  oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
  oplot,!x.crange,[250,250],color=fsc_color('darkgreen')
  yrange=!y.crange
  plothist,param[2,*],xrange=yrange,xstyle=5,/autobin,position = [0.75,0.35,0.95,0.15]

  xyouts,0.5,0.97,filename,/normal,alignment=0.5
  xyouts,0.35,0.05,'Best '+strjoin(paraname,' ')+' : '+string(param[0,pos], param[1,pos], param[2,pos],param[3,pos],format='(D6.3,1X,D5.2,2X,D6.1,2x,D6.3)'),/normal,alignment=0.5
  device,/close


  stop
end
