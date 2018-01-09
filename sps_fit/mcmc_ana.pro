pro mcmc_ana,xval,xname,tryname
  realFeH = -0.1
  realAGe = 3.
  realVdisp = 250.
  realz     = 0.55

  filenames = file_search('fitdata*.sav',count=cfile)
  ;rearrange the files according to the name (1-24)
  length=strpos(filenames,'.')-8
  filenum = intarr(cfile)
  for i=0,cfile-1 do filenum[i] = fix(strmid(filenames[i],8,length[i]))
  filenames = filenames(sort(filenum))
  filenum = filenum(sort(filenum))

  bestchisqarr = fltarr(cfile)
  bestparam = fltarr(4,cfile)
  bestparamerr = fltarr(4,cfile)
  medparam = fltarr(4,cfile)
  medparamerr = fltarr(4,cfile)

  for i=0,cfile-1 do begin
     restore,filenames[i]
     nloop = n_elements(str.chisq)
     chisq=str.chisq
     param = str.param
     perror = str.perror
     paraname= str.paraname
     dof = str.dof
     chisq = chisq*dof
     
     minchisq=min(chisq,pos)
     if minchisq eq 0 then begin
        chisq = chisq[0:pos-1]
        param = param[*,0:pos-1]
        perror = perror[*,0:pos-1]
        nloop = pos
        minchisq=min(chisq,pos)
     endif
     loopnum = indgen(nloop)+1
     
     goodchisq = where(chisq ge minchisq-1. and chisq le minchisq+1.,cgood)
     
     file_mkdir,'mcmcana_out'
     set_plot,'ps'
     psname = 'mcmcana_out/'+filenames[i]+'.eps'
     device, filename = psname,xsize = 30,ysize = 18, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
     !p.multi = [0,1,1]
     !p.font=0
     !p.charsize=0
     
     plot,loopnum,chisq/dof,yrange=minmax(chisq[1:nloop-1])/dof,position=[0.1,0.75,0.8,0.95],ytitle='chisq'
     vsym,4,rot=45
     oplot,!x.crange,([minchisq,minchisq]-1.)/dof,color=fsc_color('red'),linestyle=2
     oplot,!x.crange,([minchisq,minchisq]+1.)/dof,color=fsc_color('red'),linestyle=2
     oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
     yrange=!y.crange
     xrange=[0,cgood]
     plothist,chisq(goodchisq)/dof,yrange=yrange,bin=0.25/dof,position = [0.8,0.75,0.95,0.95],/rotate,ytickformat="(A1)",xrange=xrange,xtickformat="(A1)",/noerase
     
     ploterror,loopnum,param[0,*],perror[0,*],psym=1,ytitle='[Fe/H]',position=[0.1,0.55,0.8,0.75],/noerase,xtickformat="(A1)"
     oplot,loopnum(goodchisq),param[0,goodchisq],psym=8,color=fsc_color('red')
     oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
     oplot,!x.crange,[-0.1,-0.1],color=fsc_color('darkgreen')
     yrange=!y.crange
     plothist,param[0,goodchisq],yrange=yrange,bin=0.025,position = [0.8,0.55,0.95,0.75],/rotate,ytickformat="(A1)",xrange=xrange,xtickformat="(A1)",/noerase
     
     ploterror,loopnum,param[1,*],perror[0,*],psym=1,ytitle='age',position=[0.1,0.35,0.8,0.55],/noerase,xtickformat="(A1)"
     oplot,loopnum(goodchisq),param[1,goodchisq],psym=8,color=fsc_color('red')
     oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
     oplot,!x.crange,[3,3],color=fsc_color('darkgreen')
     yrange=!y.crange
     plothist,param[1,goodchisq],yrange=yrange,bin=0.25,position = [0.8,0.35,0.95,0.55],/rotate,ytickformat="(A1)",xrange=xrange,xtickformat="(A1)",/noerase
     
     ploterror,loopnum,param[2,*],perror[0,*],psym=1,ytitle='Vel disp',position=[0.1,0.15,0.8,0.35],/noerase,xtitle='loop continuum fit'
     oplot,loopnum(goodchisq),param[2,goodchisq],psym=8,color=fsc_color('red')
     oplot,[loopnum(pos),loopnum(pos)],!y.crange,linestyle=1,color=fsc_color('red')
     oplot,!x.crange,[250,250],color=fsc_color('darkgreen')
     yrange=!y.crange
     plothist,param[2,goodchisq],yrange=yrange,bin=5,position = [0.8,0.15,0.95,0.35],/rotate,ytickformat="(A1)",xrange=xrange,/noerase,charsize=0
     
     xyouts,0.5,0.97,filenames[i],/normal,alignment=0.5
     xyouts,0.15,0.06,'min(chisq)                   : '+strjoin(paraname,', ')+' : '+string(param[0,pos], param[1,pos], param[2,pos],param[3,pos],format='(D6.3,1X,D5.2,2X,D6.1,2x,D6.3)'),/normal
     medp=[median(param[0,goodchisq]), median(param[1,goodchisq]), median(param[2,goodchisq]),median(param[3,goodchisq])]
     if cgood gt 1 then medperr=[stdev(param[0,goodchisq]), stdev(param[1,goodchisq]), stdev(param[2,goodchisq]),stdev(param[3,goodchisq])] else medperr=fltarr(4)

     xyouts,0.15,0.03,'Median(min(chisq)+-1): '+strjoin(paraname,', ')+' : '+string(medp[0],medp[1],medp[2],medp[3],format='(D6.3,1X,D5.2,2X,D6.1,2x,D6.3)'),/normal
     restore,'obj'+strtrim(string(filenum[i]),2)+'loopinfo.sav'
     xyouts,0.15,0,'initial guess      :'+strjoin(paraname,', ')+' : '+string(firstguess[0],firstguess[1],firstguess[2],firstguess[3],format='(D6.3,1X,D5.2,2X,D6.1,2x,D6.3)'),/normal

     bestchisqarr[i] = minchisq/dof
     bestparam[*,i] = param[*,pos]
     bestparamerr[*,i] = perror[*,pos]
     medparam[*,i] = medp
     medparamerr[*,i] = medperr
     
     device,/close
  endfor

  psname='mcmcana_out/all_objs.eps'
  device, filename = psname,xsize = 20,ysize = 18, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
                                ;chisq
  plot,xval,bestchisqarr,ytitle='best chisq',position=[0.1,0.75,0.9,0.95],xtickformat="(A1)",yrange=minmax(bestchisqarr)
                                ;Age
  ploterror,xval,bestparam[1,*],bestparamerr[1,*],ytitle='Age (Gyr)',position=[0.1,0.55,0.9,0.75],/noerase,xtickformat="(A1)"
  oploterror,xval,medparam[1,*],medparamerr[1,*],color=fsc_color('blue'),errcolor=fsc_color('blue')
  oplot,!x.crange,[realAge,realAge                                                                                                                                                                                   ],linestyle=1,thick=2,color=fsc_color('red') 
                                ;vdisp
  ploterror,xval,bestparam[2,*],bestparamerr[2,*],ytitle='vdisp',position=[0.1,0.35,0.9,0.55],/noerase,xtickformat="(A1)"
  oploterror,xval,medparam[2,*],medparamerr[2,*],color=fsc_color('blue'),errcolor=fsc_color('blue')
  oplot,!x.crange,[realVdisp,realVdisp],linestyle=1,thick=2,color=fsc_color('red') 
                                ;FeH
  ploterror,xval,bestparam[0,*],bestparamerr[0,*],ytitle='[Fe/H] (all wl)',position=[0.1,0.15,0.9,0.35],/noerase,xtitle=xname
  oploterror,xval,medparam[0,*],medparamerr[0,*],color=fsc_color('blue'),errcolor=fsc_color('blue')
  oplot,!x.crange,[realFeh,realFeh],linestyle=1,thick=2,color=fsc_color('red') 
  xyouts,0.5,0.05,tryname,alignment=0.5,/normal
  device,/close
  


  stop
end
