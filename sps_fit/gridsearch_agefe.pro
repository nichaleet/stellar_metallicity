pro gridsearch_agefe,science,nstepage,nstepfe,rangeage,rangefe
  common get_sps, dlam
  clight = 299792.458
  mask  = science.fitmask
  vdisp = science.vdisp
  zfit  = science.zfit
  pi = [science.feh,science.age,vdisp,zfit]
  restlambda = science.lambda
  reallambda = science.lambda*(1.+science.zfit)

  agestepsize = (rangeage[1]-rangeage[0])/nstepage
  festepsize  = (rangefe[1]-rangefe[0])/nstepfe
  agearr= rangeage[0]+agestepsize*findgen(nstepfe)
  fearr = rangefe[0]+festepsize*findgen(nstepfe)
  chisqarr = dblarr(nstepage,nstepfe)-99.

  won   = where(mask eq 1,cwon)
  xin   = reallambda[won]
  yin   = science.contdiv[won]/science.spscont[won]
  yinerr  = (science.contdivivar[won])^(-0.5) / science.spscont[won]
  fediff = 1.
  agediff = 1.
  print, '* * * * * * * * * * * * * * * * * * * *'
  print, strtrim(science.objname, 2)
  print, '* * * * * * * * * * * * * * * * * * * *'
  print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
  print, '--- ------- ----- ------- ---------  -------- ----'
  count=0
  ;set_plot,'x'
  dlam = reallambda*350./clight
  nlam = float(n_elements(won))
  dof  = nlam-2.-1.
  for ai = 0, nstepage-1 do for fi=0,nstepfe-1 do begin
     
     spsnow = get_sps_obs(xin,[fearr[fi],agearr[ai],vdisp,zfit])
     chisqarr[ai,fi] = total((yin-spsnow)^2/yinerr^2)
     if fi mod 10 eq 0 then begin
        print, ai,fi
        print, 'age:'+strtrim(string(agearr[ai]),2)+' Fe/H:'+strtrim(string(fearr[fi]),2)+'Reduced Chisq:'+strtrim(string(chisqarr[ai,fi]/dof),2)
                                ;plot,xin,yin
                                ;oplot,xin,spsnow,color=fsc_Color('green')
     endif
  endfor

  Likelihoodarr = -0.5*chisqarr
  likelihoodarr = likelihoodarr-max(likelihoodarr)
  probarr  = exp(likelihoodarr)
  
                                ;normallize
  agearr2d = rebin(agearr,nstepage,nstepfe)
  fearr2d  = rebin(transpose(fearr),nstepage,nstepfe)
  volume   = int_tabulated_2d(agearr2d,fearr2d,probarr)
  probarr  = probarr/volume
 ; stop
                                ;marginalize
  p_age = dblarr(nstepage)
  p_fe  =  dblarr(nstepfe)
  for i=0,nstepage-1 do p_age(i) = int_tabulated(fearr,probarr[i,*])
  for i=0,nstepfe-1 do p_fe(i)= int_tabulated(agearr,probarr[*,i])
  agegrid = confidence_interval(agearr,p_age)
  fehgrid  = confidence_interval(fearr,p_fe) 
                                ;peak,lowerbound,median,upperbound

  dlam = restlambda*350./clight
  spsbestfit = get_sps_rest(restlambda,pi)

  ;save best results
  science.spsspecgrid = spsbestfit
  science.spscontgrid = science.spscont
  science.agegrid = agegrid
  science.fehgrid = fehgrid
  ;plot
  set_plot,'ps'
  psname='/scr2/nichal/workspace2/sps_fit/choiout/grid/'+science.objname+'_grid.eps'
  device, filename = psname,xsize = 30,ysize = 10, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  cgLoadCT, 33, CLIP=[30,255]
  !p.multi=[0,3,1]
  position =   [0.2, 0.125, 0.9, 0.85]
  cgImage, probarr, Stretch=1, MinValue=min(probarr), MaxValue=max(probarr), $
           /Axes, XTitle='Age(Gyr)',YTitle='[fe/H]',title=science.objname, Position=position,XRange=rangeage, YRange=rangefe,font=0
  ell = covar2ellipse(science.covar[0:1,0:1],nsigma=1.)
  tvellipse,ell.minor,ell.major,science.age,science.feh,180-ell.angle,linestyle=1,thick=2,color=fsc_color('yellow'),/data
  ell = covar2ellipse(science.covar[0:1,0:1],nsigma=2.)
  tvellipse,ell.minor,ell.major,science.age,science.feh,180-ell.angle,linestyle=1,thick=2,color=fsc_color('yellow'),/data
  ell = covar2ellipse(science.covarfix[0:1,0:1],nsigma=1.)
  tvellipse,ell.minor,ell.major,science.agefix,science.fehfix,180-ell.angle,linestyle=1,thick=2,color=fsc_color('red'),/data
  ell = covar2ellipse(science.covarfix[0:1,0:1],nsigma=2.)
  tvellipse,ell.minor,ell.major,science.agefix,science.fehfix,180-ell.angle,linestyle=1,thick=2,color=fsc_color('red'),/data

  oploterror,[science.age],[science.feh],[science.ageerr],[science.feherr],color=fsc_color('yellow'),thick=2,errcolor=fsc_color('yellow')
  oploterror,[science.agefix],[science.fehfix],[science.agefixerr],[science.fehfixerr],color=fsc_color('red'),thick=2,errcolor=fsc_color('red')
  cgplot,agearr,p_age,xtitle='Age(Gyr)',ytitle='P(Age)',font=0
  mpfitage = exp(-0.5*(agearr-science.agefix)^2/(science.agefixerr)^2)/(sqrt(2.*!pi)*science.agefixerr)
  oplot,agearr,mpfitage,linestyle=2,color=fsc_color('dark green')
  agechoi = exp(-0.5*(agearr-science.agechoi)^2/(science.agechoierr)^2)/(sqrt(2.*!pi)*science.agechoierr)
  oplot,agearr,agechoi,linestyle=2,color=fsc_color('red')
  vline,agegrid,color=fsc_color('pink')
  pcormpfit = science.pcor[0,1]
  xyouts,rangeage[0]+0.05,!y.crange[1]*0.9,'MPFIT correlation: '+strtrim(string(pcormpfit,format='(F5.2)'),2),/data


  cgplot,fearr,p_fe,xtitle='[Fe/H]',ytitle='P([Fe/H])',font=0,color=fsc_color('black')
  mpfitfeh = exp(-0.5*(fearr-science.fehfix)^2/(science.fehfixerr)^2)/(sqrt(2.*!pi)*science.fehfixerr)
  oplot,fearr,mpfitfeh,linestyle=2,color=fsc_color('dark green')
  fehchoi = exp(-0.5*(fearr-science.fehchoi)^2/(science.fehchoierr)^2)/(sqrt(2.*!pi)*science.fehchoierr)
  oplot,fearr,fehchoi,linestyle=2,color=fsc_color('red')
  vline,fehgrid,color=fsc_color('pink')
  al_Legend,['Grid Search','MPFIT SSP','Choi14'],psym=[0,0,0],linestyle=[0,2,2],color=[fsc_color('black'),fsc_color('dark green'),fsc_color('red')],box=0,thick=2,symsize=1.5,/right,/top,font=0

  device,/close
end
