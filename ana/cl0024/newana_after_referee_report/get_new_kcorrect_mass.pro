pro get_new_kcorrect_mass,makefastcat=makefastcat
   sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_cl0024_ana.fits',1)
   cgals = n_elements(sci)
   color = dblarr(7,cgals)
   ;magnitudes are in AB mag
   color[0,*] = sci.B
   color[1,*] = sci.V
   color[2,*] = sci.R
   color[3,*] = sci.I
   color[4,*] = sci.J
   color[5,*] = sci.K
   color[6,*] = sci.F814W
   badcol = where(color eq 99. or color eq 0., cbadcol)
   color[0,*] = sci.B-0.09 
   color[1,*] = sci.V+0.02
   color[2,*] = sci.R+0.21
   color[3,*] = sci.I+0.45
   color[4,*] = sci.J+0.91
   color[5,*] = sci.K+1.85


   color_err = dblarr(7,cgals)      ;b,v,r,i,j,k,f814
   color_err[0,*] = sci.BERR 
   color_err[1,*] = sci.VERR
   color_err[2,*] = sci.RERR
   color_err[3,*] = sci.IERR
   color_err[4,*] = sci.JERR
   color_err[5,*] = sci.KERR
   color_err[6,*] = sci.F814WERR
   color_ivar = 1/color_err^2

   ;fast cat
   if keyword_set(makefastcat) then begin
      badcol = where(color eq 99. or color eq 0,cbadcol)
      flux = 10.^((25.-color)/2.5)
      flux_err = abs(flux*alog(10.)*color_err/2.5)
      if cbadcol gt 0 then flux(badcol) = -99.
      if cbadcol gt 0 then flux_err(badcol) = -99.
      idarr = lindgen(cgals)+1
      writecol,'/scr2/nichal/FAST/FAST_v1.0/leethochawalit17/cl0024.cat',$
        idarr,sigfig(flux[0,*],4,/sci),sigfig(flux_err[0,*],4,/sci),sigfig(flux[1,*],4,/sci),sigfig(flux_err[1,*],4,/sci),$
        sigfig(flux[2,*],4,/sci),sigfig(flux_err[2,*],4,/sci),sigfig(flux[3,*],4,/sci),sigfig(flux_err[3,*],4,/sci),$
        sigfig(flux[4,*],4,/sci),sigfig(flux_err[4,*],4,/sci),sigfig(flux[3,*],4,/sci),sigfig(flux_err[5,*],4,/sci),$
        sigfig(flux[6,*],4,/sci),sigfig(flux_err[6,*],4,/sci),sigfig(sci.zfit,4,/sci),$
        fmt='(I,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A,2x,A)'
      stop
      ;#id F136 E136 F137 E137 F138 E138 F139 E139 F161 E161 F163 E163 F16 E16 z_spec

   endif
   ;finish fast

   filterlist=['cfh12k_B.par','cfh12k_V.par','cfh12k_R.par','cfh12k_I.par','lco_wirc_J.par','lco_wirc_Ks.par','wfpc2_f814w.par']

   zlist = sci.zfit

   kcorrect,color,color_ivar,zlist,kcorrect1,chi2=chi2,filterlist=filterlist,/magnitude,mass=mass,b300=b300
   mass = mass/0.52 ;fix for the hubble constant h=0.7
   mass = alog10(mass)

   datout = create_struct(sci[0],'kcorrect_mass',-999d,'b300_new',-999d)
   datout = replicate(datout,cgals)
   Struct_Assign, sci, datout
   datout.kcorrect_mass = mass
   datout.b300_new = b300
   mwrfits,datout,'sci_cl0024_ana_afterrev.fits',/create,/silent
   ;plot
   window,0,xsize=1000,ysize=400
   !p.multi = [0,2,1]
   !p.charsize = 2.
   plot,sci.logmstar,mass,psym=1,xrange=[8,11],yrange=[8,11],xtitle='old kcorrect',ytitle='new kcorrect',title='Cl0024'
   oplot,[8,11],[8,11],color=fsc_color('red')
   good = where(finite(mass))
   fitparam = linfit(sci(good).logmstar-10.,mass(good),sigma=sigma)
   oplot,!x.crange,(!x.crange-10.)*fitparam(1)+fitparam(0),color=fsc_color('yellow')
   print, fitparam,sigma
   ;read fast cl0024.fout.AB.nosubtract
   readcol,'/scr2/nichal/FAST/FAST_v1.0/leethochawalit17/cl0024.fout',id,z,ltau,metal,lage,Av,lmass_fast,lsfr,lssfr,la2t,chi2,comment='#'
   plot,mass,lmass_fast,psym=1,xrange=[8,11],yrange=[8,11],xtitle='new k correct',ytitle='fast'
   oplot,!x.crange,!x.crange,color=fsc_color('red')
   stop
end
