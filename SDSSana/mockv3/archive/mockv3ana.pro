pro mockv3ana
  realFeH = -0.1
  realAGe = 3.
  realVdisp = 250.
  realz     = 0.55
  inputsn = 10.^(findgen(25)/30.+0.7)
  dir = '/scr2/nichal/workspace2/sps_fit/data/mockv3/'
  filenames = ['sps_fit_mcmc.fits.gz']

  for i=0,n_elements(filenames)-1 do begin
     sci = mrdfits(dir+filenames[i],1,hdr)
  
     set_plot,'ps'
     psname = 'mockv3_mcmc.eps'
     device, filename = psname,xsize = 25,ysize = 14, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
     !p.multi = [0,3,2]
     !p.font=0
     !p.charsize=2
     
                                ;SN
     plot,inputsn,sci.sn,ytitle='measured SN',xtitle='input SN'
     oplot,[0,35],[0,35],linestyle=1,color=fsc_color('red')
                                ;Age
     ploterror,sci.sn,sci.age,sci.ageerr,psym=1,ytitle='Age (Gyr)',xtitle='measured SN'
     oplot,!x.crange,[realAge,realAge],linestyle=1,thick=2,color=fsc_color('red') 
                                ;vdisp
     ploterror,sci.sn,sci.vdisp,sci.vdisperr,psym=1,ytitle='Vdisp (Gyr)',xtitle='measured SN'
     oplot,!x.crange,[realVdisp,realVdisp],linestyle=1,thick=2,color=fsc_color('red') 
                                ;FeH
     ploterror,sci.sn,sci.feh,sci.feherr,psym=1,ytitle='[Fe/H] (all wl)',xtitle='measured SN'
     oplot,!x.crange,[realFeh,realFeh],linestyle=1,thick=2,color=fsc_color('red') 
                                ;FeH from Fe bands
     ploterror,sci.sn,sci.fe,sci.feerr,psym=1,ytitle='[Fe/H] (Fe bands)',xtitle='measured SN'
     oplot,!x.crange,[realFeH,realFeH],linestyle=1,thick=2,color=fsc_color('red') 
                                ;FeH from Mg bands
     ploterror,sci.sn,sci.mg,sci.mgerr,psym=1,ytitle='[Fe/H] (Mg bands)',xtitle='measured SN'
     oplot,!x.crange,[realFeH,realFeH],linestyle=1,thick=2,color=fsc_color('red') 
  
     device,/close
  endfor

  stop
end
