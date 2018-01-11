pro mkplottest,sci,sdss,ageform,sdss_ageform,fehideal,sdss_fehideal

  psname='FEH_deviation_fromz0line_obs_color_EW.eps'

   scioii = sci.oiiew*(-1.) ;change from negative=emission line to positive = emission line
   rangeoii = [-10,5.]
   color_cl = BytScl(scioii,min=rangeoii[0],max=rangeoii[1],top=23) ;the real range is -10,7
   sdssha = sdss.haew*(-1.)
   rangeha = [-3,1.]
   color_sdss = BytScl(scioii,min=rangeha[0],max=rangeha[1],top=23) ;the real range is -3,2
   device, filename = psname,xsize = 15,ysize = 15, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,sdss_ageform,sdss.feh-sdss_fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle='delta [Fe/H]',xrange=[0,14],xstyle=9,yrange=[-0.4,0.4],/nodata,position=[0.15,0.15,0.95,0.7]
      axis,xaxis=0,xstyle=1,xrange=[0,14]
      axis,xaxis=1,xticks=4,xtickv=[1.558,3.316, 5.903, 8.628,13.712],xtickn=['4','2','1','0.5','0'],xtitle='z'
      cgLoadCT,33,Ncolors=24
      for i=0,23 do begin
          inrange = where(color_cl eq i,cinrange)
          if cinrange gt 0 then cgplot,ageform(inrange),sci(inrange).feh-fehideal(inrange),psym=14,/overplot,color=i
          inrange = where(color_sdss eq i,cinrange)
          if cinrange gt 0 then cgplot,sdss_ageform(inrange),sdss(inrange).feh-sdss_fehideal(inrange),$
                                       psym=16,/overplot,symsize=0.5,color=i
      endfor
      cgColorbar, NColors=24,divisions=5, Range=rangeoii,tickinterval=3, $
        /Discrete,Position=[0.125, 0.86, 0.9, 0.90],charsize=1,ticklen=0.25

      cgColorbar, NColors=24,divisions=4, Range=rangeha,tickinterval=1, $
        /Discrete,Position=[0.125, 0.86, 0.9, 0.90],charsize=1,ticklen=0.25,/top

   device,/close

end
