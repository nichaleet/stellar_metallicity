pro spec_combination_ana
   ;get science fits
   sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/combined_uniq_ms0451/sps_fit.fits.gz',1)
   nspec = n_elements(sci)
   set_plot,'ps'
   psname='snr_spec_combination.eps'
        device, filename = psname,xsize = 30,ysize = 10, $
                        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,[0,nspec],[0,20],/nodata,xtitle='spec number',ytitle='SNR'
   for i=0,nspec-1 do begin
      ;read original files
      specnow = sci[i]
      oplot,[i],[specnow.sn],psym=cgsymcat(46),color=fsc_color('red')
      ori = mrdfits(strtrim(specnow.spec1dfile,2),2)
      nfiles = n_elements(ori)
      oplot,fltarr(nfiles)+i,fltarr(nfiles)+ori.sn,psym=cgsymcat(15),color=fsc_color('blue')
   endfor
   device,/close
end
