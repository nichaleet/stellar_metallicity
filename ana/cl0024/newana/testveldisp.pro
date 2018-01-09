pro testveldisp

   sci = mrdfits('sci_cl0024_ana.fits',1)
   probsci = mrdfits('cl0024_feh_age_probdist.fits',1)
   cat =mrdfits('moran05cat.fits',1)
   ngals=n_elements(sci)
   matchsci = []
   matchcat = []
   for i=0,ngals-1 do begin 
      aa = strmatch(cat.name,strtrim(sci(i).objname,2))
      good = where(aa eq 1)
      if good[0] ne -1 then begin
	matchsci = [matchsci,i]
        matchcat = [matchcat,good]
      endif
   endfor

   sci = sci(matchsci)
   cat = cat(matchcat)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   set_plot,'ps'
   !p.font=0
   psname='testvdisp.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      plot,cat.sigma0,sci.vdisp,/nodata,xtitle='moran05 Vdisp',ytitle='NL17 Vdisp',psym=1
      oploterror,cat.sigma0,sci.vdisp,cat.E_sigma0,sci.vdisperr,psym=1
      cgplot,cat.sigma0,sci.vdisp,psym=14,/overplot
      oplot,[0,500],[0,500],color=fsc_Color('red'),linestyle=2
   device,/close
   stop
end
