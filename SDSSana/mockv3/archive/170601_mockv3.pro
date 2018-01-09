pro mockv3,dirname
  dirraw = '/scr2/nichal/workspace2/mock_data/'+dirname+'/'
  restore,dirraw+'variables.sav'
  realFeH = params_save.z
  realage = params_save.age	

  dir = '/scr2/nichal/workspace2/sps_fit/data/'+dirname+'/'
  filenames = ['sps_fit.fits.gz']
 
  for i=0,n_elements(filenames)-1 do begin
     sci = mrdfits(dir+filenames[i],1,hdr)
     if tag_exist(params_save,'SN') then inputsn = params_save.sn/sqrt(0.65) else inputsn = 10.^(findgen(n_elements(sci))/30.+0.7)/sqrt(0.65)
     set_plot,'ps'

     psname = dirname+'_small.eps'
     device, filename = psname,xsize = 12,ysize = 12, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
     !p.multi = [0,1,1]
     !p.font=0
     !p.charsize=0
     angstrom = cgSymbol("angstrom") 
     ploterror,inputsn,(sci.age-realage)/realage*100,sci.ageerr/realage*100.,psym=1,ytitle='Age inaccuracy (%)',xrange=[0,30],position=[0.2,0.55,0.95,0.95],xstyle=1,xtickformat='(A1)',yrange=[-50,50]
     oplot,!x.crange,[0,0],linestyle=1,thick=2,color=fsc_color('red') 
;     ploterror,inputsn,sci.feh-realfeh,sci.feherr,psym=1,ytitle='[Fe/H] inaccuracy (dex)',xtitle='input S/N ratios',xrange=[0,30],position=[0.2,0.15,0.95,0.55],/noerase,xstyle=1,yticks=6,ytickv=[-0.1,-0.05,0,0.05,0.1,0.15],yrange=[-0.2,0.2]
     ploterror,inputsn,sci.feh-realfeh,sci.feherr,psym=1,ytitle='[Fe/H] inaccuracy (dex)',xtitle='Input S/N ratios ('+angstrom+'!u-1!n)',xrange=[0,30],position=[0.2,0.15,0.95,0.55],/noerase,xstyle=1,yrange=[-0.2,0.2],yticks=7,ytickv=[-0.15,-0.1,-0.05,0.05,0.1,0.15]
     oplot,!x.crange,[0,0],linestyle=1,thick=2,color=fsc_color('red') 

     device,/close
  endfor

  stop
end
