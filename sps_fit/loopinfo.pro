pro loopinfo,file
  restore,file
  WINDOW, 0, XSIZE=1000, YSIZE=1100, TITLE='Loop Info'
  !p.font = 0
  !p.multi=[0,1,1]
  !p.charsize=2
  for i=0,4 do begin
     loop = loopstr[i]
     won = loop.won
     yup = 0.95-(0.9/5.)*i
     ydown = yup-(0.9/5.)
     if i eq 0 then noerase=0 else noerase = 1
     plot,loop.lambda(won),loop.contdiv(won)/loop.fakecont(won),position=[0.15,ydown,0.95,yup],/normal,noerase=noerase
     oplot,loop.lambda(won),loop.spsbestfit(won),color=fsc_color('red')
     xyouts,0.03,(yup+ydown)/2.,'loop '+strtrim(string(loop.loop),2),/normal
     xyouts,7000.,0.2,strjoin(string(loop.param),' ')
  endfor
  xyouts,0.5,0.97,file,/normal,alignment=0.5
  ;oplot,loopstr[0].lambda(won),loopstr[0].spsbestfit(won),color=fsc_Color('green')
  stop
end
