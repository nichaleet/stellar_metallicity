pro myplan_onchoi


set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
psname='myplan_onchoi.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

zcat = [0.1,0.2,0.3,0.4,0.55,0.7,1.0]
color  = ['deep pink','gold','orange red','dark green','blue','purple']

vsym,5,/star
z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
z04 = {zlow:0.3,zhigh:0.4,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
z06 = {zlow:0.3,zhigh:0.4,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}

plot,z01.mass,z01.feh,xtitle='Log(M/Msun)',ytitle='[Fe/H]',/nodata,xrange=[9.5,11.5],xstyle=5,yrange=[-.8,0.3],ystyle=5

oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color(color(0)),linethick=2,errcolor=fsc_color(color(0))
oplot,z01.mass,z01.feh,psym=8,color=fsc_color(color(0)),symsize=1.5
oploterror,z02.mass,z02.feh,z02.feherr,color=fsc_color(color(1)),linethick=1.5,errcolor=fsc_color(color(1))
oplot,z02.mass,z02.feh,psym=8,color=fsc_color(color(1)),symsize=1.5
oploterror,z04.mass,z04.feh,z04.feherr,color=fsc_color(color(3)),linethick=1.5,errcolor=fsc_color(color(3))
oplot,z04.mass,z04.feh,psym=8,color=fsc_color(color(3)),symsize=1.5
oploterror,z06.mass,z06.feh,z06.feherr,color=fsc_color(color(4)),linethick=1.5,errcolor=fsc_color(color(4))
oplot,z06.mass,z06.feh,psym=8,color=fsc_color(color(4)),symsize=1.5

oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color(color(2)),linethick=1.5,errcolor=fsc_color(color(2))
oplot,z03.mass,z03.feh,psym=8,color=fsc_color(color(2)),symsize=1.5

;add my plan

;Labelling
zarr_str = strarr(n_elements(zcat)-1)
for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
vsym,5,/star
al_Legend,[zarr_str[0:5],'Choi et al. (2014)'],psym=[15,15,15,15,15,15,8],color=[color[0:5],'black'],box=0,thick=2,charsize=1,symsize=[1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5],/right,/bottom,font=0
plot,z01.mass,z01.feh,xtitle='Log(M!D*!N/M!Isun!N)',ytitle='[Fe/H]',/nodata,xrange=[9.5,11.5],xstyle=1,yrange=[-1.,0.3],ystyle=1,/noerase

device,/close
stop

end
