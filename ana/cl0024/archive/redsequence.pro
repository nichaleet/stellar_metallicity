pro gauss_emission,x,A,F,pder
  if A[0] gt 0 and A[4] lt 0 and A[4] gt -1. and A[2] lt 700./300000.*A[1] and A[5] lt 400./300000.*A[1] and a[3] lt 0.1 and a[3] gt 0. then begin
     emission = A[0]*EXP(-0.5*((x-A[1])/A[2])^2)+A[3]
     absorption= A[4]*EXP(-0.5*((x-A[1])/A[5])^2)+1.-A[3]
     F = emission+absorption
  endif else F=fltarr(n_elements(x))-99.
  
  if n_params() ge 4 then $
     pder = [[EXP(-0.5*((x-A[1])/A[2])^2)],[(emission-A[3])*(x-A[1])/(A[2])^2],[(emission-A[3])*(x-A[1])^2/(A[2])^3],[replicate(1.0, N_ELEMENTS(X))],[EXP(-0.5*((x-A[1])/A[5])^2)],[A[4]*EXP(-0.5*((x-A[1])/A[5])^2)*(x-A[1])^2/(A[5])^3]]
end

pro redsequence
;1) red sequence
;2) ssfr
;3) UVJ
outdir = '/scr2/nichal/workspace2/ana/cl0024/redsequence_output/'
;read in sps_measurement
spsfile = '/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz'
spsstr = mrdfits(spsfile,1,/silent)

;kcorrect the colors
absmagfile=file_search('/scr2/nichal/workspace2/ana/cl0024/cl0024_absmag_kcorrected.sav',count=cabsfile) 

if cabsfile eq 0 then begin
   goodobj = where(spsstr.good eq 1 and spsstr.goodfit eq 1, cgals)
   absmagarr = fltarr(7,cgals) ;output ;U,B,V,R,I,J,K
   if cgals gt 0 then begin
      zarr  = spsstr[goodobj].zfit
      color = dblarr(7,cgals)   ;b,v,r,i,j,k,f814
      ;;assume magnitudes were in vega system (Treu2003). 
      ;;Converting to AB using http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
      color[0,*] = spsstr[goodobj].b-0.09
      color[1,*] = spsstr[goodobj].v+0.02
      color[2,*] = spsstr[goodobj].r+0.21
      color[3,*] = spsstr[goodobj].i+0.45
      color[4,*] = spsstr[goodobj].j+0.91
      color[5,*] = spsstr[goodobj].k+1.85
      ;; convert to AB mag using http://www.stsci.edu/hst/wfc3/phot_zp_lbn using 0.4" aperture
      color[6,*] = spsstr[goodobj].f814w+0.4182
      
      color_err = dblarr(7,cgals) ;b,v,r,i,j,k,f606,f814
      color_err[0,*] = spsstr[goodobj].berr
      color_err[1,*] = spsstr[goodobj].verr
      color_err[2,*] = spsstr[goodobj].rerr
      color_err[3,*] = spsstr[goodobj].ierr
      color_err[4,*] = spsstr[goodobj].jerr
      color_err[5,*] = spsstr[goodobj].kerr
      color_err[6,*] = spsstr[goodobj].f814werr
      color_ivar = 1/color_err^2
      nanpos = where(~finite(color_ivar),cnan)
      if cnan gt 0 then color_ivar(nanpos) = -999d
      
      mask = spsstr(goodobj).mask
      objname = spsstr(goodobj).objname
      slit = spsstr(goodobj).slit
      ;;somehow the kcorrect code stop working on drum. Save the parameters to run on Mahler 
      save,filterlist,color,color_ivar,zarr,mask,objname,slit,filename='redsequenceformahler.sav'
      ;go to mahler /raid/nichal/redsequence.pro

   endif

endif else restore,absmagfile ;mask,objname,slit,absmagarr,bands
   
set_plot,'x'
!p.multi = [0,1,2]
!p.font=0
!p.charsize=1
;;Calculate EW(Hb), EW(OIII)
goodobj = where(spsstr.good eq 1 and spsstr.goodfit eq 1, cgals)
EWHb= fltarr(cgals)
EWOIII = fltarr(cgals)

for i=0,cgals-1 do begin
   contdiv = spsstr(goodobj(i)).contdiv
   contdivivar = spsstr(goodobj(i)).contdivivar
   lambda = spsstr(goodobj(i)).lambda/(1.+spsstr(goodobj(i)).zspec)
   Hbreg = where(lambda gt 4830 and lambda lt 4890)
   widehbreg = [-1.*reverse(indgen(5)+1)+min(hbreg),hbreg,indgen(5)+1+max(hbreg)]
   a = [0.5,4862,2.,0.05,-0.3,3]
   if max(contdiv(hbreg)) gt 2. then a[0]=max(contdiv(hbreg))-1
   fita = [1,1,1,1,1,1]
   hbfit = curvefit(lambda(hbreg),contdiv(hbreg),contdivivar(hbreg),a,fita=fita,function_name='gauss_emission',itmax=1e5,status=status)
   ewhb(i) = int_tabulated(lambda(hbreg),1.-hbfit)
   plot,lambda(widehbreg),contdiv(widehbreg),psym=10
   oplot,lambda(hbreg),hbfit,color=fsc_color('orange')
   hba = a
   if a[0] gt 0.05 and abs(a[1]-4862) lt 2. then begin
      emission = A[0]*EXP(-0.5*((lambda(hbreg)-A[1])/A[2])^2)+A[3]
      absorption = A[4]*EXP(-0.5*((lambda(hbreg)-A[1])/A[5])^2)+1.-A[3]
      oplot,lambda(hbreg),emission,color=fsc_color('green')            
      oplot,lambda(hbreg),absorption,color=fsc_color('green')            
      ;EWhb(i) = sqrt(2)*a[0]*a[2]*sqrt(!pi)
      xyouts,0.2,0.7,'EW = '+strtrim(string(EWhb(i),format='(F4.1)'),2),/normal
   endif
   
   
   OIIIreg = where(lambda gt 5000 and lambda lt 5016)
   wideOIIIreg = [-1.*reverse(indgen(5)+1)+min(OIIIreg),OIIIreg,indgen(5)+1+max(OIIIreg)]
   a = [0.5,5009,2.,0.05,-0.3,3]
   if max(contdiv(oiiireg)) gt 2. then a[0]=max(contdiv(oiiireg))-1
   fita = [1,1,1,1,1,1]
   OIIIfit = curvefit(lambda(OIIIreg),contdiv(OIIIreg),contdivivar(OIIIreg),a,fita=fita,function_name='gauss_emission',itmax=1e5,status=status)
   ewoiii(i) = int_tabulated(lambda(oiiireg),1.-oiiifit)
   plot,lambda(wideOIIIreg),contdiv(wideOIIIreg),psym=10
   oplot,lambda(OIIIreg),oIIIfit,color=fsc_color('orange')

   OIIIa = a
   if a[0] gt 0.05 and abs(a[1]-5009) lt 2. then begin
      emission = A[0]*EXP(-0.5*((lambda(OIIIreg)-A[1])/A[2])^2)+A[3]
      absorption = A[4]*EXP(-0.5*((lambda(OIIIreg)-A[1])/A[5])^2)+1.-A[3]
      oplot,lambda(OIIIreg),emission,color=fsc_color('green')            
      oplot,lambda(OIIIreg),absorption,color=fsc_color('green')            
      xyouts,0.2,0.3,'EW = '+strtrim(string(EWoiii(i),format='(F4.1)'),2),/normal
   endif
   
   ;wait,1
endfor

;;PLOT

SF = where(ewhb lt -1.,csf)
set_plot,'ps'
!p.multi=[0,2,2]
psname=outdir+'cl0024_redsequence.eps'
device, filename = psname,xsize = 20,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;;1) Red sequence
;;absmagarr are in UBVRIJK
vsym,4,rot=25,/fill
plot,absmagarr[3,*],absmagarr[1,*]-absmagarr[3,*],yrange=[0,2],psym=8,xtitle='R',ytitle='B-R',title='Cl0024';B-R vs R
oplot,absmagarr[3,SF],absmagarr[1,SF]-absmagarr[3,SF],color=fsc_color('blue'),psym=8
al_legend,['EW(Hb)<-1'],psym=8,color=fsc_color('blue'),/bottom,/right
;;2) UVJ diagram
plot,absmagarr[2,*]-absmagarr[5,*],absmagarr[0,*]-absmagarr[2,*],yrange=[0,2.5],xrange=[0,2],psym=8,xtitle='V-J',ytitle='U-V'
oplot,absmagarr[2,SF]-absmagarr[5,SF],absmagarr[0,SF]-absmagarr[2,SF],color=fsc_color('blue'),psym=8
oplot,[0,0.75,1.5,1.5],[1.3,1.3,1.9,2.5],color=fsc_color('red')
;;3) Hb EW
plothist,ewhb,bin=0.3,xtitle='EW(Hb)'
;;3) Hb vs OIII eq
sfsure = where(ewhb lt -0. and ewoiii lt 0.,csure)
kewley_sii_ha = findgen(100)/35.-2.6
kewley_oiii_hb=0.72/(kewley_sii_ha-0.32)
xrange=[-2.,0.5]
plot,kewley_sii_ha,kewley_oiii_hb,xtitle='log([SII]/Ha)',ytitle='log([OIII]/Hb)',/nodata,xrange=xrange,xstyle=1,yrange=[-1.5,1.5],ystyle=1
oplot,kewley_sii_ha,kewley_oiii_hb,color=fsc_color('red')
oiii_hb = alog10(ewoiii(sfsure)/ewhb(sfsure))
xarr = findgen(csure)/csure*(xrange[1]-xrange[0]-0.2)+xrange[0]+0.1
oplot,xarr,oiii_hb,psym=8,color=fsc_color('cyan')
device,/close
stop
end
