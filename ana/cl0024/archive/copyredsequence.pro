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
   absmagarr = fltarr(8,cgals)+99.
   if cgals gt 0 then begin
      color = dblarr(8,cgals)   ;b,v,r,i,j,k,f606,f814
      ;;assume magnitudes were in vega system (Treu2003). 
      ;;Converting to AB using http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
      color[0,*] = spsstr[goodobj].b-0.09
      color[1,*] = spsstr[goodobj].v+0.02
      color[2,*] = spsstr[goodobj].r+0.21
      color[3,*] = spsstr[goodobj].i+0.45
      color[4,*] = spsstr[goodobj].j+0.91
      color[5,*] = spsstr[goodobj].k+1.85
      ;; convert to AB mag using http://www.stsci.edu/hst/wfc3/phot_zp_lbn using 0.4" aperture
      color[6,*] = spsstr[goodobj].f606w+0.0825
      color[7,*] = spsstr[goodobj].f814w+0.4182
      
      color_err = dblarr(8,cgals) ;b,v,r,i,j,k,f606,f814
      color_err[0,*] = spsstr[goodobj].berr
      color_err[1,*] = spsstr[goodobj].verr
      color_err[2,*] = spsstr[goodobj].rerr
      color_err[3,*] = spsstr[goodobj].ierr
      color_err[4,*] = spsstr[goodobj].jerr
      color_err[5,*] = spsstr[goodobj].kerr
      color_err[6,*] = spsstr[goodobj].f606werr
      color_err[7,*] = spsstr[goodobj].f814werr
      color_ivar = 1/color_err^2
      nanpos = where(~finite(color_ivar),cnan)
      if cnan gt 0 then color_ivar(nanpos) = -999d
      
      filterlist=['cfh12k_B.par','cfh12k_V.par','cfh12k_R.par','cfh12k_I.par','lco_wirc_J.par','lco_wirc_Ks.par','wfpc2_f606w.par','wfpc2_f814w.par']
      massarr = dblarr(cgals)
      b300arr = dblarr(cgals)-999
   
      for i=0,cgals-1 do begin
         colornow = color[*,i]
         color_ivarnow = color_ivar[*,i]  
         filternow = filterlist
         badpos = where(colornow gt 50. or colornow lt 5 or color_ivarnow eq -999,cbad,complement=goodpos)
         if cbad lt n_elements(colornow) then begin
            if cbad gt 0 then remove,badpos,colornow,color_ivarnow,filternow
            znow = spsstr[goodobj[i]].zfit
            if cbad le 5 then begin
               kcorrect,colornow,color_ivarnow,znow,kcorrect1,chi2=chi2,filterlist=filternow,absmag=absmag,/magnitude
            endif
            absmagarr[goodpos,i] = absmag
         endif
      endfor
      mask = spsstr(goodobj).mask
      objname = spsstr(goodobj).objname
      slit = spsstr(goodobj).slit
      bands='B,V,R,I,J,K,f606w,f814w'
      save,mask,objname,slit,absmagarr,bands,filename='cl0024_absmag_kcorrected.sav'
   endif
endif else restore,absmagfile
   
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
   plot,lambda(widehbreg),contdiv(widehbreg),psym=10
   hba = a
   if a[0] gt 0.05 and abs(a[1]-4862) lt 2. then begin
      emission = A[0]*EXP(-0.5*((lambda(hbreg)-A[1])/A[2])^2)+A[3]
      absorption = A[4]*EXP(-0.5*((lambda(hbreg)-A[1])/A[5])^2)+1.-A[3]
      oplot,lambda(hbreg),emission,color=fsc_color('green')            
      oplot,lambda(hbreg),absorption,color=fsc_color('green')            
      oplot,lambda(hbreg),hbfit,color=fsc_color('orange')
      EWhb(i) = sqrt(2)*a[0]*a[2]*sqrt(!pi)
      xyouts,0.2,0.7,'EW = '+strtrim(string(EWhb(i),format='(F4.1)'),2),/normal
   endif
   
   OIIIreg = where(lambda gt 5000 and lambda lt 5016)
   wideOIIIreg = [-1.*reverse(indgen(5)+1)+min(OIIIreg),OIIIreg,indgen(5)+1+max(OIIIreg)]
   a = [0.5,5009,2.,0.05,-0.3,3]
   if max(contdiv(oiiireg)) gt 2. then a[0]=max(contdiv(oiiireg))-1
   fita = [1,1,1,1,1,1]
   OIIIfit = curvefit(lambda(OIIIreg),contdiv(OIIIreg),contdivivar(OIIIreg),a,fita=fita,function_name='gauss_emission',itmax=1e5,status=status)
   plot,lambda(wideOIIIreg),contdiv(wideOIIIreg),psym=10
   OIIIa = a
   if a[0] gt 0.05 and abs(a[1]-5009) lt 2. then begin
      emission = A[0]*EXP(-0.5*((lambda(OIIIreg)-A[1])/A[2])^2)+A[3]
      absorption = A[4]*EXP(-0.5*((lambda(OIIIreg)-A[1])/A[5])^2)+1.-A[3]
      oplot,lambda(OIIIreg),emission,color=fsc_color('green')            
      oplot,lambda(OIIIreg),absorption,color=fsc_color('green')            
      oplot,lambda(OIIIreg),oIIIfit,color=fsc_color('orange')
      EWoiii(i) = sqrt(2)*a[0]*a[2]*sqrt(!pi)
      xyouts,0.2,0.3,'EW = '+strtrim(string(EWoiii(i),format='(F4.1)'),2),/normal
   endif
   stop
endfor

;;PLOT
set_plot,'ps'

!p.multi=[0,2,2]
psname=outdir+'cl0024_redsequence.eps'
device, filename = psname,xsize = 12,ysize = 8, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;;1) Red sequence
plot,absmagarr[2,*],absmagarr[0,*]-absmagarr[2,*],yrange=[0,2],psym=1,xtitle='R',ytitle='B-R',title='Cl0024';B-R vs R

;;2) UVJ diagram
;;2) Hb EW
device,/close

end
