pro get_mass,mask,redo=redo

directory = '/scr2/nichal/workspace2/sps_fit/data/'+mask+'/'
sciencefits = directory+'sps_fit.fits.gz'
scienceall = mrdfits(sciencefits,1,/silent)

;goodobj = where(scienceall.good eq 1 or scienceall.sn gt 5. or scienceall.goodfit eq 1 or scienceall.feh ne -999., cgals)
goodobj = where(scienceall.z ge 0.38 and scienceall.z le 0.42, cgals)

if cgals gt 0 then begin
   color = dblarr(8,cgals);b,v,r,i,j,k,f606,f814
   ;assume magnitudes were in vega system
   color[0,*] = scienceall[goodobj].b-0.09
   color[1,*] = scienceall[goodobj].v+0.02
   color[2,*] = scienceall[goodobj].r+0.21
   color[3,*] = scienceall[goodobj].i+0.45
   color[4,*] = scienceall[goodobj].j+0.91
   color[5,*] = scienceall[goodobj].k+1.85
   color[6,*] = scienceall[goodobj].f606w
   color[7,*] = scienceall[goodobj].f814w

   color_err = dblarr(8,cgals)      ;b,v,r,i,j,k,f606,f814
   color_err[0,*] = scienceall[goodobj].berr
   color_err[1,*] = scienceall[goodobj].verr
   color_err[2,*] = scienceall[goodobj].rerr
   color_err[3,*] = scienceall[goodobj].ierr
   color_err[4,*] = scienceall[goodobj].jerr
   color_err[5,*] = scienceall[goodobj].kerr
   color_err[6,*] = scienceall[goodobj].f606werr
   color_err[7,*] = scienceall[goodobj].f814werr
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
      badpos = where(colornow gt 50. or colornow lt 5 or color_ivarnow eq -999,cbad)
      if cbad lt n_elements(colornow) then begin
         if cbad gt 0 then remove,badpos,colornow,color_ivarnow,filternow
         znow = scienceall[goodobj[i]].zfit
         if znow eq 0 then znow = scienceall[goodobj[i]].zspec
         massnow = scienceall[goodobj[i]].logmstar
         if massnow eq 0 or keyword_set(redo) then begin
            if cbad le 5 then begin
               kcorrect,colornow,color_ivarnow,znow,kcorrect1,chi2=chi2,filterlist=filternow,/magnitude,mass=mass,b300=b300
               mass = 0.52*mass ;fix for the hubble constant h=0.72
               mass = alog10(mass)
               if ~finite(mass) then mass = -999d 
               if ~finite(b300) then b300 = -999d
               massarr(i) = mass
               b300arr(i) = b300
            endif
         endif else massarr(i) = massnow
      endif
   endfor
   stop
   scienceall[goodobj].logmstar = massarr
   if ~tag_exist(scienceall,'b300') then begin
      newscienceall = []
      for i=0,n_Elements(scienceall)-1 do begin
         science = scienceall[i]
         newscience = create_struct(science,'b300',-999d)
         newscienceall = [newscienceall,newscience]
      endfor
   endif 
   scienceall[goodobj].b300 = b300arr
   ;mwrfits,scienceall,sciencefits,/create,/silent
endif   

end
