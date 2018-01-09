function thomas_interp,res
;Trilinear interpolate Thomas11b table at resolution of res(dex) usually res= 0.025 dex
;retern grids of indices prediction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Read the Thomas11b table
thomasfile = '/scr2/nichal/Thomass11b/tmj_padova.dat'
readcol,thomasfile,age,ZH,AlphaFe,HdA,HdF,CN1,CN2,Ca4227,G4300,HgA,HgF,Fe4383,Ca4455,Fe4531,C24668,Hb,Fe5015,Mg1,Mg2,Mgb,Fe5270,Fe5335,Fe5406,Fe5709,Fe5782,NaD,TiO1,TiO2,format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D'
thomas_orig = {age:age,ZH:ZH,AlphaFe:alphafe,HdA:HdA,HdF:HdF,CN1:CN1,CN2:CN2,Ca4227:Ca4227,G4300:G4300,HgA:HgA,HgF:HgF,Fe4383:Fe4383,Ca4455:Ca4455,Fe4531:Fe4531,C24668:C24668,Hb:Hb,Fe5015:Fe5015,Mg1:Mg1,Mg2:Mg2,Mgb:Mgb,Fe5270:Fe5270,Fe5335:Fe5335,Fe5406:Fe5406,Fe5709:Fe5709,Fe5782:Fe5782,NaD:NaD,TiO1:TiO1,TiO2:TiO2}
thomas_orig.age = alog10(thomas_orig.age)+9.
;Find how many output combinations there are
thomasage = age(uniq(age,sort(age))) ;gyr
thomasage = alog10(thomasage)+9.     ;log(age)
thomasz   = ZH(uniq(ZH,sort(ZH)))
thomasalpha = alphafe(uniq(alphafe,sort(alphafe)))
maxage = max(thomasage)
maxz   = max(thomasz)
maxalpha = max(thomasalpha)

noutage = ceil((maxage-min(thomasage))/res)
noutz   = ceil((maxz-min(thomasz))/res)
noutalpha = ceil((maxalpha-min(thomasalpha))/res)

nout = noutage*noutz*noutalpha

;make output structure
outstr = {age:1.d,ZH:1.d,AlphaFe:1.d,HdA:1.d,HdF:1.d,CN1:1.d,CN2:1.d,Ca4227:1.d,G4300:1.d,HgA:1.d,HgF:1.d,Fe4383:1.d,Ca4455:1.d,Fe4531:1.d,C24668:1.d,Hb:1.d,Fe5015:1.d,Mg1:1.d,Mg2:1.d,Mgb:1.d,Fe5270:1.d,Fe5335:1.d,Fe5406:1.d,Fe5709:1.d,Fe5782:1.d,NaD:1.d,TiO1:1.d,TiO2:1.d}
thomas_fine = replicate(outstr,noutage,noutz,noutalpha)

;Interpolate into finer grids at uniform interval res dex
;Trilinear_interpolation
ind = 0L
for indage=0,noutage-1 do begin
   for indz=0,noutz-1 do begin
      for indalpha = 0,noutalpha-1 do begin
         age = thomasage[0]+float(indage)*res<maxage
         z = thomasz[0]+float(indz)*res<maxz
         alpha = thomasalpha[0]+float(indalpha)*res<maxalpha
         thomas_fine[indage,indz,indalpha].age = age
         thomas_fine[indage,indz,indalpha].zh  = z
         thomas_fine[indage,indz,indalpha].alphafe = alpha

         wage = value_locate(thomasage,age)>0
         if wage eq noutage-1 then wage -= 1
         iage = [wage,wage+1]
         dage = (age-thomasage[iage[0]])/(thomasage[iage[1]]-thomasage[iage[0]])

         wz = value_locate(thomasz,z)>0
         if wz eq noutz-1 then wz -= 1
         iz = [wz,wz+1]
         dz = (z-thomasz[iz[0]])/(thomasz[iz[1]]-thomasz[iz[0]])

         walpha = value_locate(thomasalpha,alpha)>0
         if walpha eq noutalpha-1 then walpha -= 1
         ialpha = [walpha,walpha+1]
         dalpha = (alpha-thomasalpha[ialpha[0]])/(thomasalpha[ialpha[1]]-thomasalpha[ialpha[0]])

         ;Use Wiki page for trilinear interpolation         
         ;first find the position of the 8 corners around the point
         wijk = intarr(2,2,2)
         for i=0,1 do begin
            for j=0,1 do begin
               for k=0,1 do begin
                  w = where(thomas_orig.age eq thomasage[iage[i]] and thomas_orig.zh eq thomasz[iz[j]] and thomas_orig.alphafe eq thomasalpha[ialpha[k]],c)
                  if c eq 0 then message, 'Thomas table not found'                  
                  wijk[i,j,k] = w
               endfor
            endfor
         endfor

         ;loop over each lick index
         lickind = dblarr(28)
         for il=3,27 do begin
            ;interpolate along age
            c00 = (thomas_orig.(il)[wijk[0,0,0]]*(1.-dage))+(thomas_orig.(il)[wijk[1,0,0]]*dage)
            c01 = (thomas_orig.(il)[wijk[0,0,1]]*(1.-dage))+(thomas_orig.(il)[wijk[1,0,1]]*dage)
            c10 = (thomas_orig.(il)[wijk[0,1,0]]*(1.-dage))+(thomas_orig.(il)[wijk[1,1,0]]*dage)
            c11 = (thomas_orig.(il)[wijk[0,1,1]]*(1.-dage))+(thomas_orig.(il)[wijk[1,1,1]]*dage)
            c0 = c00*(1.-dz)+c10*dz
            c1 = c01*(1.-dz)+c11*dz
            c = c0*(1.-dalpha)+c1*dalpha
            thomas_fine[indage,indz,indalpha].(il)=c
         endfor
         ind += 1
         if ind mod 10000 eq 0 then print, 'doing '+strtrim(string(ind),2)+'/'+strtrim(string(nout),2)
      endfor
   endfor
endfor

return, thomas_fine
end

