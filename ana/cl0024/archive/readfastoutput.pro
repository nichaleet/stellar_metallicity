function readfastoutput,file
  dir = '/scr2/nichal/workspace2/fast/'
  readcol,dir+file,id,redshift,ltau,metal,lage,Av,lmass,lsfr,lssfr,la2t,chi2
  feh = alog10(metal/0.019)
  foutall = []
  for i=0,n_elements(id)-1 do begin
     if metal[i]eq 0.0008 then feherr=0.3
     if metal[i]eq 0.0031 then feherr=0.25
     if metal[i]eq 0.0096 then feherr=0.15
     if metal[i]eq 0.019 then feherr=0.1
     if metal[i]eq 0.03 then feherr=0.1
     fout={id:id[i],z:redshift[i],ltau:ltau[i],feh:feh[i],feherr:feherr,age:10.^(lage[i]-9.),ageerr:0.7,Av:Av[i],lmass:lmass[i],lsfr:lsfr[i],lssfr:lssfr[i],la2t:la2t[i],chi2:chi2[i]}
     foutall = [foutall,fout]
  endfor
  return,foutall
end
