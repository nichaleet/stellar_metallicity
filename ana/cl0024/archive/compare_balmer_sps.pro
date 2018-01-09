pro compare_balmer_sps

if ~file_test('merge_moran_sps2.sav') then begin

   directory = '/scr2/nichal/workspace2/sps_fit/data/original_results/'
   masks=['0024-1B','0024-2B','0024-3B','0024-4','0024-6']
   cat = mrdfits('/scr2/nichal/workspace2/sps_fit/data/moran05_table.fit',1)
   cat = rename_tags(cat,'z','z_Moran')
   cat = rename_tags(cat,'R','R_Moran')
   setplot,12
   plot,cat.raj2000,cat.dej2000,psym=4,background=255,color=0,yrange=[16.9,17.35],xtitle='RA',ytitle='Dec',ystyle=1
   write_ds9_regionfile,cat.raj2000,cat.dej2000,filename='imagefile/Moran05obj.reg',symbol='diamond',color='cyan'
   newscience = []
   count=0
   for i=0,n_elements(masks)-1 do begin
      science_All  = mrdfits(directory+masks(i)+'/sps_fit.fits.gz',1)
      oplot,science_all.ra,science_all.dec,psym=2,color=100+20*i
      write_ds9_regionfile,science_all.ra,science_all.dec,filename='imagefile/'+masks(i)+'.reg',symbol='diamond',color='red'

      goodfit      = where(science_all.sn gt 1.)
      good_science = science_All(goodfit)
      ;oplot,good_science.ra,good_science.dec,psym=2,color=100+20*i
      for j=0,n_elements(goodfit)-1 do begin
         gcirc,2,good_science(j).ra,good_science(j).dec,cat.raj2000,cat.dej2000,dis
                                ;print, min(dis)
         if min(dis) lt 1. then begin
            count +=1
            pos = where(dis eq min(dis))
            newscience=[newscience,create_Struct(good_science(j),cat(pos))]
         endif
      endfor
   endfor
   save,newscience,filename='merge_moran_sps.sav'
   oplot,newscience.ra,newscience.dec,psym=6,color=20,symsize=2
   write_ds9_regionfile,newscience.ra,newscience.dec,filename='imagefile/matchedobjects.reg',symbol='circle',color='green'

   al_Legend,['Moran2005','Objects in 5 masks','matched found'],psym=[4,2,6],color=[0,120,40],position=[6.89,16.95],/right,box=0,font=0,thick=2
stop,count
endif else restore,'merge_moran_sps.sav'

;Find age/metallicity from Balmer lines. Interpolating from worthey97 table 6
restore,'/scr2/nichal/workspace2/wortheytable.sav'
worthey = worthey_table  
objs_num = n_elements(newscience)
agearr = worthey(uniq(worthey.age,sort(worthey.age))).age
low_age = where(agearr lt 8.)
high_age = where(agearr ge 8.)
feharr = worthey(uniq(worthey.feh,sort(worthey.feh))).feh
!p.multi=[0,3,2]
for i=0, objs_num-1 do begin
   HDA = newscience(i).HDA
   HGA = newscience(i).HGA
   chisq = (HDA-worthey.HDA)^2+(HGA-worthey.HGA)^2
   minloc = where(chisq eq min(chisq))
   first_age = worthey(minloc).age
   first_feh = worthey(minloc).feh
   first_chisq = chisq(minloc)
   print,first_age
   ;marginalize
   ;marginalize over feh
   prob = max(chisq)-chisq
   pdf_age = []
   for ii=0,n_Elements(agearr)-1 do begin
      loc     = where(worthey.age eq agearr(ii))
      feh_now = worthey(loc).feh
      prob_now= prob(loc) 
      pdf_age = [pdf_age,int_tabulated(feh_now,prob_now)]
   endfor
   plot,agearr,pdf_age,psym=10   
   minloc = where(dis_Sq eq min(dis_Sq))
   first_age = worthey(minloc).age
   first_feh = worthey(minloc).feh
   first_dis = sqrt(dis_sq(minloc))
   plot,agearr(low_age),pdf_age(low_age),psym=10
   plot,agearr(high_age),pdf_age(high_age),psym=10

   ;marginalize over age
   pdf_feh = []
   for ii=0,n_Elements(agearr)-1 do begin
      loc     = where(worthey.age eq agearr(ii))
      feh_now = worthey(loc).feh
      prob_now= prob(loc) 
      pdf_age = [pdf_age,int_tabulated(feh_now,prob_now)]
   endfor
   plot,agearr,pdf_age,psym=10
   plot,agearr(low_age),pdf_age(low_age),psym=10
   plot,agearr(high_age),pdf_age(high_age),psym=10
   stop
endfor
stop
end


 ;  remove,minloc,dis_sq
 ;  minloc = where(dis_sq eq min(dis_sq))
 ;  second_age = worthey(minloc).age
 ;  second_feh = worthey(minloc).feh 
 ;  second_dis = sqrt(dis_sq(minloc))
 ;  remove,minloc,dis_sq
 ;  minloc = where(dis_sq eq min(dis_sq))
 ;  third_age = worthey(minloc).age
 ;  third_feh = worthey(minloc).feh 
 ;  third_dis = sqrt(dis_sq(minloc))
 ;  print, first_age,first_feh,first_dis
 ;  print, second_age,second_feh,second_dis
 ;  print, third_age,third_feh,third_dis
