pro SEDv2_abell2744
path='/scr2/nichal/RichardInterestingObjs/'
file_list = ['abell2744_IR.cat']

catid = replicate({file:0,id:0},7)
catid.file = [0,0,0,0,0,0,0]
catid.id  = [1701,1681,2126,348,268,124,58]

count = 0
for i=0,n_elements(file_list)-1 do begin
   numobject = n_elements(where(catid.file eq i))
   
  
   file = path+file_list(i)
   name = strmid(file_list(i),0,strpos(file_list(i),'_IR.cat'))
   rows= file_lines(file)
   data=dblarr(74,rows)

   openr,lun,file,/get_lun
   Readf,lun,data
   free_lun,lun
  ; if name eq 'abell2744' then begin
  ;    ra_pos=[3.5828812,3.5975483,3.5862243,3.5719463,3.5850454,3.5797133,3.5765449]
  ;    dec_pos =[-30.402875,-30.403941,-30.410009,-30.383024,-30.381825,-30.378425,-30.378758]
  ;    find_id_from_pos,ra_pos,dec_pos,data,raind=1,decind=2,idind=0,id,distance
  ;    print, id
  ;    print, distance
  ;    stop
  ; endif

   for j=0,numobject-1 do begin
      lambda = [435,606,814,1050,1250,1400,1600]
      id = catid(count).id
      row_obj = where(data(0,*) eq id)

      mag_obj = [data(13,row_obj),data(20,row_obj),data(27,row_obj),data(34,row_obj),data(41,row_obj),data(48,row_obj),data(55,row_obj)]
      mag_obj_err = [data(14,row_obj),data(21,row_obj),data(28,row_obj),data(35,row_obj),data(42,row_obj),data(49,row_obj),data(56,row_obj)]
      nodata = where(abs(mag_obj) eq 99.)
      if nodata(0) ne -1 then remove,nodata,mag_obj,mag_obj_err,lambda
      ;plot=errorplot(lambda,mag_obj,mag_obj_err,$
      ;               xtitle='$\lambda$(nm)',ytitle='magnitude',$
      ;               title=name+' id:'+string(id),$
      ;               symbol = 'X',linestyle = 'none',sym_thick=3)
      
      
      z_bpz = [data(64,row_obj),data(65,row_obj),data(66,row_obj)]
      ra    = data(1,row_obj)
      dec   = data(2,row_obj)
      mag   = data(27,row_obj)

      print, name, 'id:'+string(id)
      print, 'z:',z_bpz
      print, 'ra, dec:',ra,dec
 
      for lambi=0,n_elements(lambda)-1 do begin
         print, lambda(lambi), mag_obj(lambi)
      endfor
      count=count+1
   endfor
; v                                                    stop
endfor
stop
end
