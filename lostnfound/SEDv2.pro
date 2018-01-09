pro SEDv2
path='/scr2/nichal/RichardInterestingObjs/'
file_list = ['abell_2261_IR.cat',$
             'abell_383_IR.cat',$
             'macs0717_IR.cat',$
             'macs1115_IR.cat',$
             'macs1149_IR.cat',$
             'macs1206_IR.cat',$
             'macs1931_IR.cat',$
             'ms2137_IR.cat',$
             'rxj1532_IR.cat',$
             'rxj2248_IR.cat',$
             'clj1226_IR.cat']

catid = replicate({file:0,id:0},17)
catid.file = [0,0,1,2,3,4,5,5,5,6,6,6,7,7,8,9,10]
catid.id  = [1923,1424,1131,432,1120,1758,429,1858,552,139,2162,2137,324,400,542,99,928]

count = 0
for i=0,n_elements(file_list)-1 do begin
   numobject = n_elements(where(catid.file eq i))
   
  
   file = path+file_list(i)
   name = strmid(file_list(i),0,strpos(file_list(i),'_IR.cat'))
   rows= file_lines(file)
   data=dblarr(125,rows)

   openr,lun,file,/get_lun
   Readf,lun,data
   free_lun,lun
;Find the object that is next to ms2137 id 324
;   if name eq 'ms2137' then begin
;      ra_pos=[325.05375,325.0524,325.0534]
;      dec_pos =[-23.6473,-23.647933,-23.647833]
;      find_id_from_pos,ra_pos,dec_pos,data,raind=1,decind=2,idind=0,id,distance
;      print, id
;      print, distance
;      for aa=0,n_elements(id)-1 do begin
;         print, data(1,id(aa)-1),data(2,id(aa)-1),data(115,id(aa)-1),data(116,id(aa)-1),data(117,id(aa)-1)
;      endfor
;      stop
;   endif

   for j=0,numobject-1 do begin
      lambda = [225,275,336,390,435,475,555,606,625,775,814,850,1050,1100,1250,1400,1600]
      id = catid(count).id
      row_obj = where(data(0,*) eq id)
      print, name, 'id:'+string(id)
      mag_obj = [data(13,row_obj),data(19,row_obj),data(25,row_obj),data(31,row_obj),data(37,row_obj),data(43,row_obj),data(49,row_obj),data(55,row_obj),data(61,row_obj),data(67,row_obj),data(73,row_obj),data(79,row_obj),data(85,row_obj),data(91,row_obj),data(97,row_obj),data(103,row_obj),data(109,row_obj)]
      mag_obj_err = [data(14,row_obj),data(20,row_obj),data(26,row_obj),data(32,row_obj),data(38,row_obj),data(44,row_obj),data(50,row_obj),data(56,row_obj),data(62,row_obj),data(68,row_obj),data(74,row_obj),data(80,row_obj),data(86,row_obj),data(92,row_obj),data(98,row_obj),data(104,row_obj),data(110,row_obj)]
      flux_obj = [data(16,row_obj),data(22,row_obj),data(28,row_obj),data(34,row_obj),data(40,row_obj),data(46,row_obj),data(52,row_obj),data(58,row_obj),data(64,row_obj),data(70,row_obj),data(76,row_obj),data(82,row_obj),data(88,row_obj),data(94,row_obj),data(100,row_obj),data(106,row_obj),data(112,row_obj)]
      flux_obj_err = [data(17,row_obj),data(23,row_obj),data(29,row_obj),data(35,row_obj),data(41,row_obj),data(47,row_obj),data(53,row_obj),data(59,row_obj),data(65,row_obj),data(71,row_obj),data(77,row_obj),data(83,row_obj),data(89,row_obj),data(95,row_obj),data(101,row_obj),data(107,row_obj),data(113,row_obj)]

      nodata = where(abs(mag_obj) eq 99.)
      if nodata(0) ne -1 then remove,nodata,mag_obj,mag_obj_err,lambda,flux_obj,flux_obj_err
      ;if name eq 'ms2137' then plot=errorplot(lambda,flux_obj,flux_obj_err,$
      ;               xtitle='$\lambda$(nm)',ytitle='flux',$
      ;               title=name+' id:'+string(id),$
      ;               symbol = 'X',linestyle = 'none',sym_thick=3)
      
      z_bpz = [data(115,row_obj),data(116,row_obj),data(117,row_obj)]
      ra    = data(1,row_obj)
      dec   = data(2,row_obj)
      mag   = mag_obj(11)
      
      
      print, 'z:',z_bpz
      print, 'ra, dec:',ra,dec
      print, 'mag at 850 nm:',mag

      if id eq 324 or id eq 400 then begin
         for lambi=0,n_elements(lambda)-1 do begin
            print, lambda(lambi), mag_obj(lambi), mag_obj_err(lambi)          
         endfor
       ;  stop
      endif


      count=count+1
   endfor
; v                                                    stop
endfor
stop
end
