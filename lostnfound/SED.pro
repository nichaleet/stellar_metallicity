pro SED
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
             'clj1226_IR.cat', $
            'abell2744.cat']
id = [1424,1131,432,1120,1758,552,2137,324,542,99,928,2126]


for i=11,n_elements(file_list)-1 do begin
   lambda = [225,275,336,390,435,475,555,606,625,775,814,850,1050,1100,1250,1400,1600]
   file = path+file_list(i)
   name = strmid(file_list(i),0,strpos(file_list(i),'_IR.cat'))
   rows= file_lines(file)
   data=dblarr(73,rows)
   openr,lun,file,/get_lun
   Readf,lun,data
   free_lun,lun
   
   row_obj = where(data(0,*) eq id(i))
   stop
   mag_obj = [data(13,row_obj),data(19,row_obj),data(25,row_obj),data(31,row_obj),data(37,row_obj),data(43,row_obj),data(49,row_obj),data(55,row_obj),data(61,row_obj),data(67,row_obj),data(73,row_obj),data(79,row_obj),data(85,row_obj),data(91,row_obj),data(97,row_obj),data(103,row_obj),data(109,row_obj)]
   mag_obj_err = [data(14,row_obj),data(20,row_obj),data(26,row_obj),data(32,row_obj),data(38,row_obj),data(44,row_obj),data(50,row_obj),data(56,row_obj),data(62,row_obj),data(68,row_obj),data(74,row_obj),data(80,row_obj),data(86,row_obj),data(92,row_obj),data(98,row_obj),data(104,row_obj),data(110,row_obj)]
   nodata = where(abs(mag_obj) eq 99.)
   if nodata(0) ne -1 then remove,nodata,mag_obj,mag_obj_err,lambda
   ;plot=errorplot(lambda,mag_obj,mag_obj_err,$
   ;               xtitle='$\lambda$(nm)',ytitle='magnitude',$
   ;               title=name+' id:'+string(id(i)),$
   ;              symbol = 'X',linestyle = 'none',sym_thick=3)
; v                                                    stop
stop
endfor
end
