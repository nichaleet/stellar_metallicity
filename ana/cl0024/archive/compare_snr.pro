pro compare_snr
;compare Moran reduced data with Evan reduced data
moran_dir = '/scr2/nichal/workspace2/sps_fit/data/original_results/'
evan_dir  = '/scr2/nichal/workspace2/sps_fit/data/'
masks     = file_search(evan_dir+'00*',/test_directory,/mark_directory)
nmasks    = n_elements(masks)
goodmasks = bytarr(nmasks)
namearr   = strarr(nmasks)
masknum   = findgen(nmasks)
colors = ['red','tomato','dark green','blue', 'royal blue']
esnrarr = []
msnrarr = []
numarr  = []
for nm=0,nmasks-1 do begin
   maskname = strsplit(masks[nm],'/',/extract)
   maskname = maskname[n_elements(maskname)-1]
   namearr[nm] = maskname
   efile = file_search(masks[nm]+'sps_fit.fits.gz',count=ce)
   mfile = file_search(moran_dir+maskname+'/sps_fit.fits.gz',count=cm)
   esnr = []
   msnr = []
   num = []
   if ce ge 1 and cm ge 1 then begin
      escience = mrdfits(efile[0],1)
      mscience = mrdfits(mfile[0],1)
      for i=0, n_elements(escience)-1 do begin
         objname = escience[i].objname
         imatch = where(mscience.objname eq objname,cmatch)
         if cmatch ne 0 then begin
            esnr = [esnr,escience[i].sn]
            msnr = [msnr,mscience[imatch].sn]
            num = [num,nm]
         endif
      endfor
      if nm eq 0 then begin
         set_plot,'ps'
         !p.multi = [0,2,1]
         !p.font = 0
         psname='SNR_comparison.eps'
         device, filename = psname,xsize = 20,ysize = 8, $
                 xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
         plot,esnr,msnr,xtitle='Evan Reduction SNR',ytitle='Moran Reduction SNR',/nodata,xrange=[0,20],yrange=[0,20]
      endif else oplot,esnr,msnr,color=fsc_color(colors[nm]),psym=1
      goodmasks[nm] = 1
      esnrarr = [esnrarr,esnr]
      msnrarr = [msnrarr,msnr]
      numarr = [numarr,num]
   endif
endfor
oplot,!x.crange,!y.crange,linestyle=2
goodm = where(goodmasks eq 1,cgoodmasks)
al_Legend,namearr(goodm),psym=1,color=colors(goodm),box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0

for i=0,cgoodmasks-1 do begin
   curnum = masknum(goodm(i))
   goodobj = where(numarr eq curnum)
   if i eq 0 then plothist,esnrarr(goodobj),xtitle='SNR',/nodata
   plothist,esnrarr(goodobj),color=fsc_color(colors(curnum)),/overplot
endfor

device,/close
stop
end
