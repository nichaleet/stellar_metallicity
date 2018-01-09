pro combine_spsstr,mask,combinename,zrange=zrange
;propose: 
;        1) combine and stack spectra if observed in multiple masks
;        2) select only those in zrange if zrange specified
;example:
;combine_spsstr,['rse1','rse2','rse3','rse4','rse5','rse6','rse7','rse8','rse9','rse10','rse11','rse12','rse14','0024_1B','0024_2B','0024_3B','0024_4'],'all_cl0024',zrange=[0.35,0.45]
;combine_spsstr,['ms1','ms2','ms3','ms4','ms5','ms6','ms7','ms8','ms9','ms10','ms11','0451_1','0451_2','0451_3','0451_4'],'all_ms0451',zrange=[0.4,0.65]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nmasks  = n_elements(mask)
allstr = []
directory = '/scr2/nichal/workspace2/sps_fit/data/'
if ~keyword_Set(zrange) then zrange=[0,2.]
for i=0, nmasks-1 do begin
   sciencefits = mrdfits(directory+mask[i]+'/sps_fit.fits.gz',1)
   good = where(sciencefits.zspec gt min(zrange) and sciencefits.zspec lt max(zrange),cgood)
   if cgood gt 0 then allstr = [allstr,sciencefits(good)]
endfor
nobjs = n_elements(allstr)

;find if there is any multiple visits
checkmark = bytarr(nobjs)
megastr = []
for i=0,nobjs-1 do begin
   if checkmark(i) eq 0 then begin
      pos = where(strmatch(allstr.objname,strtrim(allstr(i).objname,2)+'*'),cdup)
      if cdup le 1 then begin ;register spectra
         megastr = [megastr,allstr(i)]
         checkmark(i)=1
      endif else begin ;combine the spectra
         objnow = allstr(pos)
         ;check if they are actually the same object
         if total(objnow.b-objnow(0).b) ne 0 or $
            total(objnow.v-objnow(0).v) ne 0 or $
            total(objnow.r-objnow(0).r) ne 0 or $
            total(abs(objnow.zspec-objnow(0).zspec)) gt 0.005 then $
               stop, 'Diff objects might have the same name'
         locmax = where(objnow.sn eq max(objnow.sn),complement=otherloc)
         bestobj = objnow(locmax)
         otherobj = objnow(otherloc)
         ;interpolate lambda of objects with less SN to the best object
         for ii=0,n_elements(otherobj)-1 do begin
            ;continuum normalized spectra
            otherobj(ii).contdiv = interpol(otherobj(ii).contdiv,otherobj(ii).lambda,bestobj.lambda,/nan)
            otherobj(ii).contdivivar = interpol(otherobj(ii).contdivivar,otherobj(ii).lambda,bestobj.lambda,/nan)
         endfor
         objnow(otherloc) = otherobj
         badpix = where(~finite(objnow.contdivivar),cbadpix)
         if cbadpix gt 0 then objnow.contdivivar(badpix) = 0.

         ;weighted average
         avg_contdiv = total(objnow.contdiv*objnow.contdivivar,2)/total(objnow.contdivivar,2)
         avg_contdivivar = total(objnow.contdivivar,2)

         bestobj.contdiv = avg_contdiv
         bestobj.contdivivar = avg_contdivivar
         bestobj.telldiv = avg_contdiv*bestobj.continuum
         bestobj.telldivivar = avg_contdivivar/(bestobj.continuum)^2
         bestobj.mask = strjoin(objnow.mask,';')
         
         ;calculate new signal to noise
         contmask = bestobj.contmask
         n = n_elements(bestobj.lambda)
         wcont = where(contmask[3:n-4] eq 1)+3
         wcont = wcont[where(finite(bestobj.telldiv[wcont]) and finite(bestobj.continuum[wcont]) and bestobj.continuum[wcont] ne 0)]
         dev = abs((bestobj.telldiv[wcont] - bestobj.continuum[wcont]) / bestobj.continuum[wcont])
         avgdev = mean(dev)
         w = where(dev lt 3.0*avgdev, c)
         if c gt 0 then bestobj.sn = 1.0/mean(dev[w])
         print, ';;;;;;;;;;;;;;;;;;;'
         print, bestobj.objname
         print,'old SN',allstr(pos).sn
         print, 'new SN', bestobj.sn

         ;keep the new structure
         megastr = [megastr,bestobj]
         checkmark(pos) = 1
      endelse
   endif
endfor
;finish checking for multiple visits
file_mkdir,directory+combinename
mwrfits,megastr,directory+combinename+'/sps_fit.fits.gz',/create,/silent
stop
end
