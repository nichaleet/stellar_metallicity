;+
; NAME:
;   QUICK_CALIBS
; PURPOSE:
;   check for calibration files and start quicklook reductions
; CALLING SEQUENCE:
;   quick_calibs
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;   run in data directory
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;
;----------------------------------------------------------------------

function deimos_hdr_struct, N

  str = {fname:   '', $
         mask:    '', $
         exptime: 0., $
         obstype: '', $
         lamps:   '', $
         grating: '', $
         wave: 0., $
         object:  ''}

  if n_elements(N) NE 0 AND N ne 0 then str = replicate(str, N)
  return, str
end

pro quick_plan, s, maskname, rawdatadir
  w = where(s.mask eq maskname)
  a = s[w]
  flatind    = where(strpos(a.lamps,'Qz') ne -1, nflats)
  arcind     = where( ( (strpos(a.lamps, 'Ar') ne -1) + $
                      (strpos(a.lamps, 'Kr') ne -1) + $
                      (strpos(a.lamps, 'Xe') ne -1) + $
;                      (strpos(a.lamps, 'Cd') ne -1) or $
;                      (strpos(a.lamps, 'Hg') ne -1) or $
                      (strpos(a.lamps, 'Ne') ne -1) ge 3), narcs)
  scienceind = where(strpos(a.lamps, 'Off') ne -1 AND (a.exptime GT 250) AND (a.wave gt 6500), nscience)



;  if (nflats GE 1) AND (narcs GE 1) AND (nscience GE 1) then begin 
  if (nflats GE 1) AND (narcs GE 1)  then begin 

     planname = maskname+'.plan'
     openw, wlun, planname, /get_lun
     
     printf, wlun, '# Plan file auto-generated by deimos_planfile.pro ', $
       systime()
     printf, wlun, '# Grating: ', a[0].grating, '       Wavelength: ', $
       a[arcind[0]].wave
     printf, wlun, maskname, format='("MASK: ",A)'
     printf, wlun, rawdatadir, format='("RAWDATADIR: ",A)'
     printf, wlun, 'polyflag   - use polyflag for fitting lambda'
     printf, wlun, a[flatind].fname, format='("FLATNAME: ",A)'
;     printf, wlun, a[arcind[narcs/2.]].fname, format='("ARCNAME: ",A)'
     printf, wlun, a[arcind[n_elements(arcind)-1]].fname, format='("ARCNAME: ",A)'

;     printf, wlun, a[scienceind].fname, format='("SCIENCENAME: ",A)'
     printf, wlun, 'SIMPLEST: 15'
;     printf, wlun, 'CHIPS: 2,6'
     free_lun, wlun
     print, planname
;     spawn, 'more '+planname

  endif else begin 
     
     print, 'Cannot create plan file for mask ', maskname

  endelse 

  return
end



pro quick_calibs
  spawn, 'obslog', logout
  nlog = n_elements(logout)
  nlogfiles = n_elements(logout)-10
  fnames = findfile('d*.fits')
  fileinfo = deimos_hdr_struct(nlogfiles)

print,'obslog done'

  cd,current=cwd
  deimos_data = getenv('DEIMOS_DATA')
  if strpos(cwd, deimos_data) ge 0 then $
    rawdatadir = strmid(cwd, strlen(deimos_data), $
                        strlen(cwd)-strlen(deimos_data)) $
        else rawdatadir = cwd

  for i=0, nlogfiles-1 do begin 
        fileinfo[i].fname = fnames[i] 
        fileinfo[i].lamps=strmid(logout[i+10], 29, 10) 
;        fileinfo[i].mask = strcompress(strmid(logout[i+10], 40, 7), /REMOVE) 
        if strlen(fnames[i]) gt 0 then begin
		header=headfits(fnames[i],/sil,err=err)		
	 if strlen(err) eq 0 then $
		fileinfo[i].mask = strcompress(sxpar(header,'SLMSKNAM'),/REMOVE)
	endif
        fileinfo[i].grating = strmid(logout[i+10], 49, 6) 
        fileinfo[i].wave = float(strmid(logout[i+10], 57, 4)) 
        fileinfo[i].exptime = float(strmid(logout[i+10], 77, 4)) 
  endfor   


  masks = fileinfo[sort(fileinfo.mask)].mask
	masks=masks[where(strlen(masks) gt 0)]



  masks=strcompress(masks,/remove)
  masks = masks[uniq(masks)]
  whok = where(strlen(masks) ge 1 AND strpos(masks, 'None') lt 0, nmasks)




  if nmasks gt 0 then begin
     masks = masks[whok]


     for i=0, nmasks-1 do begin
;  for i=0, 2 do begin
       maskprocessed = (n_elements(findfile(masks[i]+'/arc*')) ge 4 $
        AND n_elements(findfile(masks[i]+'/cal*')) ge 10) $
        OR strpos(findfile(masks[i]+'/buckledmask.txt'),'buckle') ge 0 $
        OR strpos(masks[i],  'GOS') eq 0 $
        OR strpos(masks[i],  'GOH') eq 0 $
        OR strpos(masks[i],  'Long') eq 0 $
        OR strpos(masks[i],  'Lng') eq 0 

       if maskprocessed eq 0  then begin

 
          nflats = total(fileinfo.mask eq masks[i] AND $
                       strpos(fileinfo.lamps, 'Qz') ge 0)
          narcs = total(fileinfo.mask eq masks[i] AND $
                       strpos(fileinfo.lamps, 'Ne') ge 0)

          if nflats ge 2 and narcs gt 0 then begin
              print, 'processing mask ', masks[i]
              spawn, 'mkdir '+masks[i]
              quick_plan, fileinfo, masks[i], rawdatadir

              spawn, '/usr/bin/cp '+masks[i]+'*.plan '+masks[i]
              cd, masks[i]
           
              make_bintab_file, masks[i]+'*.plan', /quick

; write calibSlit files
              spawn, '/usr/bin/rm calib.log'
              string = 'set clobber ; echo "deimos_mask_calibrate, /noplot, /quick" | idl >& calib.log'
              spawn, string
               quick_arcqa
             cd,cwd
          endif

       endif    
  endfor

  endif

  recentfileinfo = fileinfo[(nlogfiles-10) > 0:*]
  doingscience = total(recentfileinfo.exptime gt 1000. AND $
                       strpos(recentfileinfo.lamps, 'Off') ge 0) gt 3 $
             AND n_elements(findfile('*/arc*')) gt (nmasks-1)*2 
    
  if doingscience eq 0 then quick_calibs

return
end





