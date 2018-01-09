
pro writeegsselection,zcat

; construct list of masks and read into one large mask file
  mask = fix(strmid(zcat.maskname, 0, 4)) 
  mask = mask[sort(mask)] ;sort
  mask = mask[uniq(mask)] ;just this list of uniq. masks
  for imask=0, n_elements(mask)-1 do begin ;loop over observed masks
     filen = string(mask[imask], format='(i4)')
     mask_in = mrdfits(getenv('IDLDEEPPHOTOMETRY_DIR')+'/targetselection/masks'+ $
       '/1HSmask.' +filen+ '.fits', 2,  /silent) 
     if imask eq 0 then masks_keep = mask_in else $  ;start list
                     masks_keep= [masks_keep, mask_in]
  endfor
  masks = masks_keep


; for each object in the zcat get the predicted f_m, magi, magr and 
;   whether it passes the photocut
  f_m = fltarr(n_elements(zcat)) 
  magi = fltarr(n_elements(zcat))
  magr = fltarr(n_elements(zcat))
  highz = fltarr(n_elements(zcat))
  sizes=highz

; note values used in maskmaking don't match today's

magb=magi
pgal=magi
  mask_index = photocut(masks, /index) 
  for i = 0, n_elements(zcat)-1 do begin
     a = where(masks.objno eq long(zcat[i].objname), ct)
     if ct gt 0 then begin
        f_m[i] = masks[a].f_m 
        magi[i] = zcat[i].magi
        magb[i]=zcat[i].magb
        magr[i] = zcat[i].magr
        sizes[i] = masks[a].rg*2.30*.207
        pgal[i]=masks[a].pgal
        tmp = where(mask_index eq a[0], ct2)
        if ct2 gt 0 then highz[i] = 1
     endif   
 endfor

highz= (highz + (magr lt 21.5))<1



; better code
probgt7=photocut(masks,prob=0.05,nom=0.05)


pgal_as_used=( (pgal)>0. ) < 1.

weight=f_m*pgal_as_used

f_m=0.333979  +  0.426868 * weight

zgood = zcat

completeness=(highz*0.9+(1-highz)*0.825)/(1+exp( 1.4*(magr-24.4)))
;f_m=f_m*completeness 

openw,2,'egsselection.txt'

z=translate_objno(zcat)

nz=n_elements(z)
for i=0,nz-1 do begin 
    printf,2,z[i].objno,f_m[i],completeness[i],format='(I15,2F15.5)'
endfor
close,2

return
end
