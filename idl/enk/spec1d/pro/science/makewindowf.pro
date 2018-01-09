;+
; NAME:
;   makewindowf
;
; PURPOSE:
;   create a window function of regions observed, based on mask coords.
;
; CATEGORY:
;   correlation analysis
;
; CALLING SEQUENCE:
;   windowf =makewindowf,  nra,ndec, masks=masks, ralim=ralim, $
;               declim=declim, nozcorr=nozcorr
;
; INPUTS:
;   nra -- number of pixels in RA direction (default =3000)
;   ndec -- number of pixels in dec direction (default=1000)
;   masks--     list of masks data structures
;
; KEYWORDS:
;   ralim  -- if set, limits on RA, otherwise gleaned roughly from central mask
;   declim -- if set, limits on dec
;   nozcorr -- if set, do NOT apply the z-completeness correction per mask 
;
; OUTPUTS:
;   windowf -- window function
;
; MODIFICATION HISTORY:
;  md 09dec02
;  alc oct05 - added gaps, dark corners, fixed math for overlap regions for masks designed on different dates
;-


pro mask_ccr,  x, y,  ra,  dec, rac, decc, pa, bow,  height
; rotate from mask coords. into ra,dec (degrees)
  sinpa = sind(pa) &  cospa= cosd(pa)
  x = x -bow*((2.*y/height)^2-.5)  ;maximum bow at y=height/2. and middle 
  ra =  (x*cospa + y*sinpa)/3600./cosd(decc) + rac
  dec = (-x*sinpa + y*cospa)/3600. + decc

  return
end

function makewindowf,  nra, ndec, masks=masks,  ralim=ralim, declim=declim, nozcorr = nozcorr

;   sfraction = 0.65    ; - this was the value in the code before Oct05
   sfraction = 0.6

   if n_params() eq 0 then begin
     nra = 3000
     ndec =  1000
   endif
   if n_params() eq 1 then ndec = 1000 ;default array sizes

   if NOT keyword_set(ralim) then  ralim = $ ;degrees
            [min(masks.ra)-.1, max(masks.ra)+.1]
   if NOT keyword_set(declim) then declim = $ ;degrees
            [min(masks.dec)-.14, max(masks.dec)+.14]
   ra_inc = (ralim[1]-ralim[0])/nra
   dec_inc = (declim[1]-declim[0])/ndec


; read in redshift completeness
   zfile='$DEEP2PRODUCTS/zspec_archive/zcomp.dat'
;   zfile='$DEEP2PRODUCTS/zcomp.dat'
   readcol, zfile, maskname2, zcomp, zcompblue, zcompred, format = 'I,F,F,F'

; HAVE TO MAKE SEPARATE WINDOW FUNCTIONS FOR EACH MASK DESIGN DATE AND COMBINE AT END

; read in the mask design dates
   readcol,'$DEEP2PRODUCTS/zcat_archive/maskdates.txt',maskname,filepath,desdate,MJDdesdate,MJDobsdate,fracinzcat,format="(I8,A50,A15,D10.2,D10.2,F5.2)"  
; only keep masks relevant here:
   match,masks.masknumber,maskname,submasks,submaskname
   maskname=maskname[submaskname]
   MJDdesdate=MJDdesdate[submaskname]

; add 0.5 days to any masks in 'bottom' rows - so they won't have the same date as top row masks
   bottom=where(maskname-100*(maskname/100) ge 40) ; bottom row
   if bottom[0] ne -1 then MJDdesdate[bottom]=MJDdesdate[bottom]+0.5
   filler=where(maskname-100*(maskname/100) ge 80) ; filler masks
   if filler[0] ne -1 then MJDdesdate[filler]=MJDdesdate[filler]+0.25

   desdate=dblarr(n_elements(masks))
   for i=0,n_elements(masks)-1 do desdate[i]=MJDdesdate(where(maskname eq masks[i].masknumber))
   sortdate=desdate(sort(desdate))
   uniqdate=sortdate[uniq(sortdate)]  

   windowf = 1.d0
   
   for j=0,n_elements(uniqdate)-1 do begin ; loop through the unique mask design dates
       tmasks=masks[where(MJDdesdate eq uniqdate[j])]
; designate if these masks are 2002 masks or not 
       date=0
       if (uniqdate[j]+2400000) lt 2452275.5 then date=1

       Y_ang_offset = ((187. + tmasks[0].width +187.)/2 -4.5*60.)/3600.
 ;offset of ra1 from pointing center, in arcsec, converted to degrees
;start first pass on selection
       nmasks = n_elements(tmasks)

       tempwindowf = dblarr(nra, ndec) ;declare window array
       
       for i=0,nmasks-1  do begin ;plot boundaries of masks
           
           rac  = tmasks[i].ra - Y_ang_offset/cosd(tmasks[i].dec)*cosd(tmasks[i].pa) 
           decc = tmasks[i].dec+ Y_ang_offset*sind(tmasks[i].pa)
       
; get redshift completeness for this mask
           thismask = where(maskname2 eq tmasks[i].masknumber, ct)
       
; if mask wasn't observed use z completeness of 100% (this allows us
          ; to make a cut in chi-squared instead of only using galaxies which
          ; have been pspec'd
           if ct gt 0 then zfactor = zcompred[thismask] else zfactor = 1.


; get mask coords core region
;    maskx = [0., .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, 0.0]
;    masky = [-.5, -.5, -.4, -.3, -.2, -.1, 0., .1,  .2,  .3, .4, .5, .5]
;    maskx = [maskx, -maskx]     ;concatenate array
;    masky = [masky, reverse(masky)] ;concatenate again
    
           maskx1=[-0.5,0.5,0.5,-0.5,-0.5]
           masky1=[-0.5,-0.5,-0.259,-0.259,-0.5]
           maskx2=[-0.5,0.5,0.5,-0.5,-0.5]
           masky2=[-0.248,-0.248,-0.001,-0.001,-0.248]
           maskx3=[-0.5,0.5,0.5,-0.5,-0.5]
           masky3=[0.01,0.01,0.257,0.257,0.01]
           maskx4=[-0.5,0.5,0.5,-0.5,-0.5]
           masky4=[0.268,0.268,0.5,0.5,0.268]
           if date eq 1 then begin
               masky1[2:3]=masky1[2:3]-0.005
               masky2=masky2-0.005
               masky3=masky3-0.005
               masky4[0:1]=masky4[0:1]-0.005
               masky4[4]=masky4[4]-0.005
           endif
           maskx1=-1*maskx1
           maskx2=-1*maskx2
           maskx3=-1*maskx3
           maskx4=-1*maskx4

           mask_ccr,  maskx1*tmasks[i].cwidth, masky1*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, tmasks[i].bow, $ ;include bow
             tmasks[i].height 
           xwindow = (ra_mask-ralim[0])/ra_inc ;get mask x,y
           ywindow = (dec_mask-declim[0])/dec_inc
           jcore1 = polyfillv(xwindow, ywindow, nra, ndec)
           
           mask_ccr,  maskx2*tmasks[i].cwidth, masky2*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, tmasks[i].bow, $ ;include bow
             tmasks[i].height 
           xwindow = (ra_mask-ralim[0])/ra_inc ;get mask x,y
           ywindow = (dec_mask-declim[0])/dec_inc
           jcore2 = polyfillv(xwindow, ywindow, nra, ndec)
       
           mask_ccr,  maskx3*tmasks[i].cwidth, masky3*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, tmasks[i].bow, $ ;include bow
             tmasks[i].height 
           xwindow = (ra_mask-ralim[0])/ra_inc ;get mask x,y
           ywindow = (dec_mask-declim[0])/dec_inc
           jcore3 = polyfillv(xwindow, ywindow, nra, ndec)
       
           mask_ccr,  maskx4*tmasks[i].cwidth, masky4*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, tmasks[i].bow, $ ;include bow
             tmasks[i].height 
           xwindow = (ra_mask-ralim[0])/ra_inc ;get mask x,y
           ywindow = (dec_mask-declim[0])/dec_inc
           jcore4 = polyfillv(xwindow, ywindow, nra, ndec)
           
           jcore=[jcore1,jcore2,jcore3,jcore4]


; get mask coords, outer region
;    maskxf = [-0.5,  0.5, 0.5,-0.5]
;    maskyf = [-0.5, -0.5, 0.5, 0.5]
           maskxf1=[-0.5, 0.1,  0.5,   0.5,   -0.5,  -0.5]
           maskyf1=[-0.5, -0.5, -0.404,-0.259,-0.259,-0.5]
           maskxf2=[-0.5, 0.5,  0.5,   0.367, -0.5,  -0.5]
           maskyf2=[-0.248,-0.248,-0.1,-0.001,-0.001,-0.248]
           maskxf3=[-0.5, 0.375,0.5,   0.5,   -0.5,  -0.5]
           maskyf3=[0.01, 0.01, 0.12,  0.257, 0.257, 0.01]
           maskxf4=[-0.5, 0.5,  0.5,   0.1,   -0.5,  -0.5]
           maskyf4=[0.268,0.268,0.404, 0.5,   0.5,   0.268]
           ; if a 2002 mask then the gaps shift by 0.005
           if date eq 1 then begin
               maskyf1[3:4]=maskyf1[3:4]-0.005
               maskyf2[0:1]=maskyf2[0:1]-0.005
               maskyf2[3:5]=maskyf2[3:5]-0.005
               maskyf3[0:1]=maskyf3[0:1]-0.005
               maskyf3[3:5]=maskyf3[3:5]-0.005
               maskyf4[0:1]=maskyf4[0:1]-0.005
               maskyf4[5]=maskyf4[5]-0.005
           endif
           maskxf1=-1*maskxf1
           maskxf2=-1*maskxf2
           maskxf3=-1*maskxf3
           maskxf4=-1*maskxf4


           mask_ccr,  maskxf1*tmasks[i].width, maskyf1*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, 0. , tmasks[i].height ;no bow
           xwindow = (ra_mask-ralim[0])/ra_inc
           ywindow = (dec_mask-declim[0])/dec_inc
           jouter1 = polyfillv(xwindow, ywindow, nra, ndec)
           
           mask_ccr,  maskxf2*tmasks[i].width, maskyf2*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, 0. , tmasks[i].height ;no bow
           xwindow = (ra_mask-ralim[0])/ra_inc
           ywindow = (dec_mask-declim[0])/dec_inc
           jouter2 = polyfillv(xwindow, ywindow, nra, ndec)

           mask_ccr,  maskxf3*tmasks[i].width, maskyf3*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, 0. , tmasks[i].height ;no bow
           xwindow = (ra_mask-ralim[0])/ra_inc
           ywindow = (dec_mask-declim[0])/dec_inc
           jouter3 = polyfillv(xwindow, ywindow, nra, ndec)
       
           mask_ccr,  maskxf4*tmasks[i].width, maskyf4*tmasks[i].height , ra_mask, $
             dec_mask, rac, decc, tmasks[i].pa, 0. , tmasks[i].height ;no bow
           xwindow = (ra_mask-ralim[0])/ra_inc
           ywindow = (dec_mask-declim[0])/dec_inc
           jouter4 = polyfillv(xwindow, ywindow, nra, ndec)
           
           jouter=[jouter1,jouter2,jouter3,jouter4]

;; get rid of any elements in jcore that are NOT in jouter:
           match,jcore,jouter,subcore,subouter
           jcore=jcore[subcore]


; fill in the window function itself
           if keyword_set(nozcorr) then tempwindowf[jcore] = tempwindowf[jcore] + $
             0.634*sfraction else tempwindowf[jcore] = tempwindowf[jcore] + $
             0.634*sfraction*zfactor[0] ;set to 3/4 full value

           if keyword_set(nozcorr) then tempwindowf[jouter] = tempwindowf[jouter] + $
             .216*sfraction else tempwindowf[jouter] = tempwindowf[jouter] + .216*sfraction*zfactor[0]  
;given another chance- adds to core as well giving the core 4* the
     ;sel. probability as the edge, as observed in real masks

       endfor  ; this is the loop over all masks with this one date

       windowf = windowf*(1.d0- (tempwindowf < 1)) ;deal with probabilities properly

   endfor ; this is the loop over different mask dates

   windowf =  1.d0 -windowf     ;limit overlapping zones
   whlow = where(abs(windowf) lt 1.e-6, lowct) ; just to make sure we don't have numerical problems
   if lowct gt 0 then windowf[whlow] = 0.
   return, float(windowf)
end





