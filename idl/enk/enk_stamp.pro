;+
; NAME:
;    enk_stamp
;
; PURPOSE:
;    to create an ds9 plot with postage stamps of a given list of objects
;      - note: all objects have to be in the same pointing (ie. field '12')
;
; CATEGORY:
;    photometric analysis
;
;
; CALLING SEQUENCE:
;	enk_stamp, struct,  ds9, stampsz = stampsz, $
;		outfile = outfile, mask=mask
;
; INPUTS:
;    struct = input structure containing list of objects you want postage
;                    stamps of - ex: photo or targ
; 
; OPTIONAL INPUTS:
;    none
;	
; KEYWORD PARAMETERS:
;    stampsz = sets the pixel width of each postage stamp - default 50 pix
;    outfile = string name of fits file to save the composite large image to
;    mask = Maskname, string, used to determine PA of slit
;
; OUTPUTS:
;    creates an atvplot of one large image comprised of all the
;            postage stamps, with optional ellipses drawn for each object
;
;
; MODIFICATION HISTORY:
;    enk 28jul05 - fixed slit length and position
;    npk 15jan05 - Significant rewrite to work with ds9 + provide slitimage
;    alc 28feb02
;    alc 07mar02 - loop through quilts only once
;    alc 22mar02 - fixed sorting of quilts - now works for "targets.fits"
;    alc 08dec02 - added the color and allcolor keywords
;-

pro drawslit, maskname, struct, ds9, yoffset=yoffset
	;catch, error_status
	if not keyword_set(maskname) then return
	if not keyword_set(struct) then return
	error_status = 0
	if error_status ne 0 then begin
		message, 'Could not draw slit for some reason!', /info
		return
	endif
	mn = strcompress(maskname,/rem)
	ip = mrdfits(getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + '.fits', 1)
	mm = mrdfits(getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + '.fits', 2)

	r = where(mm.objno eq struct.objno, count)

	if count eq 0 then begin
		message, 'Could not find ' + string(struct.objno) + ' in photometry file.  Can not overplot slits',/info
		return
	endif

	ra_obj = mm[r].ra & dec_obj = mm[r].dec
	o = mm[r]

	dx = o.slitwidth/3600.
	dy = (o.topdist+o.botdist)/3600.
        pa = o.slitpa+ip.pa
        angle = string(pa,format='(F15.9)')
        pa = pa * !DTOR
        dec_obj_rad = dec_obj * !DTOR

        ra = ra_obj + 0.5*(o.topdist-o.botdist)/3600.*sin(pa)/cos(dec_obj_rad)
        dec = dec_obj + 0.5*(o.topdist-o.botdist)/3600.*cos(pa)
	sn = string(o.slitn)

        ;x1 = ra-0.5*dy*sin(pa)/cos(dec_obj_rad)
        ;y1 = dec-0.5*dy*cos(pa)
        ;x2 = ra+0.5*dy*sin(pa)/cos(dec_obj_rad)
        ;y2 = dec+0.5*dy*cos(pa)
        ;slope_a = (y2-y1)/(x2-x1)
        ;x1 = string(x1,format='(F15.9)')
        ;y1 = string(y1,format='(F15.9)')
        ;x2 = string(x2,format='(F15.9)')
        ;y2 = string(y2,format='(F15.9)')
        ;print, x1 + y1 + x2 + y2
        ;print, 'slope_a: ' + string(slope_a,format='(F15.9)')
        ;ds9->pipecmd, 'fk5; line' + x1 + y1 + x2 + y2 + '# color=blue', 'regions'

        ;x3 = ra+0.5*dx*cos(pa)/cos(dec_obj_rad)
        ;y3 = dec-0.5*dx*sin(pa)
        ;x4 = ra-0.5*dx*cos(pa)/cos(dec_obj_rad)
        ;y4 = dec+0.5*dx*sin(pa)
        ;slope_b = (y4-y3)/(x4-x3)
        ;x3 = string(x3,format='(F15.9)')
        ;y3 = string(y3,format='(F15.9)')
        ;x4 = string(x4,format='(F15.9)')
        ;y4 = string(y4,format='(F15.9)')
        ;print, x3 + y3 + x4 + y4
        ;print, 'slope_b: ' + string(slope_b,format='(F15.9)')
        ;print, 'slope_a * slope_b: ' + string(slope_a*slope_b,format='(F15.9)')
        ;ds9->pipecmd, 'fk5; line' + x3 + y3 + x4 + y4 + '# color=blue', 'regions'

        ra = string(ra,format='(F15.9)')
        dec = string(dec,format='(F15.9)')        
        dx = string(dx,format='(F15.9)')
        dy = string(dy,format='(F15.9)')
        pa = string(pa,format='(F15.9)')
        ;print, 'angle: ' + angle
        ;print, 'pa: ' + pa

	ds9->pipecmd, 'fk5; box' + ra + dec + dx + dy + angle + '# text={' + sn + '} color=red', 'regions'
	;ds9->pipecmd, 'fk5; box' + ra + dec + dx + dy + string(29) + '# color=blue', 'regions'

        ra_obj = string(ra_obj,format='(F15.9)')
        dec_obj = string(dec_obj,format='(F15.9)')
	ds9->pipecmd, 'fk5; point' + ra_obj + dec_obj + '# point=cross color=red', 'regions'
	;ds9->pipecmd, 'fk5; point' + ra + dec + '# point=cross color=green', 'regions'
        
        if keyword_set(yoffset) then begin
            ytoradec, mn, o.objno, yoffset, yra, ydec
            yra = string(yra,format='(F15.9)')
            ydec = string(ydec,format='(F15.9)')
            ds9->pipecmd, 'fk5; point' + yra + ydec + '# point=cross color=green', 'regions'
        endif

END

pro update_astrometry, a, d, ast, ds9
	if (size(ast))[0] eq 0 then begin
		return
	endif
	ast.crpix = [0,0]
	ast.crval = [a,d]
	hdr = strarr(1)
	putast, hdr, ast
	s = strcompress(string(randomu(s,1)),/rem)
	openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
	for i = 0, n_elements(hdr)-1 do begin
		 printf, 1, hdr[i]
	endfor
	close,1
	ds9->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
	file_delete, '/tmp/' + s[0] + 'wcs.tmp'
END

;pro enk_stamp, struct,  ds9, stampsz = stampsz, $
;   outfile = outfile, mask=mask, yoffset=yoffset

pro enk_stamp, objno, ds9, stampsz=stampsz, mask=mask, yoffset=yoffset, outfile=outfile

objnos = strcompress(objno, /rem)
pp = getphot(strmid(objnos,0,2))
struct = pp[where(pp.objno eq objno)]

if NOT keyword_set(stampsz) then stampsz = 100
; pixel width of each postage_stamp


num = n_elements(struct)
nx = fix(sqrt(num))       ; number of objects in x direction
if nx*nx eq num then ny = float(num)/nx else ny = fix(num/nx)+1              
                          ; number of objects in y direction


; loop through objects and find which quilts they are on
quilts = strarr(num)                          ; string array of quilt names
field = strcompress(fix(struct[0].objno/1e6),/remove_all) ; pointing number

for i = 0, num-1 do quilts(i) = quilt_number(struct[i].xs, field, x, /string)

; get array of unique quilt numbers 
quiltsorted = quilts(sort(quilts))
quiltuniq = quiltsorted(uniq(quiltsorted))

; loop through each quilt and make postage stamps of each object on
                     ; that quilt and add to the composite image     

image = intarr(3,stampsz, stampsz) 

k = 0                                         ; counter for object number 

for i = 0, n_elements(quiltuniq)-1 do begin   ; loop through quilts

   	; read in R band quilt image
   	field1 = strmid(field, 0, 1)
   	field2 = strmid(field, 1, 1)
   	path = 'field'+field1+'/'+field2
   	directory=getenv('IDLDEEPPHOTOMETRY_DIR')  + path  + '/quilt/'        

   	rimage=mrdfits(directory +'R.'+  quiltuniq(i)+'.fits.Z',/silent, 0, header)
   	extast, header, rast

	bimage=readfits(directory + 'B.'+quiltuniq(i)+'.fits.Z',/silent, 0, header)
      	extast, header, bast

      	iimage=readfits(directory + 'I.'+quiltuniq(i)+'.fits.Z',/silent, 0, header)
      	extast, header, iast

	; get postage stamps 
   	thisquilt = where(quilts eq quiltuniq(i))
   	objects = struct(thisquilt)                ; objects on this quilt
   	sz = stampsz
   	qltsz = 1024

	for j = 0, n_elements(objects)-1 do begin
                ds9->cmd, 'frame delete all'
                objname = strcompress(struct[j].objno, /rem)
                vfilename = '/net/marc/raid2/egs_stamps/deep2/'+strmid(objname,0,5)+'/'+objname+'.acs.v.fits.gz'
                print, vfilename
                if file_test(vfilename) then begin
                    print, 'here'
                    stamp = mrdfits(vfilename, 0, hdr)
                    xstart = 0
                    ystart = 0
                    ds9->image,stamp,/frame
                    ds9->cmd, 'scale mode zscale'
                    s = strcompress(string(randomu(s,1)),/rem)
                    openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
                    for i = 0, n_elements(hdr)-1 do begin
                        printf, 1, hdr[i]
                    endfor
                    close,1
                    ds9->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
                    file_delete, '/tmp/' + s[0] + 'wcs.tmp'

                    if keyword_set(mask) and keyword_set(yoffset) then drawslit, mask, struct, ds9, yoffset=yoffset $
                    else if keyword_set(mask) then drawslit, mask, struct, ds9, yoffset=yoffset
                endif

                ifilename = '/net/marc/raid2/egs_stamps/deep2/'+strmid(objname,0,5)+'/'+objname+'.acs.i.fits.gz'
                if file_test(ifilename) then begin
                    stamp = mrdfits(ifilename, 0)
                    xstart = 0
                    ystart = 0
                    ds9->image,stamp,/frame
                    ds9->cmd, 'scale mode zscale'
                    s = strcompress(string(randomu(s,1)),/rem)
                    openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
                    for i = 0, n_elements(hdr)-1 do begin
                        printf, 1, hdr[i]
                    endfor
                    close,1
                    ds9->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
                    file_delete, '/tmp/' + s[0] + 'wcs.tmp'
                endif

      		tmp = quilt_number(objects[j].xs, field, x, /string)
      		x[0] = round(x[0]) & x[1] = round(x[1])
		xp = x-sz/2-1	; xp is in the coordinates of the astrometry

      		; put quilt in a larger image to take care of cases where stamp 
             	; hits edge of quilt

         	largeim = fltarr(qltsz+sz,qltsz+sz)
         	largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = bimage
		
      		; extract the B stamp
         	stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
         	xstart = 0 & ystart = 0
	 	ds9->image,stamp,/frame
                ds9->cmd, 'scale mode zscale'
                ds9->cmd, 'tile'
                ds9->cmd, 'zoom to 4'
                ds9->cmd, 'wcs fk5'        
                ds9->cmd, 'wcs align yes'
                ds9->cmd, 'wcs format degrees'
		xy2ad, xp[0]-1, xp[1]-1, bast, a, d
		update_astrometry, a, d, bast, ds9

         	largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = rimage
      		; extract the R stamp
         	stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
         	xstart = 0 & ystart = 0
	 	ds9->image,stamp,/frame
                ds9->cmd, 'scale mode zscale'
		xy2ad, xp[0]-1, xp[1]-1, rast, a, d
		update_astrometry, a, d, rast, ds9
		if keyword_set(mask) and keyword_set(yoffset) then drawslit, mask, struct, ds9, yoffset=yoffset $
                else if keyword_set(mask) then drawslit, mask, struct, ds9, yoffset=yoffset

         	largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = iimage
      		; extract the I stamp
         	stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
         	xstart = 0
         	ystart = 0
	 	ds9->image,stamp,/frame
                ds9->cmd, 'scale mode zscale'
		xy2ad, xp[0]-1, xp[1]-1, iast, a, d
		update_astrometry, a, d, iast, ds9                
                ds9->cmd, 'match frames wcs'

         	k = k+1
   	endfor
endfor

	k = 0  
	for i = 0, n_elements(quiltuniq)-1 do begin   ; loop through quilts
		thisquilt = where(quilts eq quiltuniq(i))
         	objects = struct(thisquilt)     

         	; loop through objects on this quilt
         	for j = 0, n_elements(objects)-1 do begin
            		xoffset = (sz-1)/2.+sz*(k-nx*fix(k/float(nx)))
            		yoffset = (sz-1)/2.+sz*fix(k/float(nx))
	    		a = 3*objects[j].rg
	    		e2 = objects[j].e2
	    		b = (1-e2)*a
            		ds9->ellipse, xoffset, yoffset, theta=objects[j].pa, b,a
            	k = k+1
		endfor
	endfor      


	if keyword_set(outfile) then writefits, outfile, image
end

