
;
;	ds9 uses the XPA library to communicate with the excellent
;		ds9 display program.  This is a replacement for ATV
;
;	Requirements:
;		1) A working XPA library, binaries must be in the shell's
;			PATH.
;		2) A working DS9 binary.  
;
;	Code
;		1) Checks to see if ds9 has been started
;		1b) If not started, spawns ds9 and loops for 15sec until
;			ds9 stars or timeouts.  Does not handle the case
;			of ds9 not starting!!  (BUG)
;		2) Stores the old shared-memory segment into a temporary
;			variable.  Then, allocates a new shared-memory
;			segment.
;		3) Via XPA, tells ds9 the id number of the shared mem
;			segment.
;		4) Unmaps the old shared memory segment.
;

PRO ds9enk::set_frame, ix
	self->cmd, 'frame frameno ' + string(ix)
END 

PRO ds9enk::update_wcs, hdr
	s = strcompress(string(randomu(s,1)),/rem)
	openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
	for i = 0, n_elements(hdr)-1 do begin
		printf, 1, hdr[i]
	endfor
	close,1
	self->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
	file_delete, '/tmp/' + s[0] + 'wcs.tmp'
END

PRO ds9enk::new_frame
	self->cmd, 'frame new'
END 

PRO ds9enk::cleanup
	catch, error_status
	if error_status ne 0 then begin
		catch, /cancel
		return
	endif
	self->cmd, 'quit'
	if self.segname ne '' then $
		shmunmap, self.segname
END

PRO ds9enk::delete_all_frames
	self->cmd, 'frame delete all'
END 

PRO ds9enk::zoom_to_fit
	self->cmd, 'zoom to fit'
END

PRO ds9enk::line, x0, y0, x1, y1, color=color, text=text
	
	if self.os eq -1 then return
	if keyword_set(text) eq 0 then text = ""

	if keyword_set(color) eq 0 then color = 'blue'
	cmd = 'physical; line ' + string(x0) + ' ' + string(y0)  + ' ' + string(x1) + $
		string(y1) + ' # color = ' + color + " text = {" + text + "}"

	self->pipecmd, cmd, 'regions'
end

PRO ds9enk::ellipse, x0, y0, rx, ry, theta=theta, color=color
	
	if keyword_set(color) eq 0 then color = 'blue'

	angle = (keyword_set(theta) ? string(theta) : '0')
	cmd = 'physical; ellipse' + string(x0) + ' ' + string(y0)  + ' ' + string(rx) + ' ' + $
		string(ry) + ' ' + angle + ' # color = ' + color

	self->pipecmd, cmd, 'regions'
end

; Take a bunch of segments  connect
; Take a bunch of segments  connect
PRO ds9enk::lineconnected, x, y, color=color
	
	if (size(self.os))[1] eq 0 then return

	if n_elements(x) ne n_elements(y) then begin
		message, 'The length of x and y vectors is not the same'
	endif

	if keyword_set(color) eq 0 then color = 'blue'

	for i = 0, n_elements(x) - 2 do begin
		x0 = x[i] & y0 = y[i]
		x1 = x[i+1] & y1 = y[i+1]
		ds9enkline, x0, y0, x1, y1, color=color
	endfor
end


PRO ds9enk::pipecmd, pipe, cmd
	spawn, 'echo "' + pipe+ '" | ' + self.path + 'xpaset ' + self.title + ' ' + cmd
END

PRO ds9enk::cmd,cmd 
	spawn, self.path + 'xpaset -p ' + self.title + ' ' + cmd
END

PRO ds9enk::spawn_ds9enk, geometry=geometry, _extra=extra
	p = self.path
	spawn, p + 'xpaaccess ' + self.title, res

	if keyword_set(geometry) eq 0 then geometry = '' $
	else geometry = ' -geometry ' + geometry
	if res[0] eq 'no' then begin
		print, 'Spawning ds9enk with title: ' + self.title
		spawn, p+'xpans&'
		; -- 1b
		spawn, p + 'ds9 -title ' + self.title + ' -port 0' + geometry + '&'
		for i = 0, 30 do begin
			wait, .5
			spawn, p + 'xpaaccess ' + self.title, res
			if res[0] eq 'yes' then break
		endfor
	endif
END



FUNCTION ds9enk::INIT, _EXTRA=_extra
	self.path = '/opt/NewDS9/3.0.2/saods9/bin/'
	self.os = -1
	self.segname = 'none'
	self.title = 'ds9_' + strcompress(string(systime(/sec),format='(i13)'),/rem)

	; -- 1
	self->spawn_ds9enk, _extra=_extra
	return,1
END

PRO ds9enk::rgbimage, $
	red, $
	green, $
	blue

	self->cmd, 'rgb'
	self->image, red
	self->cmd, 'rgb green'
	self->image, green 
	self->cmd, 'rgb blue'
	self->image, blue 
END

PRO ds9enk::image, $
	img, $	; The image to display 
	preserve=preserve, $	; Preserve scale during load
	frame=frame 	; Keyword set if image to display in new frame

	self->spawn_ds9enk

	p = self.path
	; -- 2
	oldseg = '0'
	if self.os ne -1 then begin
		oldseg = self.segname
	endif
	shmmap, template=img, get_name=segname,/sysv, get_os_handle=os
	self.segname = segname
	self.os = os
	var = shmvar(self.segname)

	var[0,0] = img
	sz = size(img,/dim)

	; -- 3
	if keyword_set(frame) then begin
		spawn, p + 'xpaset -p ' + self.title + ' frame new'
	endif
	if keyword_set(preserve) then begin
		spawn, p + 'xpaset -p ' + self.title + ' preserve scale yes'
	endif else begin
		spawn, p + 'xpaset -p ' + self.title + ' preserve scale no'
	endelse
	spawn, p + 'xpaset -p ' +  self.title + ' shm array shmid ' + strtrim(self.os,2) + ' "[xdim='+$
		strtrim(sz[0],2)+',ydim='+strtrim(sz[1],2)+',bitpix=-32]"&'

	; -- 4
	if oldseg ne '0' then begin
		shmunmap, oldseg
	endif

END

PRO ds9enk__define
	params = {DS9ENK, $
		title: '', $
		os: -1L, $
		segname: '', $
		path: '' $
	}
END
