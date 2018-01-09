
PRO enk_specs::show_message, msg
	WIDGET_CONTROL, self.statusid, SET_VALUE=msg
END

PRO enk_specs::describe_spectrum

	o = (*self.objects)[self.ix]

	fluxing = ''
	if self.fluxing then fluxing = 'fluxed, '
	type = 'F_counts, '
	if self.flambda then type = 'F_lambda, '
	restframe = 'observed lambda, '
	if self.restframe then restframe = 'restframe lambda, '
	smoothing = ''
	if self.smooth1d[0] then smoothing = 'smoothed, '
	persistent = ''
	if self.persistent then persitent = 'persistent'

	i = strcompress(string(self.ix+1),/remove_all)
	n = strcompress(string(n_elements((*self.objects))),/remove_all)
	obj = '  ['+i+'/'+n+']'
	
	no = ' -- ' + strcompress(o.objno,/rem)

	self->show_message, fluxing+type+restframe+smoothing+persistent+obj+no

	case self.spec_type of
		0: self->obj_info
		1: self->obj_info_treu
		2: self->obj_info_ktrs
	endcase
END

PRO enk_specs::obj_info_ktrs

	o = (*self.objects)[self.ix]
	num_Errors = 0

	catch, error_status
	if error_status ne 0 then begin
		catch, /cancel
		num_errors= num_errors+1
		if num_errors eq 1 then widget_control,self.objinfoid,set_value=str
		return
	endif

	str = ['['+num(o.objno)+','+num(o.maskname)+$
		','+num(o.slitname)+']']
	str = [str, 'ra: ' + dstring(o.ra,6) + ' dec: ' + dstring(o.dec,6)]
	str = [str,'z: '+dstring(o.z,4) + ' B: ' + $
		dstring(o.magb) + ' R: ' + dstring(o.magr) + $
		' I: ' + dstring(o.magI)]
	str = [str, ' vdisp: ' + dstring(o.vdisp) + '+-' + $
		dstring(o.vdisp_err)]
	str = [str, '(U-B)_0: ' + dstring(o.ub,3) + ' M_B ' + $
		dstring(o.mb)]
	str = [str, 'zq: ' + dstring(o.zquality,0) + ' cmt: ', o.comment]

	WIDGET_CONTROL, self.objinfoid, SET_VALUE=str
END

PRO enk_specs::obj_info
	
	o = (*self.objects)[self.ix]

	num_errors=0
	catch, error_status
	if error_status ne 0 then begin
		catch, /cancel
		num_errors= num_errors+1
		if num_errors eq 1 then widget_control,self.objinfoid,set_value=str
		return
	endif

	str = ['['+num(o.objno)+','+num(o.maskname)+$
		','+num(o.slitname)+']']
	str = [str,'z: '+dstring(o.z,4) + ' B: ' + $
		dstring(o.magb) + ' R: ' + dstring(o.magr) + $
		' I: ' + dstring(o.magI)]
	str = [str, 'zq: ' + dstring(o.zquality,0) + ' cmt: ', o.comment]
	str = [str, 'KRB: ' + dstring(o.krb) + ' vdisp: ' + dstring(o.vdisp) + '+-' + $
		dstring(o.vdisp_err)]
	str = [str, '(U-B)_0: ' + dstring(o.ub0,3) + ' M_B ' + $
		dstring(o.abm)]

	oii = measure_simple_ew(*self.ss,$
		[[3696.3,3716.3],[3738.3,3758.3]],$
		[3716.3,3738.3])

	if (size(oii))[0] ne 0 then $
		str = [str, '[OII] ew: '+dstring(oii[0],2)+'+-' + dstring(oii[1],2)]


	WIDGET_CONTROL, self.objinfoid, SET_VALUE=str
	
END

PRO enk_specs::obj_info_treu

        o = (*self.objects)[self.ix]

return
        str = ['['+num(o.objno)+','+num(o.maskname)+$
                ','+num(o.slitname)+']']
        str = [str,'z: '+dstring(o.z,4) + ' Zq: ' + $
                dstring(o.zquality)]
        str = [str,'R: ' + dstring(o.magr) + ' Z: ' + dstring(o.magz)]
        str = [str, '(U-B)_0: ' + dstring(o.ub0,3) + ' M_B ' + $
                dstring(o.abm)]

        str = [str, ' Hb: '+dstring(o.hbeta[0],2) + '+-' + dstring(o.hbeta[1],2)+$
                ' Hg: '+dstring(o.hgamma[0],2) +'+-' + dstring(o.hgamma[1],2)]
        str = [str, ' Hd: '+dstring(o.hdelta[0],2) + '+-' + dstring(o.hdelta[1],2)+$
                ' [OII] ew: '+dstring(o.oii[0],2)+'+-' + dstring(o.oii[1],2)]
	str = [str, ' Dn4000: ' + dstring(o.dn4000[0],2)+ '+-' + dstring(o.dn4000[1],2)]

        oii = measure_simple_ew(*self.ss, [[3696.3,3716.3],[3738.3,3758.3]],$
                [3716.3,3738.3])
        str = [str, ' [OII]: ' + dstring(oii[0]) + '+-'+dstring(oii[1])]

        bi  = (o.hdelta[0]+o.hgamma[0])/2
        str = [str, ' BI: ' + dstring(bi)]

        WIDGET_CONTROL, self.objinfoid, SET_VALUE=str
END
	
PRO enk_specs::show_twod

	if not self.show_twod then return

	o = (*self.objects)[self.ix]
	if not ptr_valid(self.s2d) then begin

		case self.spec_type of
			0: begin	
				slit = read_slit(o,2000,9100) 
			   end
			1: begin
				slit = kp_readslit(o.treuid)
			   end
			2: begin
				slit = read_slit(o, 2000, 10500, /ktrs)
			   end
		endcase

		ptr_free, self.s2d
		self.s2d = ptr_new(slit)
	endif

	slit = (*self.s2d)

	if self.restframe then begin
		slit.lambda = slit.lambda/(1+o.z)
	endif

	nx = (size(slit.flux))[1]
	ny = (size(slit.flux))[2]

	if self.smooth2d[0] then $
		slit.flux = smooth(slit.flux, self.smooth2d[1])

	width = nx/8
	img = dblarr(width,ny*8+5*8)
	for i = 0, 7 do begin
		j = 7-i
		startx = width*i
		endx = width*(i+1)
		l = slit.lambda[startx:endx-1]
		l = [[l],[l],[l],[l],[l]]
		starty = ny*j
		endy = ny*(j+1)+4

		img[*,endy-4:endy] = l

		img[*,starty:endy-5] = slit.flux[startx:endx-1,*]
	end

	self.ds9a->image, img

	for i = 0,7 do begin
		j = 7-i
		self.ds9a->line, 0, ny*j+self.r1, 50, ny*j+self.r1
		self.ds9a->line, 0, ny*j+self.r2, 50, ny*j+self.r2
                self.ds9a->cmd, 'scale limits -50 100'
	endfor
END

PRO enk_specs::show_ktrs, ev

	o = (*self.objects)[self.ix]
	if ptr_valid(self.ktrs) eq 1 then return
	self.ds9b->delete_all_frames
	if o.slitname eq -1 then return
	mask = fix(strmid(o.maskname, 4, 2))
	ktrs_stamps, mask, o.slitname, B,V,I,Z,status=s
	if s lt 0 then return

	renormktrs,b
	renormktrs,v
	renormktrs,i
	renormktrs,z
	
	self.ktrs = ptr_new({B: B, V: V, I: I, Z: Z})

	img = *(self.ktrs)
	B = img.b
	V = img.V
	I = img.I
	Z = img.Z

	self.ds9b->image, B, /frame
	self.ds9b->image, V, /frame
	self.ds9b->image, I, /frame
	self.ds9b->image, Z, /frame
	self.ds9b->rgbimage, v, i, z
END


;                 spec_type -- 0: DEEP, 1: Ellis/Treu, 2: KTRS

PRO enk_specs::show_oned
	COMMON splot_state,state,graphkeys

	tellam = [6860,7600]

	o = (*self.objects)[self.ix]
	if not ptr_valid(self.ss) then begin
		ptr_free, self.ss

		case self.spec_type of
			0: begin
				ss = npk_readspec_2(o.objname, o.maskname, $
					o.slitname, status=s)
			   end
			1: begin
				enk_readtreuspec, o.treuid, lambda=l, spec=s, ivar=i,$
					r1=r1,r2=r2
				self.r1 = r1
				self.r2 = r2
				vactoair, l
				if (size(l))[0] ne 0 then $
					ss = {lambda: l, spec: s, ivar: i}
			   end
			2: begin
				ss = ktrs_read_spec(o.objno)
			   end
		endcase
			
		self.ss = ptr_new(ss)
	endif

	ss = (*self.ss)
	if (size(ss))[0] eq 0 then return

	if self.fluxing then begin
		case self.spec_type of
			0: begin
				dc = deimos_correction(ss.lambda)
				ss.spec = ss.spec*dc
			if has_tag(o, 'flambdafactor') then begin 
				if (size(o.flambdafactor))[0] ne 0 then ss.spec=ss.spec*o.flambdafactor
			endif
			end
			1:
			2:
		endcase
	endif

	if self.flambda then begin
		ml = ss.lambda/mean(ss.lambda)
		ss.spec = ss.spec/ml
	endif

	if self.restframe then begin
		ss.lambda = ss.lambda/(1.+o.z)
		tellam = tellam/(1+o.z)
	endif

	if self.smooth1d[0] ne 0 then $
		ss.spec = poly_smooth(ss.spec, self.smooth1d[1])


	mx = mean(ss.spec[1000:6000])+2*stddev(ss.spec[1000:6000])
	if (size(state))[0] eq 0 then begin
		tely = 1.5
		splot, ss.lambda, ss.spec, yrange=[0,mx]
	endif else begin
		if self.persistent eq 0 then $
			splot,ss.lambda,ss.spec, yrange=[0,mx],psym=10 $
		else $
			splot, ss.lambda, ss.spec, $
			xrange=state.xrange, yrange=state.yrange,psym=10
	endelse

	if self.show_ew ne 0 then $
		show_ews

	tely=state.yrange[1]
	sxyouts, tellam[0], tely, 'B', color=1,charthick=3
	sxyouts, tellam[1], tely, 'A', color=1,charthick=3

	if (size(self.model))[0] ne 0 then begin
		ms = self.model.spec
		ll = self.model.lambda
		ms = ms/median(ms)
		soplot, ll, ms*median(poly_smooth(ss.spec,111)),color=1
	endif
	if has_tag(o, 'bmlambda') then begin
		if ((size(o.bmlambda))[0] ne 0) then begin
			bms = interpol(o.bmspec,o.bmlambda,ss.lambda)
			bms = bms/mean(bms)*mean(ss.spec)
			soplot,ss.lambda,bms,color=1
		endif
	endif

	self->overplot_lines

END

PRO enk_specs::show_ews
	COMMON emission_lines, el_names, el_lambda, el_emission, el_continuum

	setup_emission_lines

	ss = (*self.ss)
	mx = mean(ss.spec[1000:6000])+2*stddev(ss.spec[1000:6000])
	n = n_elements(el_names)
	for i = 0, n-1 do begin
		l = el_emission[0, i]
		soplot, [l,l], [0,mx], color=3
		l = el_emission[1, i]
		soplot, [l,l], [0,mx], color=3
		l = el_continuum[0,0,i]
		ls = (i mod 4) + 1
		soplot, [l,l], [0,mx], color=1, linestyle=ls
		l = el_continuum[1,0,i]
		soplot, [l,l], [0,mx], color=1, linestyle=ls
		l = el_continuum[0,1,i]
		soplot, [l,l], [0,mx], color=1, linestyle=ls
		l = el_continuum[1,1,i]
		soplot, [l,l], [0,mx], color=1, linestyle=ls
	endfor
END

PRO enk_specs::overplot_lines

	n = n_elements(*self.linenames)
	ss = (*self.ss)
	mx = mean(ss.spec[1000:6000])+2*stddev(ss.spec[1000:6000])

	for i = 0, n-1 do begin
		wl = (*self.linewaves)[i]
		wn = (*self.linenames)[i]

		soplot, [wl, wl], [0, mx], color=1
	endfor
		
END


PRO enk_specs::show_web

	if not self.img_type then return
	o = (*self.objects)[self.ix]
	s = 'http://saci.ucolick.org/private/goods/web/objects/' +  $
		string(o.objno,format='(i7.7)')+ '.html'
	
	spawn, 'mozilla -remote "openurl(' + s + ')&"'
END

PRO enk_specs::overplot_sourcecat, ix, cmd, fluxcmd, off=off
	o = (*self.objects)[ix]
	dd = *self.secondary
	mt = (where(o.objno eq dd.d2z_objno,count))[0]
	mt = dd[mt]
	mnra = o.ra-1/60d
	mxra = o.ra+1/60d
	mndc = o.dec-1/60d
	mxdc = o.dec+1/60d
	result = execute('r = where(dd.'+cmd+'ra gt mnra and dd.'+cmd+'ra lt mxra and ' + $
		'dd.'+cmd+'dec gt mndc and dd.'+cmd+'dec lt mxdc, count)')

	if result eq 0 then $
	res = execute('r = where(dd.'+cmd+'_ra_j2000 gt mnra and dd.'+cmd+'_ra_j2000 lt mxra and ' + $
		'dd.'+cmd+'_dec_j2000 gt mndc and dd.'+cmd+'_dec_j2000 lt mxdc, count)')

	print, fluxcmd, count
	if count eq 0 then return

	if result eq 1 then begin
		result = execute('ras = dd[r].'+cmd+'ra')
		result = execute('decs= dd[r].'+cmd+'dec')
		result = execute('flx = dd[r].' + fluxcmd)
		result = execute('thisra = mt.'+cmd+'ra')
		result = execute('thisdec = mt.'+cmd+'dec')
		result = execute('thisflx = mt.' + fluxcmd)
	endif else begin
		result = execute('ras = dd[r].'+cmd+'_ra_j2000')
		result = execute('decs= dd[r].'+cmd+'_dec_j2000')
		result = execute('flx = dd[r].' + fluxcmd)
		result = execute('thisra = mt.'+cmd+'_ra_j2000')
		result = execute('thisdec = mt.'+cmd+'_dec_j2000')
		result = execute('thisflx = mt.' + fluxcmd)
	endelse

	if keyword_Set(off) then off = 0.005 else off = 0.
	

	x = string(thisra) & y = string(thisdec-off) & f = num(thisflx)
	rad = string(3/3600d)
	print,'this: ', x, y, f
	str = 'fk5; circle' + x + y + rad + $
		' # color=yellow text={'+f+'}\n'
	
	if 1 then begin
	for i = 0, count-1 do begin
		x = string(ras[i]) & y = string(decs[i]) 
		lab = num(flx[i]) & rad = string(3/3600d)
		str = str + 'fk5; circle' + x + y + rad+ $
			' # color=red text={'+lab +'}\n'
	endfor
	endif
	self.ds9b->pipecmd, str, 'regions'
END

PRO enk_specs::show_object, ix
	;common splot_state, state, graphkeys


	if ix lt 0 or ix gt n_elements((*self.objects))-1 then begin
		null = dialog_message(['There are not more objects!'])
		return
	endif
	
	widget_control,/hourglass

	self->show_message, 'Working...'
	; Asked to show a new slit, load into memory
	if self.ix ne ix and self.img_type eq 2 then begin
		ptr_free,self.ktrs
	endif
	if self.ix ne ix and (self.img_type eq 0 or self.img_type eq 4) then begin
		ptr_free,self.ss
		ptr_free,self.s2d

		self->show_message, 'Working... Drawing Postage Stamps...'
		o = (*self.objects)[ix]
		pp = getphot(strmid(o.objname, 0, 2))
		r = where(pp.objno eq o.objno)
	;	self.ds9b->cmd, 'frame 13'
	;	self.ds9b->cmd, 'scale limits -40 120'
	;	self.ds9b->cmd, 'scale linear'
	;	ds9_stamp, pp[r], mask=o.maskname, ds9=self.ds9b, stampsz=200, frame=14
	end
	;;; MTR Type
	if self.ix ne ix and self.img_type eq 4 then begin
		o = (*self.objects)[ix]
		print, o.ra, o.dec
	;	egs_stamp, o.objno, self.ds9b
	;	self.ds9b->cmd, 'tile'
	;	self.ds9b->cmd, 'frame 3'
	;	self.ds9b->cmd, 'match frames wcs&'


		if ptr_valid(self.secondary) then begin
		;	self.ds9b->cmd, 'frame 14'
		;	self->overplot_sourcecat, ix, 'd2z_', 'd2z_z'
		;	self.ds9b->cmd, 'frame 1'
		;	self->overplot_sourcecat, ix, 'gx', 'gx_fuv_flux_mj'
		;	self.ds9b->cmd, 'frame 2'
		;	self->overplot_sourcecat, ix, 'gx', 'gx_nuv_flux_mj'
		;	self.ds9b->cmd, 'frame 3'
		;	self->overplot_sourcecat, ix, 'acs_acs', 'acs_acs_viso'
		;	self.ds9b->cmd, 'frame 4'
		;	self->overplot_sourcecat, ix, 'acs_acs', 'acs_acs_iiso'
		;	self.ds9b->cmd, 'frame 7'
		;	self->overplot_sourcecat, ix, 'k_', 'k_kmag_auto'
		;	self.ds9b->cmd, 'frame 8'
		;	;self->overplot_sourcecat, ix, 'i8_', 'i8_irac1'
		;	self->overplot_sourcecat, ix, 'k_', 'k_isubmag_d5[0]'
		;	self.ds9b->cmd, 'frame 9'
		;	;self->overplot_sourcecat, ix, 'i8_', 'i8_irac2'
		;	self->overplot_sourcecat, ix, 'k_', 'k_isubmag_d5[1]'
		;	self.ds9b->cmd, 'frame 10'
		;	;self->overplot_sourcecat, ix, 'i8_', 'i8_irac3'
		;	self->overplot_sourcecat, ix, 'k_', 'k_icrdmag_d5[2]'
		;	self.ds9b->cmd, 'frame 11'
		;	;self->overplot_sourcecat, ix, 'i8_', 'i8_irac4'
		;	self->overplot_sourcecat, ix, 'k_', 'k_icrdmag_d5[3]'
		;	self.ds9b->cmd, 'frame 12'
		;	self->overplot_sourcecat, ix, 'm24_mips_', 'M24_MIPS24'
		endif
	endif
	if self.ix ne ix then begin
		ptr_free, self.ss
	endif


	self.ix = ix
	;self->show_message, 'Working... Showing Website...'
	;self->show_web
	;self->show_message, 'Working... Showing 1d...'
	;self->show_oned
	self->show_message, 'Working... Showing 2d...'
        self.ds9a->cmd, 'scale limits -50 100'
	if self.img_type eq 2 then begin
		self->show_ktrs
	endif else begin
		self->show_twod
	endelse

	self->describe_spectrum
END

PRO enk_specs::toggle_persistence

	self.persistent = not self.persistent
	self->show_object, self.ix
END

PRO enk_specs::toggle_show_ew

	self.show_ew = not self.show_ew
	self->show_object, self.ix
END
	
PRO enk_specs::toggle_fit_model

	if (size(params))[0] eq 0 then begin
		params = {dustn: 0.7, n_waves: 5, dustrange: [-5,20], n_dusts: 53, $
		        n_stretches: 1, stretchrange: [1.,1.]}
	endif

	ss = *self.ss

	res = determine_model(params, ss, bestmodel=bm)
	soplot, bm.lambda, bm.spec,color=1
END

PRO enk_specs::toggle_smoothing
	
	self.smooth1d[0] = not self.smooth1d[0]
	self->show_object, self.ix
END

PRO enk_specs::toggle_smooth2d
	
	self.smooth2d[0] = not self.smooth2d[0]
	self->show_object, self.ix
END

PRO enk_specs::toggle_flambda

	self.flambda = not self.flambda
	self->show_object, self.ix
END

PRO enk_specs::goto_object
	r = where((*self.objects).objno eq self.goto_objno, count)
	if count eq 0 then begin
		print, 'Could not find ', self.goto_objno, ' in the list.'
		return
	endif
	self->show_object, r
END
	
PRO enk_specs::toggle_fluxing

	self.fluxing = not self.fluxing
	self->show_object, self.ix
END

PRO enk_specs::toggle_restframe

	self.restframe = not self.restframe
	self->show_object, self.ix
END

PRO enk_specs::toggle_twod

	self.show_twod = not self.show_twod

	self->show_object, self.ix
END

PRO enk_specs::toggle_web

	self.show_img = not self.show_img

	self->show_object, self.ix
END

PRO enk_specs::toggle_showhideimages

	self.show_img = not self.show_img

	self->show_object, self.ix
END



PRO enk_specs::oned_smooth, ev

	widget_control, ev.id, get_value=val
	self.smooth1d[1] = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(self.smooth1d[1]),/remove_all)

	self->show_object, self.ix
END

PRO enk_specs::twod_smooth, ev

	widget_control, ev.id, get_value=val
	self.smooth2d[1] = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(self.smooth2d[1]),/remove_all)

	self->show_object, self.ix
END

PRO enk_specs::goto_text, ev
	widget_control, ev.id, get_value=val
	self.goto_objno = long(val)
	widget_control, ev.id, set_value = $
		strcompress(string(self.goto_objno),/remove_all)
END

PRO enk_specs::image_command, ev

	widget_control, ev.id, get_value=val
	self.ktrs_cmd = val

	self->show_object, self.ix
END

	
PRO handle_text, ev
	widget_control, ev.id, get_uvalue=uval
	WIDGET_CONTROL, ev.TOP, GET_UVALUE=obj

	case (uval) of
		'smooth1d': obj->oned_smooth, ev
		'smooth2d': obj->twod_smooth, ev
		'goto': obj->goto_text, ev
		'image': 	obj->image_command, ev
	end

END

FUNCTION enk_specs::get_ix
	return, self.ix
END

PRO handle_button, ev

	WIDGET_CONTROL, ev.ID, GET_UVALUE=uval
	WIDGET_CONTROL, ev.TOP, GET_UVALUE=obj

	case (uval) of
		'forward': obj->show_object, obj->get_ix() + 1
		'back': obj->show_object, obj->get_ix() - 1
		'Flux Calibration': obj->toggle_fluxing
		'F_lambda/F_count': obj->toggle_flambda
		'Rest Frame': obj->toggle_restframe
		'smooth2d': obj->toggle_smooth2d
		'smooth1d': obj->toggle_smoothing
		'goto': obj->goto_object
		'Persistent Splot': obj->toggle_persistence
		'Show/Hide EW Measure': obj->toggle_show_ew
		'Fit Model': obj->toggle_fit_model
		'Show/Hide 2d': obj->toggle_twod
		'Show/Hide Web': obj->toggle_web
		'Show/Hide Images': obj->toggle_showhideimg
		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
	ENDCASE
END
;
PRO enk_specs_event, ev
 
	widget_control, ev.id, get_uvalue=value
	if (n_elements(value) eq 0) then value = ''
	name = strmid(tag_names(ev, /structure_name), 7, 4)
	
	case (name) of
		'BUTT': handle_button, ev
		'TEXT': handle_text, ev
    		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
  	ENDCASE
END
PRO enk_specs_cleanup, ev
	WIDGET_CONTROL, ev, GET_UVALUE=obj
	obj_destroy, obj
END

PRO enk_specs::cleanup
	ptr_free, self.objects
	ptr_free, self.model
	ptr_free, self.ss
	ptr_free, self.s2d
	ptr_free, self.linenames
	ptr_free, self.linewaves
	obj_destroy, self.ds9a
	obj_destroy, self.ds9b
	ptr_free, self.ktrs
END

FUNCTION enk_specs::getbase
	return, self.base
END

PRO enk_specs, zz
	n = obj_new('enk_specs', zz)
END

FUNCTION enk_specs::INIT, $	
	objects, $		; Object catalog (zcat for example)
	secondary=secondary, $	; Secondary object catalog 
	Amodel=model, $		; Model Spectrum
	treu=treu,$		; TREU Spectra
	ktrsspec=ktrsspec, $	; KTRS SPEC
	ktrsimg=ktrsimg,$	; KTRS images
	mtr=mtr, $		; MTR Region
	ix=ix 			; Index value to start at

	widget_control,/hourglass
	if keyword_set(model) eq 0 then model =0
	if keyword_set(treu) then self.spec_type = 1
	if keyword_set(ktrsspec) then self.spec_type = 2
	if keyword_set(ktrsimg) then self.img_type = self.img_type or 2
	if keyword_set(ktrsimg) then self.show_img = 1
	if keyword_set(mtr) then begin
		self.img_Type = 4
		self.spec_type = 0
	endif
	if keyword_set(secondary) then begin
		self.secondary = ptr_new(secondary)
	endif
	self.objects = ptr_new(objects)
	self.ix = -1
	if keyword_set(ix) then self.ix = ix

	read_dotlinelist, names, waves
	self.linenames = ptr_new(names)
	self.linewaves = ptr_new(waves)

	; HIGHEST LEVEL
	base = WIDGET_BASE(/column, uvalue=self)
	WIDGET_CONTROL, /MANAGED, base

	; 	Level 1

	b1 = WIDGET_BASE(base, /frame, /row)
	button = WIDGET_BUTTON(b1, VALUE='Back', UVALUE='back')
	button = WIDGET_BUTTON(b1, VALUE='Forward', UVALUE='forward')
	txt = WIDGET_TEXT(b1, XSIZE=80)
	self.statusid = txt


	;	Level 2

	;	Level 2 Left
	b2 = WIDGET_BASE(base, /frame, /row)
	
	tmp = ['Flux Calibration', 'F_lambda/F_count', 'Rest Frame', $
		'Persistent Splot', 'Show/Hide 2d', 'Show/Hide EW Measure', 'Fit Model']
	if self.img_type and 1 then tmp = [tmp, 'Show/Hide Web']
	tmp = [tmp, 'Show/Hide Images']

	XMENU, tmp, b2, /FRAME, uvalue=tmp


	;		Level 2 Middle
	l2_right = WIDGET_BASE(b2, /column, /frame)
	l2_right_top = WIDGET_BASE(l2_right,/row)
	button = WIDGET_BUTTON(l2_right_top, VALUE='Smooth1d', uvalue='smooth1d')
	text2 = WIDGET_TEXT(l2_right_top, XSIZE=5, value='0', $
		/EDITABLE, uvalue='smooth1d')
	l2_right_b = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_b, value='Smooth2d', uvalue='smooth2d')
	text2 = WIDGET_TEXT(l2_right_b, XSIZE=5, value='0', $
		/EDITABLE, uvalue='smooth2d')
	
	l2_right_c = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_c, value='goto', uvalue='goto')
	text2 = WIDGET_TEXT(l2_right_c,XSIZE=9,value='', $
		/editable, uvalue='goto')

	if (self.img_type and 2) eq 2 then begin
		l2_right_c = WIDGET_BASE(l2_right, /row)
		text = WIDGET_TEXT(l2_right_c, value='Image:',XSIZE=6)
		text2 = WIDGET_TEXT(l2_right_c, XSIZE=5,value=self.ktrs_cmd, $
			/EDITABLE, uvalue='image')
	endif

	l2_right_right = WIDGET_BASE(b2, /column, /frame)
	self.objinfoid = WIDGET_TEXT(l2_right_right, xsize=50, $
		ysize=8, value='')

	;	Level 3
	button4 = WIDGET_BUTTON(base, VALUE='Done', UVALUE='DONE')


	if obj_valid(self.ds9a) eq 0 then self.ds9a = obj_new('ds9')
	if obj_valid(self.ds9b) eq 0 then self.ds9b = obj_new('ds9')

	self.ktrs_cmd = 'B'
	self.base = base
	self.restframe = 1
	self.fluxing = 1
	self.flambda = 1
	self.persistent = 0

	WIDGET_CONTROL, base, SET_UVALUE=self
	WIDGET_CONTROL, base, /REALIZE
	XMANAGER, 'enk_specs', self->getbase(), /no_block, cleanup='enk_specs_cleanup'

	self->show_object, 0
	return,1
END



PRO enk_specs__define 
	state = {enk_specs, $
		base: 0L, $		; Base Widget ID
		objects: ptr_new(), $	; List of objects to examine
		secondary: ptr_new(), $	; Secondary object catalog
		ix : -1, $		; the object # being examined
		fluxing : 1, $ 		; Flux calibrate spectrum or not
		flambda: 1, $		; F_lambda or F_counts
		restframe: 1, $		; Show spectrum in restframe
		smooth1d: [1,0], $	; [smooth yes/no, How much smoothing to do in Pixels ]
		smooth2d: [1,0], $	; [smooth yes/no, How much smoothing to do in Pixels ]
		goto_objno: 0L, $		; Object number to jump to
		persistent: 0,$		; Are splot options persistent
		show_ew: 0, $		; Show EW continua and emission regions
		show_twod: 1, $		; Show the twod ATV plot
		model: ptr_new(), $	; A model spectrum to overplot on the dispaly window
		ss: ptr_new(),$		; A pointer to the spectrum being displayed
		s2d: ptr_new(),$	; A pointer to the TWO-D spectrum being displayed
		ktrs: ptr_new(),$	; A pointer to the KTRS Image being displayed
		img_type: 0, $		; bitfield: 1: KTRS WEB, 2: KTRS Image , 4: MTR Image
		spec_type: 0, $		; 0: DEEP, 1: Ellis/Treu, 2: KTRS
		show_img: 1, $		; Show image, dependant on img_type
		objinfoid: 0, $		; Object Info text widget id
		statusid: 0L, $		; Status text ID
		linenames: ptr_new(),$	; Pointer to list of names to display lines for
		linewaves: ptr_new(),$	; Pointer to list of wavelengths to display lines for
		r1: 0, $		; R1 Extraction window
		r2: 0, $		; R2 Extraction window
		ds9a: obj_new(''), $; ds9 spectrum window
		ds9b: obj_new(''),  $; ds9 images window
		ktrs_cmd: 'B'$		; String that contains the command to pass to ATV to display images
	}
END

