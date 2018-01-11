pro ser_specs::show_message, msg
	WIDGET_CONTROL, self.statusid, SET_VALUE=msg
end


pro ser_specs::show_object, ix
    if ix lt 0 or ix gt n_elements((*self.cat))-1 then begin
        null = dialog_message(['There are not more objects!'])
        return
    endif
    
    widget_control, /hourglass

    ;self->show_message, 'Working...'
    ; Asked to show a new slit, load into memory
    if self.ix ne ix then begin
        self->show_message, 'Drawing Postage Stamps...'
        o = (*self.cat)[ix]
        ;pp = getphot(strmid(strcompress(o.tgtobjno, /rem), 0, 2))
        ;r = where(pp.objno eq o.tgtobjno)
        self.ds9->cmd, 'frame 0'
        self.ds9->cmd, 'scale limits -40 120'
        self.ds9->cmd, 'scale linear'
        ser_stamp, o, self.ds9, stampsz=200

        ;print, o.ra, o.dec
        ;WIDGET_CONTROL, self.fieldsid, get_VALUE=msg
        ;str = ""
        ;for i = 0,n_elements(msg)-1 do str = str + msg[i]
        ;egs_stamp, o.tgtobjno, self.ds9
        self.ds9->cmd, 'tile'
        junk = where(o.cfh_type eq 2 or o.cfh_type eq 3, ccfh)
        junk = where(o.acs_tile ne 0 and o.acs_id ne 0, cacs)
        ;self->show_message,'['+strcompress(ix+1,/rem)+'/'+strcompress(n_elements((*self.cat)),/rem)+'] Target: '+strcompress(o.tgtobjno, /rem)+', '+strcompress(ccfh, /rem)+' CFH serendips, '+strcompress(cacs, /rem)+' ACS serendips'
        self->show_message,'['+strcompress(ix+1,/rem)+'/'+strcompress(n_elements((*self.cat)),/rem)+'] Target: '+strcompress(o.tgtobjno, /rem)
        WIDGET_CONTROL, self.cfhtext, SET_VALUE=strcompress(ccfh, /rem)
        WIDGET_CONTROL, self.acstext, SET_VALUE=strcompress(cacs, /rem)
    endif
    self.ix = ix
end


pro ser_specs::goto1_object
	r = where((*self.cat).tgtobjno eq self.goto1_objno, count)
	if count eq 0 then begin
		print, 'Could not find target object number ' + strcompress(self.goto1_objno, /rem) + ' in the list.'
		return
	endif

	if count gt 1 then begin
		res = dialog_message('Multiple matches to this target.  Presenting the first in the catalog.',/information)
		r = r[0]
	endif
	self->show_object, r
end

pro ser_specs::goto1_text, ev
	widget_control, ev.id, get_value=val
	self.goto1_objno = long(val)
	widget_control, ev.id, set_value = $
		strcompress(string(self.goto1_objno),/remove_all)
end

pro ser_specs::goto2_object
	r = self.goto2_objno
        if r lt 1 or r gt n_elements(*self.cat) then begin
		print, 'Object index ' + strcompress(self.goto2_objno, /rem) + ' is out of range.'
		return
	endif
	self->show_object, r-1
end

pro ser_specs::goto2_text, ev
	widget_control, ev.id, get_value=val
	self.goto2_objno = long(val)
	widget_control, ev.id, set_value = $
		strcompress(string(self.goto2_objno),/remove_all)
end

	
pro handle_text, ev
	widget_control, ev.id, get_uvalue=uval
	WIDGET_CONTROL, ev.TOP, GET_UVALUE=obj

	if keyword_set(uval) eq 0 then return
	case (uval) of
		'goto1': obj->goto1_text, ev
		'goto2': obj->goto2_text, ev
		else: print, uval
	end

end

function ser_specs::get_ix
	return, self.ix
end

pro handle_button, ev

	WIDGET_CONTROL, ev.ID, GET_UVALUE=uval
	WIDGET_CONTROL, ev.TOP, GET_UVALUE=obj

	case (uval) of
		'forward': obj->show_object, obj->get_ix() + 1
		'back': obj->show_object, obj->get_ix() - 1
		'goto1': obj->goto1_object
		'goto2': obj->goto2_object
		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
	endCASE
end
;
pro ser_specs_event, ev
 
	widget_control, ev.id, get_uvalue=value
	if (n_elements(value) eq 0) then value = ''
	name = strmid(tag_names(ev, /structure_name), 7, 4)
	
	case (name) of
		'BUTT': handle_button, ev
		'TEXT': handle_text, ev
    		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
  	endCASE
end
pro ser_specs_cleanup, ev
	WIDGET_CONTROL, ev, GET_UVALUE=obj
	obj_destroy, obj
end

pro ser_specs::cleanup
	ptr_free, self.cat
	obj_destroy, self.ds9
end

function ser_specs::getbase
	return, self.base
end

pro ser_specs, zz
	n = obj_new('ser_specs', zz)
end


function ser_specs::INIT, $	
        cat, $           ; Object catalog (zcat for example)
	ix=ix 			; Index value to start at

    widget_control,/hourglass
    if keyword_set(secondary) then begin
        self.secondary = ptr_new(secondary)
    endif
    self.cat = ptr_new(cat)
    self.ix = -1
    if keyword_set(ix) then self.ix = ix

    ; HIGHEST LEVEL
    base = WIDGET_BASE(/column, uvalue=self)
    WIDGET_CONTROL, /MANAGED, base

    ; 	Level 1

    b1 = WIDGET_BASE(base, /frame, /row)
    button = WIDGET_BUTTON(b1, VALUE='Back', UVALUE='back')
    button = WIDGET_BUTTON(b1, VALUE='Forward', UVALUE='forward')
    txt = WIDGET_TEXT(b1, XSIZE=77)
    self.statusid = txt


    ;	Level 2
    b2 = WIDGET_BASE(base, /frame, /row)
    l2_goto1 = WIDGET_BASE(b2, /column, /frame)
    l2_goto1_c = WIDGET_BASE(l2_goto1, /row)
    button = WIDGET_BUTTON(l2_goto1_c, value='goto objno', uvalue='goto1')
    text1 = WIDGET_TEXT(l2_goto1_c,XSIZE=9,value='', $
                        /editable, uvalue='goto1')

    l2_goto2 = WIDGET_BASE(b2, /column, /frame)
    l2_goto2_c = WIDGET_BASE(l2_goto2, /row)
    button = WIDGET_BUTTON(l2_goto2_c, value='goto', uvalue='goto2')
    text2 = WIDGET_TEXT(l2_goto2_c,XSIZE=3,value='', $
                        /editable, uvalue='goto2')
    gotolabel = WIDGET_LABEL(l2_goto2_c, VALUE='/'+strcompress(n_elements(*self.cat),/rem), uvalue='nlabel')

    l2_cfhser = WIDGET_BASE(b2, /column, /frame)
    l2_cfhser_c = WIDGET_BASE(l2_cfhser, /row)
    cfhlabel = WIDGET_LABEL(l2_cfhser_c, VALUE='CFH serendips:', UVALUE='cfhlabel')
    cfhtext = WIDGET_TEXT(l2_cfhser_c, XSIZE=3, value='', uvalue='cfhtext')
    self.cfhtext = cfhtext

    l2_acsser = WIDGET_BASE(b2, /column, /frame)
    l2_acsser_c = WIDGET_BASE(l2_acsser, /row)
    acslabel = WIDGET_LABEL(l2_acsser_c, VALUE='ACS serendips:', UVALUE='acslabel')
    acstext = WIDGET_TEXT(l2_acsser_c, XSIZE=3, value='', uvalue='acstext')
    self.acstext = acstext

    ;	Level 3
    button4 = WIDGET_BUTTON(base, VALUE='Done', UVALUE='DONE')

    if obj_valid(self.ds9) eq 0 then self.ds9 = obj_new('ds9')

    self.ktrs_cmd = 'B'
    self.base = base

    WIDGET_CONTROL, base, SET_UVALUE=self
    WIDGET_CONTROL, base, /REALIZE
    XMANAGER, 'ser_specs', self->getbase(), /no_block, cleanup='ser_specs_cleanup'

    self->show_object, 0
    return,1
end

pro ser_specs__define 
    state = {ser_specs, $
             base: 0L, $		; Base Widget ID
             cat: ptr_new(), $	; List of objects to examine
             secondary: ptr_new(), $	; Secondary object catalog
             ix : -1L, $		; the object # being examined
             goto1_objno: 0L, $		; Object number to jump to
             goto2_objno: 0L, $		; Object number to jump to
             ktrs: ptr_new(),$ ; A pointer to the KTRS Image being displayed
             show_img: 1, $        ; Show image, dependent on img_type
             objinfoid: 0, $		; Object Info text widget id
             statusid: 0L, $		; Status text ID
             cfhtext: 0L, $		; Status text ID
             acstext: 0L, $		; Status text ID
             ds9: obj_new(''),  $   ; ds9 images window
             ktrs_cmd: 'B'$ ; String that contains the command to pass to ATV to display images
            }
end
