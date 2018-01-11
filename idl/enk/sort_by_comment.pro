;+
; NAME:
;	sort_by_comment
;
; PURPOSE:
;	To pick out entries in kzcat whose comments match a search
;	string.
;       
; CALLING SEQUENCE:
;	sort_by_comment, searchstring
;
; INPUT:
;	searchstring    string for which to search the comment fields
;                       may include wildcards
;
; OUTPUTS:
;       A kzcat-like structure with only those entries whose comment
;       fields match searchstring.  The output is sorted by mask
;       number and then slit number.
;
; PROCEDURE:
;	Use STRMATCH() to find the strings.
;       
; REVISION HISTORY:
;       enk 2005-10-05
;-

function sort_by_comment, searchstring, ucsc=ucsc

    if keyword_set(ucsc) then zc = ucsc_zcat() else zc = full_zcat()
    w = where(strmatch(zc.comment, searchstring, /fold_case), n)
    if n gt 0 then begin
        print, '# objects with ' + searchstring + ':', n_elements(w)
        list = zc[w]
        list = list[sort(list.slitname)]
        list = list[sort(list.maskname)]
        list.comment = strcompress(list.comment)
        return, list
    endif else begin
        print, 'no objects found'
        return, -1        
    endelse

    ;mn = long(zc[w].maskname)
    ;slit = long(zc[w].slitname)
    ;comment = zc[w].comment
    ;objno = zc[w].objno
    ;list = {objno: objno[0], mask: mn[0], slit: slit[0], comment: comment[0]}
    ;list = replicate(list, n)
    ;for i=1,n-1 do begin
    ;    list[i].objno = objno[i]
    ;    list[i].mask = mn[i]        
    ;    list[i].slit = slit[i]        
    ;    list[i].comment = comment[i]        
    ;endfor
    ;list = list[sort(list.slit)]
    ;list = list[sort(list.mask)]

end
