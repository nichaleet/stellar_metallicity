;+
; NAME:
;	complement
;
; PURPOSE:
;	Given an array of indices and a number of elements n, return
;	which integers (through n) are not in the array.  This is useful
;	to use with LISTMATCH.
;       
; CALLING SEQUENCE:
;	result = complement(array, n)
;
; INPUT:
;	array   integer array of which to find the complement
;	n       highest integer to return in the complement array
;
; RESULT:
;	integer array of the complement
;
; PROCEDURE:
;	Very low-tech.  Search through the input array and add each
;	element not found to the output array, one by one.
;       
; REVISION HISTORY:
;       enk 2005-09-13
;-

function complement, array, n, count
    count = 0
    go = 0
    compare = lindgen(n)
    for i=0L, n-1 do begin
        junk = where(array eq compare[i], countw)
        if countw eq 0 then begin
            count = count + 1
            if go then result = [result, compare[i]] $
            else begin
                result = compare[i]
                go = 1
            endelse
        endif
    endfor
    if go eq 0 then result = -1
    return, result
end
