;+
; NAME:
;	sex2deg
;
; PURPOSE:
;	Converts sexagesimal celestial coordiantes to decimal degree
;	coordinates.
;       
; CALLING SEQUENCE:
;	result = sex2dec([[ralist], [declist]])
;
; INPUT:
;	inascii  string array of size [N,2], where the RA list is [N,0]
;	         and the Dec list is [N,1]
;
; RESULT :
;	[N,2] array of decimal coordinates in degrees
;
; REVISION HISTORY:
;       enk 2006-11-17
;-

function sex2deg, inascii
    if (size(inascii))[2] ne 2 then begin
        message, 'Size of input array must be [N,2].', /info
        return, [-1, -1]
    endif
    asciira = strtrim(inascii[*,0])
    hr = double(strmid(asciira, 0, 2))
    min = double(strmid(asciira, 3, 2))
    sec = double(strmid(asciira, 6))
    ra = (hr + min/60d + sec/3600d) * 360d / 24d
    asciidec = strtrim(inascii[*,1])
    deg = double(strmid(asciidec, 0, 3))
    arcmin = double(strmid(asciidec, 4, 2))
    arcsec = double(strmid(asciidec, 7))
    dec = deg + arcmin/60d + arcsec/3600d
    return, [[ra], [dec]]
end
