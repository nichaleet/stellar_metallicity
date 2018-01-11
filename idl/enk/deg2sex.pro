;+
; NAME:
;	deg2sex
;
; PURPOSE:
;	Converts decimal degree celestial coordinates to sexagesimal
;	coordiantes.
;       
; CALLING SEQUENCE:
;	result = deg2sex([[ralist], [declist]])
;
; INPUT:
;	incoords  array of size [N,2], where the RA list is [N,0] and
;	          the Dec list is [N,1]
;
; RESULT:
;	[N,2] string array of sexigesimal coordinates in hours (RA)
;	and degrees (Dec)
;
; REVISION HISTORY:
;       enk 2006-11-18
;-

function deg2sex, incoords
    if (size(incoords))[2] ne 2 then message, 'Size of input array must be [N,2].'
    w = where(incoords[*,0] lt 0 or incoords[*,0] gt 360 or incoords[*,1] lt -90 or incoords[*,1] gt 90, cw)
    if cw gt 0 then message, 'RA or Dec out of allowed range.'
    coordsra = incoords[*,0] * 24d / 360d
    hr = floor(coordsra)
    min = floor((coordsra - hr)*60d)
    sec = ((coordsra - hr)*60d - min)*60d
    ra = string(hr, format=('(I02)'))+':'+string(min, format='(I02)')+':'+string(sec, format='(D06.3)')
    wp = where(incoords[*,1] ge 0, complement=wn)
    signstr = strarr((size(incoords))[1])
    signstr[wp] = '+'
    signstr[wn] = '-'
    coordsdec = abs(incoords[*,1])
    deg = floor(coordsdec)
    arcmin = floor((coordsdec - deg)*60d)
    arcsec = ((coordsdec - deg)*60d - arcmin)*60d
    dec = signstr+string(deg, format=('(I02)'))+':'+string(arcmin, format='(I02)')+':'+string(arcsec, format='(D06.3)')
    return, [[ra], [dec]]
end
