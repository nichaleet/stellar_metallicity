;+
; NAME:
;	uJtoABerr
;
; PURPOSE:
;       Convert an error in microJanskys to an error in
;       AB magnitude.
;       
; CALLING SEQUENCE:
;	aberr = ujtoaberr(ujerr, uj)
;
; INPUT:
;	ujerr    error in microJanskys
;	uj       measurement in microJanskys
;
; RESULT:
;	error in AB magnitude
;-

function uJtoABerr, uJerr, uJ
    ABerr = 2.5d/alog(10d) * uJerr/uJ
    return, ABerr
end
