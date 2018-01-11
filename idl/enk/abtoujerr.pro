;+
; NAME:
;	ABtouJerr
;
; PURPOSE:
;       Convert an error in AB magnitudes to an error in
;       microJanskys.
;       
; CALLING SEQUENCE:
;	ujerr = abtoujerr(aberr, ab)
;
; INPUT:
;	aberr    error in AB magnitudes
;	ab       measurement in AB magnitudes
;
; RESULT:
;	error in uJ
;-

function ABtouJerr, ABerr, AB
    uJerr = alog(10d)/2.5d * ABtouJ(AB) * ABerr
    return, uJerr
end
