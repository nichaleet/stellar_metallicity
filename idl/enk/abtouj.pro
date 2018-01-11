;+
; NAME:
;	ABtouJ
;
; PURPOSE:
;       Convert AB magnitudes to microJanskys.
;       
; CALLING SEQUENCE:
;	uj = abtouj(ab)
;
; INPUT:
;	ab       measurement in AB magnitudes
;
; RESULT:
;	measurement in uJ
;-

function ABtouJ, AB
    uJ = 10^(-(AB+48.6d)/2.5d) * 1d29
    return, uJ
end
