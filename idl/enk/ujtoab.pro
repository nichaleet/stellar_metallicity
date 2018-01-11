;+
; NAME:
;	uJtoAB
;
; PURPOSE:
;       Convert microJanskys to AB magnitudes.
;       
; CALLING SEQUENCE:
;	ab = ujtoab(uj)
;
; INPUT:
;	uj       measurement in microJanskys
;
; RESULT:
;	measurement in AB magnitudes
;-

function uJtoAB, uJ
    AB = 2.5d*(29d - alog10(uJ)) - 48.6d
    return, AB
end
