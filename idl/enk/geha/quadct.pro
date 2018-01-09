;--------------------------------------------------------------------
;+
; NAME:
;
; 
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; MODIFICATION HISTORY:
;~~~mg
;- 
;--------------------------------------------------------------------

pro quadct, x, y, xx, yy, nn, fa, fb, fc, fd

    na=0
    nb=0
    nc=0
    nd=0

    wg = where(yy gt y, cwg, complement=wl)
    cwl = nn - cwg
    if cwg gt 0 then begin
       ww = where(xx[wg] gt x[wg], na)
       nb = cwg - na
       ;ww = where(xx[wg] le x[wg], nb)
    endif
    if cwl gt 0 then begin
       ww = where(xx[wl] gt x[wl], nd)
       nc = cwl - nd
       ;ww = where(xx[wl] le x[wl], nc)
    endif

    ff=1.0/nn
    fa=ff*na
    fb=ff*nb
    fc=ff*nc
    fd=ff*nd

 end
