;+
; NAME: enk_mosfire_add_truitime
;
; PURPOSE: Add TRUITIME header to MOSFIRE images taken June 9-15,
; 2014.
;
; CALLING SEQUENCE: mosfire_add_truitime
;
; INPUTS: none
;
; OUTPUTS: none
;
; PROCEDURE: Run within the raw data directory.  It will add the
; correct header to all m*.fits files.  TRUITIME is calculated from
; the ITIME header keyword.
;
; MODIFICATION HISTORY:
;       2014-06-24 ENK created
;       2014-06-26 ENK fixed bugs
;-

pro enk_mosfire_add_truitime
    f = file_search('m*.fits', count=c)
    for i=0,c-1 do begin
        hdr = headfits(f[i])
        itime = sxpar(hdr, 'ITIME') / 1000.
        truitime = floor(itime/1.45479)*1.45479
        fxhmodify, f[i], 'TRUITIME', truitime
    endfor
end
