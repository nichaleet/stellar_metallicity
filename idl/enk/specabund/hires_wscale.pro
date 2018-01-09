function hires_wscale, $
                      chip=chip    ; options
                      
; Returns a generic wavelength solution for Keck-HIRES.
; This routine works for "pre-upgrade" spectra (k-series)
; as well as "post-upgrade" spectra (j-series).
; chip = "blue", "middle" [="iod"], "red" [="ir"], "k" (for k-series spectra)
;     blue is default
; Usage: IDL> bscale = hires_wscale(chip='blue')
;        IDL> help, bscale
;        BSCALE          FLOAT     = Array[4021, 23]
;        IDL> print,wav(0,10)
;              4061.53
;        So, the 0th pixel in order 10 has a wavelength of 4061.53 Angstroms
;        Note that a wavelength shift, Dlam, corresponds to a velocity, v, through:
;             Dlam/lam = z = v/c
;
; This wavelength scale was calibrated with one particular Th-Ar spectrum
; and should be accurate to +- 1 pixel for any post-upgrade HIRES spectrum.

dir = getenv('CALTECH')+'hires/ahoward/'
if keyword_set(chip) then begin
    case chip of
        'blue':   filename = 'keck_bwav.dat'
        'middle': filename = 'keck_rwav.dat'
        'iod':    filename = 'keck_rwav.dat'        ; same as 'middle'
        'red':    filename = 'keck_iwav.dat' 
        'ir':     filename = 'keck_iwav.dat'        ; same as 'red'
        'k':      filename = 'keck_kseries_wav.dat' ; for k-series spectra
    endcase
endif else begin
    filename = 'keck_bwav.dat'
endelse
restore,dir+filename
return,wav

end ; function

