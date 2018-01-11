;+
; NAME:
;       setplot
;
; PURPOSE:
;       Set the plotting variables to reasonable values, and set the
;       output device to PostScript.
;
; CALLING SEQUENCE:
;       setplot
;
; KEYWORDS:
;
;-

pro setplot
    cleanplot, /silent
    set_plot, 'PS'
    !P.FONT = 0
    !X.THICK = 3
    !Y.THICK = 3
    !P.THICK = 3
    !P.CHARTHICK = 3
    !P.CHARSIZE = 1.5
    vsym, 24, /fill
    device, /helvetica, font_index=5
    device, /helvetica, /oblique, font_index=18
    device, /times, /italic, font_index=19
    ;device, set_font='Times New Roman', /tt_font
end
