;+
; NAME:
;       resetplot
;
; PURPOSE:
;       Set the plotting variables to their defaults, and set the
;       output device to X Windows.
;
; CALLING SEQUENCE:
;       resetplot
;-

pro resetplot
    cleanplot, /silent
    set_plot, 'X'
end
