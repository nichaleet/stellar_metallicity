function decticks, axis, index, value
;+
; NAME:
;    raticks
; PURPOSE:
;    Tick format to plot DEC values in plots.  Interfaces with PLOT command.
;
; CALLING SEQUENCE:
;    set PLOT keyword [XY]TICKFORMAT='DECTICKS'
;
; INPUTS:
;    Handled by PLOT
;
; MODIFICATION HISTORY:
;
;       Initial Documentation -- Thu Oct 5 23:03:00 2000, Erik
;                                Rosolowsky <eros@cosmic>
;
;		
;
;-


hour = floor(value)
minute = floor((value - hour)*60.)
sec = round(((value-hour)*60.-minute)*60.)
if sec eq 60 then begin
    sec = 0
    minute++
endif
degsym='!9'+string(176B)+'!X'
minsym='!9'+string(162B)+'!X'
secsym='!9'+string(178B)+'!X'

return, string(hour,  format = '(i02)')+degsym+string(minute,$
  format = '(i02)')+minsym+string(sec, format = "(i02)")+secsym
end
