PRO WFOCUS, index
;+
; NAME:
;	WFOCUS
; PURPOSE:
;	To shift the focus from the currently active window to an
;	open, previously used window, allowing the user to overplot in
;	the previous window.  Loads up the axis parameters saved by
;	'WNEW.PRO' when it switches focus.
; CALLING SEQUENCE:	
;	WFOCUS, Index
; INPUT:
;	Index   The index number of the window to focus on
; OUTPUTS:
;	None.
; PROCEDURE:
;	Checks to see if the window structure system variable exists,
;	and if not creates one.  
;       Stores the x, y and z axis parameters of the current window.
;       Switches focus to the existing window, and loads up the
;       corresponding axis parameters from the window structure. 
;       NOTE: If 'WNEW.PRO' was not used to create the current window,
;       the axis parameters of the old window will not be stored and
;       will therefore be unavailable to load.
; COMMON BLOCKS:
;	None.
;-

MAXWINDOWS = 32          ; Maximum allowed number of windows

if N_ELEMENTS(index) eq 0 then begin
    index = 0
endif

if !D.WINDOW eq -1 then begin
    print, 'No windows!'
    return
endif

; If window structure system variable does not exist, create it.

DEFSYSV, '!WIN_STRUCT', EXISTS = wstruct_exists

if (wstruct_exists eq 0) then begin

    temp = {windowstruct, xaxis:!X, yaxis:!Y, zaxis:!Z, plotv:!P}
    DEFSYSV, '!WIN_STRUCT', REPLICATE({windowstruct}, MAXWINDOWS)

endif

; Store info from current window before shifting focus

!WIN_STRUCT[!D.WINDOW].xaxis = !X
!WIN_STRUCT[!D.WINDOW].yaxis = !Y
!WIN_STRUCT[!D.WINDOW].zaxis = !Z
!WIN_STRUCT[!D.WINDOW].plotv = !P

; Shift focus to the requested window

wset, index

; Set x and y axes on refocused window to match their stored values

!X = !WIN_STRUCT[index].xaxis
!Y = !WIN_STRUCT[index].yaxis
!Z = !WIN_STRUCT[index].zaxis
!P = !WIN_STRUCT[index].plotv

end
