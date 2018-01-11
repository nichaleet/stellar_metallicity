PRO WNEW, index
;+
; NAME:
;	WNEW
; PURPOSE:
;	Make a plot window, but store information from the previous
;	window so that it can be recalled if the focus is transferred
;	back to the original window.
;       To be used with WFOCUS.PRO
; CALLING SEQUENCE:	
;	WNEW, Index
; INPUT:
;	Index   The index number of the newly created window, by which
;	        IDL can identify it.
; OUTPUTS:
;	None.
; PROCEDURE:
;	Creates a structure in which to store information from the
;	current window, if it does not exist already.
;       Stores the x, y and z axis parameters of the current window in
;       the window structure. 
;       Creates a new window with the specified index.
; COMMON BLOCKS:
;	None.
;-

MAXWINDOWS = 32          ; Maximum allowed number of windows

; If index is not specified, assign an index for the next window: 0 if
;    there is no current window, otherwise the current index + 1

if N_PARAMS() eq 0 then begin
    index = !D.WINDOW + 1
endif

if (index gt MAXWINDOWS - 1) then begin
    print, 'Too many windows to manage!'
    return
endif

; If window structure system variable does not exist, create it.

DEFSYSV, '!WIN_STRUCT', EXISTS = wstruct_exists

if (wstruct_exists eq 0) then begin

    temp = {windowstruct, xaxis:!X, yaxis:!Y, zaxis:!Z, plotv:!P}
    DEFSYSV, '!WIN_STRUCT', REPLICATE({windowstruct}, MAXWINDOWS)

endif

; If a window already exists, store its x and y axis ranges in the
; window structure

if (!D.WINDOW ne -1) then begin

    !WIN_STRUCT[!D.WINDOW].xaxis = !X
    !WIN_STRUCT[!D.WINDOW].yaxis = !Y
    !WIN_STRUCT[!D.WINDOW].zaxis = !Z
    !WIN_STRUCT[!D.WINDOW].plotv = !P

endif

; Make new window, force IDL to provide Backing Store for the window

device, retain=2 
window, index

end
