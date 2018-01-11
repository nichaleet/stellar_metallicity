;+
; NAME:
;   dirnames
;
; PURPOSE:
;   Separate the subdirectories of a path.
;   
; CALLING SEQUENCE:
;   Result = dirnames(path)
; 
; INPUTS:
;   path       path to be separated
;
; OUTPUTS:
;   Result     array of subdirectory names; the last element in the
;              array is the filename
;
; MODIFICATION HISTORY:
;   06nov28    ENK
;-

function dirnames, path
    return, strsplit(path, path_sep(), /extract)
end     
