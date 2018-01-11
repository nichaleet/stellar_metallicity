;+
; NAME: 
;        LISTMATCH
;
; PURPOSE: 
;   find the indices where a matches b
;   works only for SETS OF UNIQUE INTEGERS, e.g., array indices
;
;   for example: suppose you have a list of people with their ages and
;   social security numbers (AGE, AGE_SS), and a partially overlapping
;   list of people with their incomes and s.s. numbers (INCOME,
;   INCOME_SS).  And you want to correlate ages with incomes in the
;   overlapping subset.  Call
;      LISTMATCH, AGE_SS, INCOME_SS, AGE_IND, INCOME_IND
;   then AGE[AGE_IND] and INCOME[INCOME_IND] will be the desired
;   pair of variables.
;
; AUTHOR:
;   Evan Kirby
;
; CALLING SEQUENCE:
;   LISTMATCH, a, b, a_ind, b_ind
;   
; INPUTS:
;   a and b are sets of unique integers (no duplicate elements)
; 
; OUTPUTS:
;   a_ind, b_ind are the indices indicating which elements of a and b
;   are in common
;
; EXAMPLE:
;   a = [2,4,6,8]
;   b = [6,1,3,2]
;   listmatch, a, b, a_ind, b_ind
;    print, a[a_ind]
;       2       6
;    print, b[b_ind]
;       2       6
;   
;
; MODIFICATION HISTORY:
;   enk 05dec12
;-

pro listmatch, a, b, a_ind, b_ind
    na = n_elements(a)
    nb = n_elements(b)
    a_ind = -1
    b_ind = -1
    for i=0, na-1 do begin
        w = where(b eq a[i], c)
        if c gt 0 then begin
            a_ind = [a_ind, i]
            b_ind = [b_ind, w[0]]
        endif
    endfor
    na_ind = n_elements(a_ind)
    nb_ind = n_elements(b_ind)
    if na_ind gt 1 then a_ind = a_ind[1:na_ind-1]
    if nb_ind gt 1 then b_ind = b_ind[1:nb_ind-1]    
end
