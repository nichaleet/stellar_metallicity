;+
; NAME:
;
;  MAKELIKE2003
;
; PURPOSE:
;
;  Take a full zcat (or STRUCTINDEX results) and return a zcat
;  equivalent to if the standard DEEP2 target selection was used in
;  all fields.
;
; CATEGORY:
;
;  Databases
;
; CALLING SEQUENCE:
;
;  result=makelike2003(zcat,[SEED=seed,/KEEPBAD,/SILENT,BADLIST=badlist])
;
; INPUTS:
;
;  zcat -- a zcat structure or structure from structindex
;
; KEYWORD PARAMETERS:
;
;  SEED -- random number seed to use in selecting the subset of
;          objects to keep.  Set this to some specific value in your
;          code if you want to ensure
;          always getting the same set of objects.
; 
;  /KEEPBAD -- if set, objects from 2002 masks will always be included
;             in the result structure (i.e., we will ignore the fact they
;             differ from 2003 selection).  For situations where we
;             don't care about spatial information and want a clean
;             dataset, the default behavior is best; when spatial
;             uniformity is an issue, we probably have to keep 2002
;             data, so set KEEPBAD then.  
;
;  /SILENT -- if set, the list of masks rejected will be suppressed. 
;
; OUTPUTS:
;
;  result -- as zcat, but we only keep a subset of the original objects:
;     1) Objects on bad masks (see RESTRICTIONS) are thrown out
;
;     2) Objects on 2003, non-EGS mask designs are always kept
;
;     3) A subset of all other objects are kept, such that the
;     probability of being placed on a mask for that subset, given an
;     object's properties (colors/mag/size) is the
;     same as the probability it would have been in 2003 maskmaking
;     (ignoring dependence of selection probability on crowding, which
;     we effectively assume is separable and has the same effect on
;     EGS maskmaking and standard maskmaking).  
;
;   The net result is that we have a zcat which can be treated as if
;   all objects were observed on 2003 masks: equivalent redshift
;   histograms, color selection, etc.
;
; OPTIONAL OUTPUTS:
;
;  BADLIST - array of the mask names whose objects were thrown out
;            entirely. 
;
; RESTRICTIONS:
;
;  Currently only does reselection for EGS fields 11-13; 2002 masks
;  are either rejected or kept, depending on the value of KEEPBAD,
;  while field 14 masks are always rejected.  When their selection
;  gets characterized, we can fix that.
;
;  There is a small slice of color-color space (around R-I=0.5,
;  B-R=0.5) where the selection probability in 2003 is _higher_ than
;  in the EGS (typically by ~20%).  These objects will be undersampled
;  compared to a 2003 dataset by that amount.  They constitute
;  ~1.5% of all EGS objects / ~3% of 2003-selected objects; median
;  redshift for such objects is z=0.75, so the effect on high-z
;  samples is smaller; we are undersampling 2% of the high-z
;  population by ~20%.  I very much doubt this should be of
;  concern to anyone.
;
; PROCEDURE:
;
;   result=makelike2003(zcat)
;
; MODIFICATION HISTORY:
;   4aug05 jan
;-



function makelike2003,zcat,seed=seed,keepbad=keepbad,silent=silent,badlist=badlist

if total(tag_names(zcat) eq 'PSEL03') gt 0 then selstructs=zcat else $
  selstructs=mrdfits(getenv('DEEP2PRODUCTS')+'/selcat.latest.fits',1,/silent)
if n_elements(keepbad) eq 0 then keepbad = 0


whbad=where(selstructs.p_onmask eq -1.,badct)

badmasks=0
if badct gt 0 then badmasks=selstructs[whbad].maskname

badmasks=badmasks[sort(badmasks)]
badmasks=badmasks[uniq(badmasks)]


; this code was needed to deal with MJD=0 objects.  Put a fix in
; writeselection, so this code probably isn't necessary now.

;badct=long(badmasks)*0
;for i=0,n_elements(badmasks) -1 do begin 
;    wh=where(selstructs.maskname eq badmasks[i] AND selstructs.p_onmask eq -1.,ct)
;    badct[i]=ct
;endfor

badmasks=badmasks[where(badct gt 100 OR badmasks lt 2000)]


if n_elements(silent) eq 0 and keepbad eq 0 then $
  print,'Warning: This routine will throw out the following masks: ',badmasks $
  else if keepbad eq 1 then print,$
     'Warning: This routine will throw out the following masks: ',$
     badmasks[where(badmasks lt 2000)]

print

badlist=badmasks


if keepbad eq 1 then begin
    whbad=where(selstructs.p_onmask lt 0 AND selstructs.maskname ge 2000)
    selstructs[whbad].p_onmask = selstructs[whbad].psel03
endif


; Just in case the zcat and selection info mismatch
if n_elements(zcat) lt n_elements(selstructs) then begin
    print,'Mismatch between zcat and selection catalog - having to rematch!'
    buffer=20
; assume both are sorted by objno
    sel=replicate(selstructs[0],n_elements(zcat))
    sel.p_onmask = -1.
    sel.psel03=1.
    idx=0
    nz=n_elements(zcat)-1
    for i=0L,n_elements(zcat)-1 do begin
        wh=where(zcat[i].objno eq $
                 selstructs[(idx-buffer)>0:(idx+buffer)<nz].objno AND $
                 zcat[i].maskname eq $
                 selstructs[(idx-buffer)>0:(idx+buffer)<nz].maskname AND $
                 zcat[i].mjd eq $
                 selstructs[(idx-buffer)>0:(idx+buffer)<nz].mjd,ct)
        idx=((idx-buffer)>0)+wh[0]
        if ct eq 0 then begin 
            print,'failed to find in search region!'
            wh=where(zcat[i].objno eq $
                 selstructs.objno AND $
                 zcat[i].maskname eq $
                 selstructs.maskname AND $
                 zcat[i].mjd eq $
                 selstructs.mjd,ct)
            idx=wh[0]
         endif

     sel[i]=selstructs[idx]

    endfor


endif else if n_elements(zcat) lt n_elements(selstructs) then $
     message,'Selection structures have fewer elements than zcat - crashing!!!'
 
if n_elements(sel) eq 0 then sel=selstructs

keep=where(randomu(seed,n_elements(zcat)) lt sel.psel03/sel.p_onmask)

return,zcat[keep]
end
