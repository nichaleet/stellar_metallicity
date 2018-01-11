;+
; NAME:
;       smooth_clip
;
; PURPOSE:
;       Return an array smoothed with sigma clipping.
;
; CALLING SEQUENCE:
;       Result = smooth_clip(array, window, clip)
;
; INPUTS:
;       array:  array to be smoothed
;       window: size of boxcar smoothing window (should be even)
;       clip:   clipping threshold in sigma
;
; OUTPUTS:
;       Result: the smoothed array
;
; MODIFICATION HISTORY:
;       Mon Apr 23 18:16:54 2007, Evan Kirby
;       
;-

function smooth_clip, x, y, window, clip, niter
    if ~arg_present(window) then window = 100
    if ~arg_present(clip) then clip = 2
    if ~arg_present(niter) then niter = 3
    n = n_elements(y)
    if n_elements(x) ne n then message, 'x and y must have the same size.'
    if window ge n then message, 'The smoothing window size is larger than the input array.'
    s = dblarr(n)
    w2 = round(window/2.0)
    for i=0L,n-1 do begin 
        ;in = lindgen((window < (i+w2)) < (n-1-i+w2)) + (i gt w2 ? i-w2 : 0)
        in = where(x ge x[i]-w2 and x le x[i]+w2, cin)
        cur = y[in]
        w = lindgen(cin)
        for j=0,niter-1 do begin
            imean = (moment(cur[w], sdev=sd))[0]
            ww = where(abs(cur[w]-imean) lt clip*sd, cww)
            if cww le 1 then continue else w = w[ww]
        endfor
        s[i] = mean(cur[w])
    endfor
    return, s
end
