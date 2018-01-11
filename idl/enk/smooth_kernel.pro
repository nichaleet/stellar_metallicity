function smooth_kernel, x, y1, dy
    n = n_elements(x)
    y2 = y1*0d
    dx = -1*ts_diff(x, 1)
    dx = dx[where(dx gt 0)]
    halfn = floor(double(n)/2.)
    halfwindow = ceil(1.1*4.*max(dy)/min(dx)) < halfn

    gauss = dblarr(n)
    gaussx = x - mean(x)
    for j=0,n_elements(dy)-1 do gauss += exp(-1.0*gaussx^2 / (2.0*dy[j]^2)) / dy[j]

    for i=0L,n-1 do begin
        low = (i-halfwindow) > 0
        high = ((i+halfwindow)-1 < (n-1)) > 0
        lowg = halfn-halfwindow-(i-halfwindow-low)
        highg = halfn+halfwindow-(i+halfwindow-high)
        if low lt n and low lt high then begin
            w = where(abs(x[low:high] - x[i]) lt max(dy)*4., cw) + low
            if cw gt 0 then y2[i] = total(y1[low:high]*gauss[lowg:highg]) / total(gauss[lowg:highg])
        endif
    endfor
    return, y2
end
