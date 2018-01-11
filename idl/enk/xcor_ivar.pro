function xcor_ivar, x, y, xivar=xivar, yivar=yivar
    n = n_elements(x)
    xmean = mean(x)
    ymean = mean(y)
    if keyword_set(xivar) then begin
        if keyword_set(yivar) then begin
            xcor = total((x-xmean)*xivar*(y-ymean)*yivar) / $
                   (sqrt(total((x-xmean)^2.0*xivar)*total(xivar)* $
                         total((y-ymean)^2.0*yivar)*total(yivar))/n)
        endif else begin
            xcor = total((x-xmean)*xivar*(y-ymean)) / $
                   sqrt(total((x-xmean)^2.0*xivar)*total(xivar)* $
                        total((y-ymean)^2.0)/n)
        endelse
    endif else begin
        if keyword_set(yivar) then begin
            xcor = total((x-xmean)*(y-ymean)*yivar) / $
                   sqrt(total((x-xmean)^2.0)* $
                        total((y-ymean)^2.0*yivar)*total(yivar)/n)
        endif else begin
            xcor = total((x-xmean)*(y-ymean)) / $
                   sqrt(total((x-xmean)^2.0)*total((y-ymean)^2.0))
        endelse
    endelse
    return, xcor
end
