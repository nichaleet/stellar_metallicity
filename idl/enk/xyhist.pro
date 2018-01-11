; NAME:
;     xyhist
; PURPOSE:
;     Create a (contour) plot and separate histograms for each axis.
; CALLING SEQUENCE
;      xyhist, x, y, xtitle=xtitle, ytitle=ytitle, xrange=xrange,
;      yrange=yrange, psym=psym, symsize=symsize, xbins=xbins,
;      ybins=ybins, xhistpercent=xhistpercent,
;      yhistpercent=yhistpercent, /noxhist, /noyhist,
;      /xhistlog, /yhistlog, /contour
; INPUTS:
;     x - vector of x data
;     y - vector of y data, same length as x
; OPTIONAL INPUTS:
;     xtitle - title of x axis
;     ytitle - title of y axis
;     xrange - force x range; default is min(x) to max(x)
;     yrange - force y range; default is min(y) to max(y)
;     psym - symbol to plot; default is 4 (diamond)
;     symsize - symbol size
;     xbins - number of bins in x axis histogram; default is 10
;     ybins - number of bins in y axis histogram; default is 10
;     xhistpercent - fraction of plot which x axis histogram
;                    consumes; default is 20%
;     yhistpercent - fraction of plot which y axis histogram
;                    consumes; default is 20%
;     noxhist - suppress x axis histogram
;     noyhist - suppress y axis histogram
;     xhistlog - scale x axis histogram logarithmically
;     yhistlog - scale y axis histogram logarithmically
;     contour - make a contour plot instead of plotting symbols; with
;               this keyword, psym and symsize keywords do nothing
; OUTPUT:
;     plot output to the current graphics device
; REVISION HISTORY:
;     2006-09-05 E. Kirby (ekirby@ucolick.org)

pro xyhist, xin, yin, xerr, yerr, xrange=xrange, yrange=yrange, psym=psym, xbins=xbins, ybins=ybins, xhistpercent=xhistpercent, yhistpercent=yhistpercent, noxhist=noxhist, noyhist=noyhist, xhistlog=xhistlog, yhistlog=yhistlog, contour=contour, _extra=extra_key, x2=x2, y2=y2, title=title
    x = xin
    y = yin
    if keyword_set(xtitle) eq 0 then xtitle = ''
    if keyword_set(ytitle) eq 0 then ytitle = ''
    if keyword_set(xrange) eq 0 then xrange = [min(x), max(x)]
    if keyword_set(yrange) eq 0 then yrange = [min(y), max(y)]
    if keyword_set(psym) eq 0 then psym = 4
    if keyword_set(symsize) eq 0 then symsize = !P.SYMSIZE
    if keyword_set(xbins) eq 0 then xbins = 10
    if keyword_set(ybins) eq 0 then ybins = 10
    if keyword_set(xhistpercent) eq 0 then xhistpercent = 0.2
    if keyword_set(yhistpercent) eq 0 then yhistpercent = 0.2
    if keyword_set(noxhist) then xhistpercent = 0.
    if keyword_set(noyhist) then yhistpercent = 0.
    topbuf = keyword_set(title) ? 0.075 : 0.02
    botbuf = 0.13
    leftbuf = 0.12
    rightbuf = 0.03
    w = where(x ge xrange[0] and x le xrange[1] and y ge yrange[0] and y le yrange[1])
    x = x[w]
    y = y[w]
    if keyword_set(noxhist) eq 0 then begin
        xhist = histogram(x, nbins=xbins, locations=xloc, min=xrange[0], max=xrange[1])
        if keyword_set(x2) then xhist2 = histogram(x2, nbins=xbins, min=xrange[0], max=xrange[1])
        xbinsize = xloc[1] - xloc[0]
        xloc = xloc + 0.5*xbinsize
        xhist = [0., xhist, 0.]
        xloc = [xloc[0]-xbinsize, xloc, xloc[n_elements(xloc)-1]+xbinsize]
        xhistmax = 1.05*(keyword_set(x2) ? max(xhist) > max(xhist2) : max(xhist))
        xhistmag = floor(alog10(xhistmax))
        case xhistmag of
            0: xhistrange = [0., xhistmax]
            1: xhistrange = [0., ceil(double(xhistmax) / 5)*5]
            else: xhistrange = [0., ceil(double(xhistmax) / 10^(xhistmag-1))*10^(xhistmag-1)]
        endcase
        if keyword_set(xhistlog) then xhistrange[0] = 1.
        yticks = 3
        ;ytickstep = ceil(xhistrange[1]/double(yticks+1))
        ;ytickv = (dindgen(yticks+1)+1)*ytickstep
        if keyword_set(x2) then plot, xloc, xhist2, position=[leftbuf, 1.-topbuf-xhistpercent, 1.-rightbuf-yhistpercent, 1.-topbuf], psym=10, xrange=xrange, yrange=xhistrange, ylog=xhistlog, xtickname=replicate(' ', 30), xstyle=5, ystyle=5, yticks=yticks, yminor=2, title=title, color=fsc_color('red')
        plot, xloc, xhist, position=[leftbuf, 1.-topbuf-xhistpercent, 1.-rightbuf-yhistpercent, 1.-topbuf], psym=10, xrange=xrange, yrange=xhistrange, ylog=xhistlog, xtickname=replicate(' ', 30), xstyle=1, ystyle=1, yticks=yticks, yminor=2, title=title, noerase=keyword_set(x2)
    endif
    if keyword_set(noyhist) eq 0 then begin
        yhist = histogram(y, nbins=ybins, locations=yloc, min=yrange[0], max=yrange[1])
        ybinsize = yloc[1] - yloc[0]
        yloc = yloc + 0.5*ybinsize
        yhist = [0., yhist, 0.]
        yloc = [yloc[0]-ybinsize, yloc, yloc[n_elements(yloc)-1]+ybinsize]
        yhistx = fltarr(2*n_elements(yhist))
        yhistx[2*indgen(n_elements(yhist))] = yhist
        yhistx[(2*indgen(n_elements(yhist)))+1] = yhist
        yhisty = fltarr(2*n_elements(yloc))
        yhisty[2*indgen(n_elements(yloc))] = yloc - 0.5*ybinsize
        yhisty[2*indgen(n_elements(yloc))+1] = yloc + 0.5*ybinsize
        if keyword_set(y2) then begin
            yhist2 = histogram(y2, nbins=ybins, min=yrange[0], max=yrange[1])
            yhistx2 = fltarr(2*n_elements(yhist2))
            yhistx2[2*indgen(n_elements(yhist2))] = yhist2
            yhistx2[(2*indgen(n_elements(yhist2)))+1] = yhist2
        endif
        yhistmax = 1.05*(keyword_set(y2) ? max(yhist) > max(yhist2) : max(yhist))
        yhistmag = floor(alog10(yhistmax))
        case yhistmag of
            0: yhistrange = [0., yhistmax]
            1: yhistrange = [0., ceil(double(yhistmax) / 5)*5]
            else: yhistrange = [0., ceil(double(yhistmax) / 10^(yhistmag-1))*10^(yhistmag-1)]
        endcase
        if keyword_set(yhistlog) then yhistrange[0] = 1.
        xticks = 2
        ;xtickstep = ceil(yhistrange[1]/double(xticks+1))
        ;xtickv = (dindgen(xticks+1)+1)*xtickstep
        if keyword_set(y2) then plot, yhistx2, yhisty, position=[1.-rightbuf-yhistpercent, botbuf, 1.-rightbuf, 1.-topbuf-xhistpercent], psym=0, xrange=yhistrange, yrange=yrange, xlog=yhistlog, ytickname=replicate(' ', 30), xstyle=5, ystyle=5, xticks=xticks, xminor=2, noerase=~keyword_set(noxhist), color=fsc_color('red')
        plot, yhistx, yhisty, position=[1.-rightbuf-yhistpercent, botbuf, 1.-rightbuf, 1.-topbuf-xhistpercent], psym=0, xrange=yhistrange, yrange=yrange, xlog=yhistlog, ytickname=replicate(' ', 30), xstyle=1, ystyle=1, xticks=xticks, xminor=2, noerase=(~keyword_set(noxhist) or keyword_set(y2))
    endif
    if keyword_set(contour) then begin
        bin1 = (xrange[1]-xrange[0])/15.0*0.9999
        bin2 = (yrange[1]-yrange[0])/15.0*0.9999
        dens = hist_2d(x, y, bin1=bin1, bin2=bin2, min1=xrange[0], max1=xrange[1], min2=yrange[0], max2=yrange[1])
        nlevels = 15
        levels = 1.5^(findgen(nlevels)-(nlevels-1))*max(dens)
        levels = [0., levels[where(levels ge 1)]]
        contour, dens, bin1*findgen(16)+xrange[0]+0.5*bin1, bin2*findgen(16)+yrange[0]+0.5*bin2, position=[leftbuf, botbuf, 1.-rightbuf-yhistpercent, 1.-topbuf-xhistpercent], xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, levels=levels
    endif else begin
        if n_params() eq 2 then begin
            plot, x, y, position=[leftbuf, botbuf, 1.-rightbuf-yhistpercent, 1.-topbuf-xhistpercent], xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, psym=psym, xstyle=1, ystyle=1, symsize=symsize, noerase=(~keyword_set(noxhist) and ~keyword_set(noyhist)), _extra=extra_key
        endif else begin
            ploterror, x, y, xerr, yerr, position=[leftbuf, botbuf, 1.-rightbuf-yhistpercent, 1.-topbuf-xhistpercent], xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, psym=psym, xstyle=1, ystyle=1, symsize=symsize, /nohat, noerase=(~keyword_set(noxhist) and ~keyword_set(noyhist)), _extra=extra_key
        endelse
    endelse
end
