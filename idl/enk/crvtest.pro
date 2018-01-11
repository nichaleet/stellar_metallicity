pro crvtest
    y = (dindgen(101) - 50) * 0.02d
    x = sqrt(1d - y^2d)
    ;print, '01 ', enk_crvlength(x, y)
    ;print, '02 ', enk_crvlength(x[0:1], y[0:1])
    ;print, '03 ', enk_crvlength(x[0:2], y[0:2])
    ;print, '04 ', enk_crvlength(x[0:50], y[0:50])
    ;print, '05 ', enk_crvlength(x[49:50], y[49:50])
    ;print, '06 ', enk_crvlength(x[0:51], y[0:51])
    ;print, '07 ', enk_crvlength(x[99:100], y[99:100])
    ;print, '08 ', enk_crvlength(x[49:51], y[49:51])
    ;print, '09 ', enk_crvlength(x[0:4], y[0:4])
    ;print, '10 ', enk_crvlength(x[48:52], y[48:52])
    ;print, '11 ', enk_crvlength(x[0:6], y[0:6])
    ;print, '12 ', enk_crvlength(x[47:53], y[47:53])

    y = (dindgen(1001) - 500) * 0.02d
    x = sin(y)
    plot, x, y
    print, enk_crvlength(x, y)
end
