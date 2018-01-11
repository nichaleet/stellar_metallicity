function helio2gsr, ra, dec, vhelio
    euler, ra, dec, l, b, 1
    v_lsr = 220.0
    vgsr = vhelio + 10.0*cos(l*!DTOR)*cos(b*!DTOR) + 7.2*sin(b*!DTOR) + (v_lsr + 5.2)*(sin(l*!DTOR)*cos(b*!DTOR))
    return, vgsr
end
