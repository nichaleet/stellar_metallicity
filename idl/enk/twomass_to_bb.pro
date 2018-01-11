;Carpenter (2001)
function twomass_to_bb, jhk_2mass, error2mass=jhkerr_2mass, errorbb=bberr
    j2m = jhk_2mass[0]
    h2m = jhk_2mass[1]
    k2m = jhk_2mass[2]
    j2merr = jhkerr_2mass[0]
    h2merr = jhkerr_2mass[1]
    k2merr = jhkerr_2mass[2]

    j2mk2m = j2m-k2m
    h2mk2m = h2m-k2m
    
    jbbkbb = 1.029*j2mk2m + 0.011
    hbbkbb = 1.004*h2mk2m - 0.028
    kbb = k2m + 0.044
    jbb = jbbkbb + kbb
    hbb = hbbkbb + kbb
    bb = [jbb, hbb, kbb]

    j2mk2merr = sqrt(j2merr^2. + k2merr^2)
    h2mk2merr = sqrt(h2merr^2. + k2merr^2)
    
    jbbkbberr = sqrt((1.029*j2mk2merr)^2.0 + (0.006*j2mk2m)^2.0 + 0.005^2.0)
    hbbkbberr = sqrt((1.004*h2mk2merr)^2.0 + (0.019*h2mk2m)^2.0 + 0.005^2.0)
    kbberr = sqrt(k2merr^2.0 + (0.005*jbbkbb)^2.0 + 0.003^2.0)
    jbberr = sqrt(jbbkbberr^2.0 - kbberr^2.0)
    hbberr = sqrt(hbbkbberr^2.0 - kbberr^2.0)

    bberr = [jbberr, hbberr, kbberr]
    return, bb
end
