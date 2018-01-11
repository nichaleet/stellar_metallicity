function transform_r09, ubvri, ubvrierr, error=cfhterr
    ;UBVRI -> ugriz (CFHT)
    ;Regnault et al. (2009)

    u = ubvri[0]
    b = ubvri[1]
    v = ubvri[2]
    r = ubvri[3]
    i = ubvri[4]
    uerr = ubvrierr[0]
    gerr = ubvrierr[1]
    rerr = ubvrierr[2]
    ierr = ubvrierr[3]
    zerr = ubvrierr[4]
    cfht = dblarr(5)
    cfhterr = dblarr(5)
    cfht[0] = u-b lt 0.4  ? u - 0.2450*(u-b) : v - 0.2450*0.4  - 0.2787*(u-b-0.4)
    cfht[1] = b-v lt 0.45 ? v + 0.4957*(b-v) : v + 0.4957*0.45 + 0.4583*(b-v-0.45)
    cfht[2] = v-r lt 0.65 ? r + 0.1654*(v-r) : r + 0.1654*0.65 + 0.2079*(v-r-0.65)
    cfht[3] = r-i lt 0.40 ? i + 0.2069*(r-i) : i + 0.2069*0.40 + 0.1702*(r-i-0.40)
    cfht[4] = r-i lt 0.35 ? i - 0.1949*(r-i) : i - 0.1949*0.35 - 0.4420*(r-i-0.35)
    cfhterr[0] = u-b lt 0.45 ? sqrt((0.7550*uerr)^2. + (0.2450*berr)^2. + (0.0108*(u-b))^2.) : sqrt((0.7213*uerr)^2. + (0.2787*berr)^2. + (0.0070*(u-b))^2.)
    cfhterr[1] = b-v lt 0.45 ? sqrt((0.4957*berr)^2. + (0.5043*verr)^2. + (0.0153*(b-v))^2.) : sqrt((0.4583*berr)^2. + (0.5417*verr)^2. + (0.0026*(b-v))^2.)
    cfhterr[2] = v-r lt 0.65 ? sqrt((0.1654*verr)^2. + (0.8346*rerr)^2. + (0.0049*(v-r))^2.) : sqrt((0.2079*verr)^2. + (0.7921*rerr)^2. + (0.0248*(v-r))^2.)
    cfhterr[3] = r-i lt 0.40 ? sqrt((0.2069*rerr)^2. + (0.7931*ierr)^2. + (0.0093*(r-i))^2.) : sqrt((0.1702*rerr)^2. + (0.8298*ierr)^2. + (0.0056*(r-i))^2.)
    cfhterr[4] = r-i lt 0.35 ? sqrt((0.1949*rerr)^2. + (1.1949*ierr)^2. + (0.0301*(r-i))^2.) : sqrt((0.4420*rerr)^2. + (1.4420*ierr)^2. + (0.0133*(r-i))^2.)
    return, cfht
end
