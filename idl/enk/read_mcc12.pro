function read_mcc12
    woo08 = {name:' ', logmstar:0d, ml:0d, logl:0d, logz:0d, feh:0d}
    woo08 = replicate(woo08, 41)
    woo08.name = ['M33', 'LMC', 'NGC 55', 'SMC', 'M32', 'NGC 205', 'IC 10', 'NGC 6822', 'NGC 147', 'NGC 185', 'NGC 3109', 'IC 1613', 'WLM', 'IC 5152', 'Sextans B', 'Sextans A', 'Fornax', 'EGB0427+63', 'Andromeda I', 'Pegasus dIrr', 'Andromeda VII', 'Leo I', 'Andromeda II', 'UKS 2323-326', 'GR 8', 'Antlia', 'Sagittarius dIrr', 'Andromeda III', 'Leo A', 'Aquarius', 'Andromeda VI', 'Leo II', 'Phoenix', 'Sculptor', 'Tucana', 'Draco', 'Sextans (I)', 'LGS 3', 'Ursa Minor', 'Carina', 'Andromeda V']
    woo08.logmstar = [9.55, 9.19, 9.04, 8.67, 8.66, 8.65, 8.30, 8.23, 8.17, 8.14, 8.13, 8.03, 7.65, 7.54, 7.54, 7.43, 7.27, 6.96, 6.85, 6.83, 6.67, 6.66, 6.63, 6.51, 6.42, 6.41, 6.36, 6.27, 6.26, 6.25, 6.17, 6.16, 6.10, 6.08, 5.97, 5.96, 5.93, 5.86, 5.75, 5.71, 5.60]
    woo08.ml = [1.2, 0.7, 0.8, 0.8, 1.3, 1.4, 0.9, 0.8, 1.6, 1.0, 0.8, 1.0, 0.9, 0.5, 0.8, 0.5, 1.2, 1.0, 1.6, 1.0, 0.9, 0.9, 1.0, 0.6, 0.7, 1.0, 0.4, 1.8, 0.5, 0.9, 0.5, 1.6, 1.8, 1.7, 1.6, 1.8, 1.6, 1.0, 1.9, 1.0, 1.1]
    woo08.logz = [-2.7, -2.3, -999, -2.9, -2.8, -2.2, -3.0, -2.9, -2.8, -2.5, -3.4, -3.1, -3.1, -3.1, -3.8, -3.6, -2.9, -999, -3.1, -3.7, -3.2, -3.1, -3.2, -999, -3.7, -3.6, -4.0, -3.4, -3.8, -3.6, -3.4, -3.3, -3.6, -3.2, -3.4, -3.7, -3.6, -3.4, -3.6, -3.5, -3.6]
    woo08.logl = alog10(10.^(woo08.logmstar)/woo08.ml)
    woo08.feh = woo08.logz - alog10(0.019)

    mar08 = {name:' ', lum:0d, lumerrupper:0d, lumerrlower:0d, mstar:0d, mstarerrupper:0d, mstarerrlower:0d}
    mar08 = replicate(mar08, 12)
    mar08.name = ['Bootes (I)', 'Bootes II', 'Canes Venatici (I)', 'Canes Venatici II', 'Coma Berenices', 'Hercules', 'Leo IV', 'Segue (I)', 'Ursa Major (I)', 'Ursa Major II', 'Willman 1', 'Draco']
    mar08.lum = [3.0d4, 1.0d3, 2.3d5, 7.9d3, 3.7d3, 3.6d4, 8.7d3, 335, 1.4d4, 4.0d3, 1000, 2.7d5]
    mar08.lumerrupper = [0.6d4, 0.8d3, 0.3d5, 3.4d3, 1.8d3, 1.1d4, 4.4d3, 235, 0.4d4, 1.8d3, 660, 0.4d5]
    mar08.lumerrlower = [0.6d4, 0.8d3, 0.3d5, 3.7d3, 1.6d3, 1.1d4, 4.7d3, 185, 0.4d4, 1.9d3, 670, 0.4d5]
    mar08.mstar = [3.4d4, 1.4d3, 3.0d5, 8.0d3, 4.8d3, 3.7d4, 8.5d3, 600, 1.9d4, 5.4d3, 1.5d3, 3.2d5]
    mar08.mstarerrupper = [0.3d4, 0.7d3, 0.2d5, 2.0d3, 0.9d3, 0.6d4, 3.0d3, 115, 0.3d4, 0.6d3, 0.3d3, 0.1d5]
    mar08.mstarerrlower = [0.3d4, 0.5d3, 0.2d5, 1.7d3, 0.9d3, 0.6d4, 2.1d3, 105, 0.3d4, 0.5d3, 0.3d3, 0.1d5]

    mcc12 = {name:' ', abbrev:' ', system:' ', $
             ra:0d, dec:0d, $
             ebv:0d, $
             dm:0d, dmerr:0d, dmerrupper:0d, dmerrlower:0d, $
             vhelio:0d, vhelioerr:0d, vhelioerrupper:0d, vhelioerrlower:0d, $
             vmag:0d, vmagerr:0d, vmagerrupper:0d, vmagerrlower:0d, $
             mv:0d, mverr:0d, mverrupper:0d, mverrlower:0d, $
             lum:0d, lumerr:0d, lumerrupper:0d, lumerrlower:0d, $
             logL:0d, logLerr:0d, logLerrupper:0d, logLerrlower:0d, $
             mstar:0d, mstarerr:0d, mstarerrupper:0d, mstarerrlower:0d, $
             logmstar:0d, logmstarerr:0d, logmstarerrupper:0d, logmstarerrlower:0d, $
             pa:0d, paerr:0d, paerrupper:0d, paerrlower:0d, $
             ellip:0d, elliperr:0d, elliperrupper:0d, elliperrlower:0d, $
             muv:0d, muverr:0d, muverrupper:0d, muverrlower:0d, $
             rhalf:0d, rhalferr:0d, rhalferrupper:0d, rhalferrlower:0d, $
             sigmav:0d, sigmaverr:0d, sigmaverrupper:0d, sigmaverrlower:0d, $
             mhalf:0d, mhalferr:0d, mhalferrupper:0d, mhalferrlower:0d, $
             ml:0d, mlerr:0d, mlerrupper:0d, mlerrlower:0d, $
             mlstar:1.6d, $
             vrot:0d, vroterr:0d, vroterrupper:0d, vroterrlower:0d, $
             mhi:0d, $
             sigmavhi:0d, sigmavhierr:0d, sigmavhierrupper:0d, sigmavhierrlower:0d, $
             vrothi:0d, vrothierr:0d, vrothierrupper:0d, vrothierrlower:0d, $
             feh:0d, feherr:0d, feherrupper:0d, feherrlower:0d, $
             flag:0, $
             refs:' '}
    mcc12 = replicate(mcc12, 1000)
    name = ' '
    dec1 = ' '
    openr, lun, getenv('M31')+'papers/mcc12/NearbyGalaxies.dat', /get_lun
    skip_lun, lun, 34, /lines
    i = 0L
    while ~eof(lun) do begin
        readf, lun, name, ra1, ra2, ra3, dec1, dec2, dec3, ebv, dm, dmerrupper, dmerrlower, vhelio, vhelioerrupper, vhelioerrlower, vmag, vmagerrupper, vmagerrlower, pa, paerrupper, paerrlower, ellip, elliperrupper, elliperrlower, muv, muverrupper, muverrlower, rhalf, rhalferrupper, rhalferrlower, sigmav, sigmaverrupper, sigmaverrlower, vrot, vroterrupper, vroterrlower, mhi, sigmavhi, sigmavhierrupper, sigmavhierrlower, vrothi, vrothierrupper, vrothierrlower, feh, feherrupper, feherrlower, flag, format='(A19,2(I2,1X),D4,1X,A3,2(1X,I2),2(1X,D5),2(1X,D4),1X,D6,3(1X,D4),2(1X,D3),1X,D5,6(1X,D4),2(1X,D3),1X,D6,2(1X,D5),10(1X,D4),1X,D5,2(1X,D4),1X,D5,2(1X,D4),1X,I1)'

        mcc12[i].name = strtrim(name, 2)
        mcc12[i].ra = (ra1 + ra2 / 60d + ra3 / 3600d) * 15d
        decsign = stregex(dec1, '-', /boolean) ? -1 : 1
        mcc12[i].dec = decsign*(abs(double(dec1)) + dec2 / 60d + dec3 / 3600d)
        mcc12[i].ebv = ebv 
        mcc12[i].dm = dm
        mcc12[i].dmerr = (dmerrupper + dmerrlower) / 2.
        mcc12[i].dmerrupper = dmerrupper
        mcc12[i].dmerrlower = dmerrlower
        mcc12[i].vhelio = vhelio
        mcc12[i].vhelioerr = (vhelioerrupper + vhelioerrlower) / 2.
        mcc12[i].vhelioerrupper = vhelioerrupper
        mcc12[i].vhelioerrlower = vhelioerrlower
        mcc12[i].vmag = vmag
        mcc12[i].vmagerr = (vmagerrupper + vmagerrlower) / 2.
        mcc12[i].vmagerrupper = vmagerrupper
        mcc12[i].vmagerrlower = vmagerrlower
        mcc12[i].mv = vmag - dm

        mcc12[i].pa = pa
        mcc12[i].paerr = (paerrupper + paerrlower) / 2.
        mcc12[i].paerrupper = paerrupper
        mcc12[i].paerrlower = paerrlower
        mcc12[i].ellip = ellip
        mcc12[i].elliperr = (elliperrupper + elliperrlower) / 2.
        mcc12[i].elliperrupper = elliperrupper
        mcc12[i].elliperrlower = elliperrlower
        mcc12[i].muv = muv
        mcc12[i].muverr = (muverrupper + muverrlower) / 2.
        mcc12[i].muverrupper = muverrupper
        mcc12[i].muverrlower = muverrlower
        mcc12[i].rhalf = rhalf
        mcc12[i].rhalferr = (rhalferrupper + rhalferrlower) / 2.
        mcc12[i].rhalferrupper = rhalferrupper
        mcc12[i].rhalferrlower = rhalferrlower
        mcc12[i].sigmav = sigmav
        mcc12[i].sigmaverr = (sigmaverrupper + sigmaverrlower) / 2.
        mcc12[i].sigmaverrupper = sigmaverrupper
        mcc12[i].sigmaverrlower = sigmaverrlower
        mcc12[i].vrot = vrot
        mcc12[i].vroterr = (vroterrupper + vroterrlower) / 2.
        mcc12[i].vroterrupper = vroterrupper
        mcc12[i].vroterrlower = vroterrlower
        mcc12[i].mhi = mhi
        mcc12[i].sigmavhi = sigmavhi
        mcc12[i].sigmavhierr = (sigmavhierrupper + sigmavhierrlower) / 2.
        mcc12[i].sigmavhierrupper = sigmavhierrupper
        mcc12[i].sigmavhierrlower = sigmavhierrlower
        mcc12[i].vrothi = vrothi
        mcc12[i].vrothierr = (vrothierrupper + vrothierrlower) / 2.
        mcc12[i].vrothierrupper = vrothierrupper
        mcc12[i].vrothierrlower = vrothierrlower
        mcc12[i].feh = feh
        mcc12[i].feherr = (feherrupper + feherrlower) / 2.
        mcc12[i].feherrupper = feherrupper
        mcc12[i].feherrlower = feherrlower
        mcc12[i].flag = flag
        i++
    endwhile
    close, lun
    free_lun, lun
    mcc12 = mcc12[0:i-1]

    mcc12.mverr = sqrt(mcc12.vmagerr^2. + mcc12.dmerr^2.)
    mcc12.mverrupper = sqrt(mcc12.vmagerrupper^2. + mcc12.dmerrupper^2.)
    mcc12.mverrlower = sqrt(mcc12.vmagerrlower^2. + mcc12.dmerrlower^2.)
    mcc12.lum = 10.^((4.83 - mcc12.mv)/2.5)
    mcc12.lumerr = mcc12.mverr*alog(10.)/2.5*mcc12.lum
    mcc12.lumerrupper = 10.^((4.83 - mcc12.mv + mcc12.mverrlower)/2.5) - mcc12.lum
    mcc12.lumerrlower = mcc12.lum - 10.^((4.83 - mcc12.mv - mcc12.mverrupper)/2.5)
    mcc12.logL = (4.83 - mcc12.mv) / 2.5
    mcc12.logLerr = mcc12.mverr / 2.5
    mcc12.logLerrupper = mcc12.mverrupper / 2.5
    mcc12.logLerrlower = mcc12.mverrlower / 2.5

    match, strtrim(mcc12.name, 2), strtrim(woo08.name, 2), w1, w2
    mcc12[w1].mlstar = woo08[w2].ml
    w = where(strtrim(mcc12.name, 2) eq 'UGC 4879' or strtrim(mcc12.name, 2) eq 'Leo T')
    mcc12[w].mlstar = 1.0

    mcc12.mstar = mcc12.lum * mcc12.mlstar
    mcc12.mstarerr = mcc12.mlstar * mcc12.lumerr
    mcc12.mstarerrupper = mcc12.mlstar * mcc12.lumerrupper
    mcc12.mstarerrlower = mcc12.mlstar * mcc12.lumerrlower
    mcc12.logmstar = alog10(mcc12.mstar)
    mcc12.logmstarerr = mcc12.mverr / 2.5
    mcc12.logmstarerrupper = mcc12.mverrupper / 2.5
    mcc12.logmstarerrlower = mcc12.mverrlower / 2.5

    match, mcc12.name, mar08.name, w1, w2
    mcc12[w1].lum = mar08[w2].lum
    mcc12[w1].lumerrupper = mar08[w2].lumerrupper
    mcc12[w1].lumerrlower = mar08[w2].lumerrlower
    mcc12[w1].logL = alog10(mcc12[w1].lum)
    mcc12[w1].logLerrupper = mcc12[w1].lumerrupper / (mcc12[w1].lum * alog(10.))
    mcc12[w1].logLerrlower = mcc12[w1].lumerrlower / (mcc12[w1].lum * alog(10.))
    mcc12[w1].mstar = mar08[w2].mstar
    mcc12[w1].mstarerrupper = mar08[w2].mstarerrupper
    mcc12[w1].mstarerrlower = mar08[w2].mstarerrlower
    mcc12[w1].logmstar = alog10(mcc12[w1].mstar)
    mcc12[w1].logmstarerrupper = mcc12[w1].mstarerrupper / (mcc12[w1].mstar * alog(10.))
    mcc12[w1].logmstarerrlower = mcc12[w1].mstarerrlower / (mcc12[w1].mstar * alog(10.))
    mcc12[w1].mlstar = mcc12[w1].mstar / mcc12[w1].lum

    G = 6.67300d-11
    pc = 3.08568025d16
    Msun = 1.98892d30

    dist = 10. * 10.^(mcc12.dm / 5.)
    disterr = mcc12.dmerr * alog(10.) / 5. * dist
    disterrupper = mcc12.dmerrupper * alog(10.) / 5. * dist
    disterrlower = mcc12.dmerrlower * alog(10.) / 5. * dist
    ;disterrupper = (10. * 10.^((mcc12.dm + mcc12.dmerrupper) / 5.)) - dist
    ;disterrlower = dist - (10. * 10.^((mcc12.dm - mcc12.dmerrlower) / 5.))

    mcc12.mhalf = 4. / G * (dist*tan(mcc12.rhalf / 60. * !DTOR)*pc) * (mcc12.sigmav*1d3)^2 / Msun
    mcc12.mhalferr = mcc12.mhalf * sqrt((disterr/dist)^2. + (mcc12.rhalferr/60.*!DTOR/(sin(mcc12.rhalf / 60. * !DTOR)*cos(mcc12.rhalf / 60. * !DTOR)))^2. + (2.*mcc12.sigmaverr/mcc12.sigmav)^2.)
    mcc12.mhalferrupper = mcc12.mhalf * sqrt((disterrupper/dist)^2. + (mcc12.rhalferrupper/60.*!DTOR/(sin(mcc12.rhalf / 60. * !DTOR)*cos(mcc12.rhalf / 60. * !DTOR)))^2. + (2.*mcc12.sigmaverrupper/mcc12.sigmav)^2.)
    mcc12.mhalferrlower = mcc12.mhalf * sqrt((disterrlower/dist)^2. + (mcc12.rhalferrlower/60.*!DTOR/(sin(mcc12.rhalf / 60. * !DTOR)*cos(mcc12.rhalf / 60. * !DTOR)))^2. + (2.*mcc12.sigmaverrlower/mcc12.sigmav)^2.)
    
    mcc12.ml = mcc12.mhalf / mcc12.lum
    mcc12.mlerr = mcc12.ml * sqrt((mcc12.mhalferr/mcc12.mhalf)^2. + (mcc12.lumerr/mcc12.lum)^2.)
    mcc12.mlerrupper = mcc12.ml * sqrt((mcc12.mhalferrupper/mcc12.mhalf)^2. + (mcc12.lumerrupper/mcc12.lum)^2.)
    mcc12.mlerrlower = mcc12.ml * sqrt((mcc12.mhalferrlower/mcc12.mhalf)^2. + (mcc12.lumerrlower/mcc12.lum)^2.)

    mwnames = ['The Galaxy', 'Canis Major', 'Sagittarius dSph', 'Segue (I)', 'Ursa Major II', 'Bootes II', 'Segue II', 'Willman 1', 'Coma Berenices', 'Bootes III', 'LMC', 'SMC', 'Bootes (I)', 'Draco', 'Ursa Minor', 'Sculptor', 'Sextans (I)', 'Ursa Major (I)', 'Carina', 'Hercules', 'Fornax', 'Leo IV', 'Canes Venatici II', 'Leo V', 'Pisces II', 'Canes Venatici (I)', 'Leo II', 'Leo I']
    abbrev1 = ['MW', 'CMa', 'Sgr', 'Seg1', 'UMaII', 'BooII', 'Seg2', 'W1', 'ComBer', 'BooIII', 'LMC', 'SMC', 'BooI', 'Dra', 'UMi', 'Scl', 'Sex', 'UMaI', 'Car', 'Her', 'For', 'LeoIV', 'CVnII', 'LeoV', 'PscII', 'CVnI', 'LeoII', 'LeoI']
    m31names = ['Andromeda', 'M32', 'Andromeda IX', 'NGC 205', 'Andromeda I', 'Andromeda XVII', 'Andromeda XXVII', 'Andromeda III', 'Andromeda XXV', 'Andromeda XXVI', 'Andromeda V', 'Andromeda XI', 'Andromeda XIX', 'Andromeda XXIII', 'Andromeda XX', 'Andromeda XIII', 'Andromeda X', 'Andromeda XXI', 'And XXXII', 'NGC 147', 'And XXX', 'Andromeda XIV', 'Andromeda XII', 'Andromeda XV', 'Andromeda II', 'NGC 185', 'Andromeda XXIX', 'Triangulum', 'Andromeda XXIV', 'Andromeda VII', 'IC 10', 'And XXXI', 'LGS 3', 'Andromeda VI', 'Andromeda XXII', 'Andromeda XVI', 'Andromeda XXVIII']
    abbrev2 = ['M31', 'M32', 'AndIX', 'N205', 'AndI', 'AndXVII', 'AndXXVII', 'AndIII', 'AndXXV', 'AndXXVI', 'AndV', 'AndXI', 'AndXIX', 'AndXXIII', 'AndXX', 'AndXIII', 'AndX', 'AndXXI', 'AndXXXII', 'N147', 'AndXXX', 'AndXIV', 'AndXII', 'AndXV', 'AndII', 'N185', 'AndXXIX', 'M33', 'AndXXIV', 'AndVII', 'IC10', 'AndXXXI', 'LGS3', 'AndVI', 'AndXXII', 'AndXVI', 'AndXXVIII']
    lgnames = ['IC 1613', 'Phoenix', 'NGC 6822', 'Cetus', 'Pegasus dIrr', 'Leo T', 'WLM', 'Andromeda XVIII', 'Leo A', 'Aquarius', 'Tucana', 'Sagittarius dIrr', 'UGC 4879', 'NGC 3109', 'Sextans B', 'Antlia', 'Sextans A', 'HIZSS 3(A)', 'HIZSS 3B', 'ESO 410- G 005', 'KKR 25', 'NGC 55', 'Leo P', 'ESO 294- G 010', 'NGC 300', 'IC 5152', 'KKH 98', 'UKS 2323-326', 'KKR 3', 'GR 8', 'UGC 9128', 'UGC 8508', 'IC 3104', 'DDO 125', 'UGCA 86', 'DDO 99', 'IC 4662', 'DDO 190', 'KKH 86', 'NGC 4163', 'DDO 113']
    abbrev3 = ['I1613', 'Phe', 'N6822', 'Cet', 'Peg', 'LeoT', 'WLM', 'AndXVIII', 'LeoA', 'Aqr', 'Tuc', 'SagDIG', 'VV124', 'N3109', 'SexB', 'Ant', 'SexA', 'HIZSS3A', 'HISZZ3B', 'ESO410', 'KKR25', 'N55', 'LeoP', 'ESO294', 'N300', 'I5152', 'KKH98', 'UKS2323', 'KKR3', 'GR8', 'U9128', 'U8508', 'I3104', 'DDO125', 'U86', 'DDO99', 'I4662', 'DDO190', 'KKH86', 'N4163', 'DDO113']
    match, strtrim(mcc12.name, 2), strtrim(mwnames, 2), w1, w2
    mcc12[w1].system = 'MW'
    mcc12[w1].abbrev = abbrev1[w2]
    match, strtrim(mcc12.name, 2), strtrim(m31names, 2), w1, w2
    mcc12[w1].system = 'M31'
    mcc12[w1].abbrev = abbrev2[w2]
    match, strtrim(mcc12.name, 2), strtrim(lgnames, 2), w1, w2
    mcc12[w1].system = 'LG'
    mcc12[w1].abbrev = abbrev3[w2]

    return, mcc12
end
