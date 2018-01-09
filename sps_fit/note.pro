pro sps_fit::getscience, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common mask_in, mask_in
    case 1 of
        strmid(mask_in, 0, 8) eq 'cluster1' or strmid(mask_in, 0, 3) eq 'rse': begin
            photfile = '/scr2/nichal/workspace2/SDSSdata/cluster1/spec/cluster1_phot.fits'  
        end
        strmid(mask_in, 0, 6) eq 'Cl1604': begin
            photfile = getenv('UCI')+'Cl1604/catalogs/Cl1604.fits.gz'
        end
        else: message, mask_in+' is not associated with a cluster that I know.'
    endcase
    phot = mrdfits(photfile, 1, /silent)

    common npixcom, npix
    npix = 3800
    if strlowcase(mask_in) eq 'cl1604deimos' then npix = 10625
    if strlowcase(mask_in) eq 'cl1604lris' then npix = 2570

    observatory, 'apo', obs
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        c = n_elements(files)
        plates = strarr(c)
        MJDs = strarr(c)
        fibers = strarr(c)
        for i=0,c-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '-.', /extract)
            plates[i] = extensions[1]
            MJDs[i] = extensions[2]
            fibers[i] = extensions[3]
        endfor
        nspec = n_elements(fibers)
        self.nspec = nspec
        speclist = plates+' '+MJDs+' '+fibers
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
        for i=0,nspec-1 do begin
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.mask = mask_in
            science.fiber = fibers[i]
            science.plate = plates[i]
            science.MJD = MJDs[i]

            data1 = mrdfits(files[i], 1,/silent)
            hdr   = mrdfits(files[i], 2,/silent)
               
            ras = hdr.plug_ra
            decs = hdr.plug_dec
            science.ra  = ras
            science.dec = decs
            science.jdobs = science.MJD + 2400000.5
            if tag_exist(hdr,'fluxobjid') then fluxobjid = hdr.fluxobjid else $  ;for sdss data
               if tag_exist(hdr,'programname') then fluxobjid = hdr.objid else $ ;for boss data
                  fluxobjid='' 
            science.objname = fluxobjid

           ;fix number of elements to match npix
            nlamb = n_elements(data1.loglam)
            if nlamb gt npix then data1=data1[0:npix-1]
            if nlamb lt npix then begin
               sdss_nan = data1[0]
               for ntag=0,n_elements(tag_names(data1))-1 do sdss_nan.(ntag) = 1./0.
               nanarr   = replicate(sdss_nan,npix-nlamb)
               data1 = [data1,nanarr]
            endif
            lambda    = 10.^data1.loglam
            spec      = data1.flux
            ivar      = data1.ivar
            
            ;find a match between phot array and the current science
            ;spherematch, phot.ra, phot.dec, science.ra, science.dec, 1.0/3600., w1, w2
            w1 = where(fluxobjid eq phot.objid,cmatch)
            if cmatch gt 0 then gcirc,2,phot.ra[w1],phot.dec[w1],science.ra,science.dec,dis
            if cmatch gt 0 and dis[0] gt 2 then stop ;there is a problem with matching phot and spec                 
            if w1[0] eq -1 then begin
            endif else begin
                w1=w1[0]
                science.zspec = hdr.z
                science.zsource = hdr.z
                science.u = phot.u[w1]
                science.g = phot.g[w1]
                science.r = phot.r[w1]
                science.i = phot.i[w1]
                science.z = phot.z[w1]
                science.uerr = phot.err_u[w1]
                science.gerr = phot.err_g[w1]
                science.rerr = phot.err_r[w1]
                science.ierr = phot.err_i[w1]
                science.zerr = phot.err_z[w1]
                science.vdisp_sdss = hdr.vdisp
                science.vdisperr_sdss = hdr.vdisp_err
                ;NOTE HERE SHOULD I USE SDSS DERED magnitude instead?
             endelse

            case 1 of
                science.u gt 10.0 and science.u lt 50.0 and science.g gt 10.0 and science.g lt 50.0: science.phot_color = 'ug'
                science.u gt 10.0 and science.u lt 50.0 and science.r gt 10.0 and science.r lt 50.0: science.phot_color = 'ur'
                science.g gt 10.0 and science.g lt 50.0 and science.i gt 10.0 and science.i lt 50.0: science.phot_color = 'gi'
                science.g gt 10.0 and science.g lt 50.0 and science.z gt 10.0 and science.z lt 50.0: science.phot_color = 'gz'
                else: science.phot_color = 'ug'
            endcase

            science.spec1dfile = files[i]
            science.good = 1
            science.goodsky = 1
            goodspec_tag = where(data1.and_mask eq 0,cgoodspec)
            if cgoodspec ne 0 then science.fitmask(goodspec_tag) = 1
            science.fitmask(where(~finite(science.spec))) = 0

            w = where(science.age le 0, c)
            if c gt 0 then begin
                science[w].age = -999d
                science[w].ageerr = -999d
            endif
            science.feh = -999d
            science.feherr = -999d
            science.vdisp = -999d
            science.vdisperr = -999d
            science.alphafe = 0.0

            self.i = i

            w = where(ivar gt 0 and finite(ivar), civar)
            if civar gt 10 then if w[civar-1] ne n_elements(ivar)-1 then ivar[w[civar-11:civar-1]] = 0
            if self.lowsn eq 1 then begin
                common random, seed
                factor = 3.0
                spec[w] += factor*ivar[w]^(-0.5)*randomn(seed, civar, /double, /normal)
                ivar[w] /= (1.0 + factor^2.)
            endif
            n = n_elements(lambda)

            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            science.lambda = lambda
            science.spec = spec
            science.ivar = ivar
            science.skyspec = data1.sky
            science.sdssmodel = data1.model
;            science.airmass = airmass

            self->specres, science
            self->continuum, science
            if array_equal(science.continuum, replicate(-999, n_elements(science.continuum))) then begin
                message, 'Failed at finding continuum for '+strtrim(objnames[i], 2)+'.  Omitting.', /info
                wgood[i] = 0
                continue
            endif
            self->mask, science
            science.spscont = 1.0
            self->indices, science, /noredraw

            self->statusbox, science=science
            scienceall[i] = science
        endfor
        self.i = 0
        wgood = where(wgood eq 1, cgood)
        scienceall = scienceall[wgood]
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        self.nspec = cgood
        self->writescience
        self->specres_mask, self.directory
        objlist = scienceall.objname
        speclist = mask_in+' '+strtrim(string(objlist[wgood]), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)

        self.nspec = n_elements(scienceall)
        speclist = mask_in+' '+strtrim(string(scienceall.objname), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
    endelse
end