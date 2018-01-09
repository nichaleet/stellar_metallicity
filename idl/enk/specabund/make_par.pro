pro make_par, parfile=parfile, linefile=linefile, atmfile=atmfile, outfile=outfile, driver=driver, minlambda=minlambda, maxlambda=maxlambda, c12c13=c12c13, stronglist=stronglist, atomic=atomic, extralines=extralines
    if ~keyword_set(parfile) then parfile = 'star.par'
    if ~keyword_set(atmfile) then atmfile = file_dirname(parfile, /mark_directory)+file_basename(parfile, '.par')+'.atm'
    if ~keyword_set(outfile) then outfile = file_dirname(parfile, /mark_directory)+file_basename(parfile, '.par')+'.out2'
    if ~keyword_set(minlambda) then minlambda = 6300.0
    if ~keyword_set(maxlambda) then maxlambda = 9100.0
    if ~keyword_set(driver) then driver = 'synth'
    if driver ne 'synth' and ~keyword_set(linefile) then linefile = file_dirname(parfile, /mark_directory)+file_basename(parfile, '.par')+'.ew'

    openw, lun2, parfile, /get_lun
    printf, lun2, driver
    printf, lun2, "terminal       'x11'"
    printf, lun2, "standard_out   '"+file_dirname(parfile, /mark_directory)+file_basename(parfile, '.par')+".out1'"
    printf, lun2, "summary_out    '"+outfile+"'"
    printf, lun2, "smoothed_out   '"+file_dirname(parfile, /mark_directory)+file_basename(parfile, '.par')+".out3'"
    printf, lun2, "model_in       '"+atmfile+"'"
    if driver eq 'synth' then begin
        ;printf, lun2, "lines_in       '"+getenv('ahome')+"m31/solspec/linelists/valdnistkirby.moog'"
        if keyword_set(stronglist) then printf, lun2, "stronglines_in '"+stronglist+"'"
        ;printf, lun2, "lines_in       '"+getenv('CALTECH')+"linelist/laboratory.dat'"
        ;printf, lun2, "lines_in       '"+getenv('M31')+"espec/synths/linelist.4391'"
        printf, lun2, "lines_in       '"+linefile+"'"
        ;printf, lun2, "stronglines_in '"+getenv('M31')+"espec/synths/stronglist.4391'"
        printf, lun2, "strong        "+(keyword_set(stronglist) ? '1' : '0')
    endif else begin
        printf, lun2, "lines_in       '"+linefile+"'"
    endelse
    if driver eq 'blends' then begin
        if ~keyword_set(atomic) then message, 'You must specify an atomic number with the blend driver.'
        printf, lun2, "blenlimits"
        printf, lun2, "     0.150    0.005  "+strtrim(atomic, 2)
    endif
    printf, lun2, "atmosphere    1"
    printf, lun2, "molecules     1"
    printf, lun2, "damping       1"
    printf, lun2, "trudamp       0"
    printf, lun2, "lines         1"
    if keyword_set(c12c13) then begin
        printf, lun2, "isotopes      2         1"
        printf, lun2, " 106.00112    "+string(1d + (1d / c12c13), format='(D6.2)')
        printf, lun2, " 106.00113    "+string(1d + c12c13, format='(D6.2)')
    endif
    printf, lun2, "flux/int      0"
    printf, lun2, "plot          0"
    if driver eq 'synth' then begin
        printf, lun2, "synlimits"
        printf, lun2, "  "+string(minlambda, format='(D8.3)')+" "+string(maxlambda, format='(D8.3)')+"  0.02  1.00"
    endif

    if keyword_set(extralines) then begin
        for i=0,n_elements(extralines)-1 do printf, lun2, extralines[i]
    endif

    printf, lun2, "obspectrum    0"
;    printf, lun2, "plotpars      0", format='(A15,$)'

    close, lun2
    free_lun, lun2    
end
