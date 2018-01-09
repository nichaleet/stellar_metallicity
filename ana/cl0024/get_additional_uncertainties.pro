pro get_additional_uncertainties,feh,sn,deltaage,deltafeh
        ;get the uncertainties files
        delstr = mrdfits('/scr2/nichal/workspace2/ana/age_metal_degeneracy/uncertainty_fn_feh_sn.fits',1)
        nobjs = n_Elements(feh)
        deltaage = fltarr(nobjs);in dex
        deltafeh = fltarr(nobjs)
        for i=0,nobjs-1 do begin
           boobee = min(delstr.sn-sn(i),loc)
           deltaage(i) = interpol(delstr(loc).sigmaage,delstr(loc).feh,feh(i))
           deltafeh(i) = interpol(delstr(loc).sigmafeh,delstr(loc).feh,feh(i))
        endfor
end
~        
