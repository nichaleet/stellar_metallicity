pro get_ha_ew_sdss
        sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass2/sps_fit.fits.gz',1,/silent)
	nsci = n_elements(sci)

        ha_ew_arr = fltarr(nsci)
        ha_ew_err_arr = fltarr(nsci)
        for i=0,nsci-1 do begin
                zin = sci(i).zspec
                ha_ew_arr(i)=Ha_ew(sci(i).contdiv,sci(i).lambda,sci(i).contdivivar,zin,widtherr=widtherr)
                ha_ew_err_arr(i) = widtherr
        endfor
        save,ha_ew_arr,ha_ew_err_arr,filename='gallazzi_allmass2_Ha_ew.sav'


end
