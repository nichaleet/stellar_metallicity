pro separate_hires, imlist
    n = n_elements(imlist)
    for i=0,n-1 do begin
        hdr0 = headfits('hires'+string(imlist[i], format='(I04)')+'.fits')
        for j=1,3 do begin
            a = mrdfits('hires'+string(imlist[i], format='(I04)')+'.fits', j, hdr, /silent)
            sxdelpar, hdr, 'XTENSION'
            sxaddpar, hdr, 'XDANGL', sxpar(hdr0, 'XDANGL')
            sxaddpar, hdr, 'ECHANGL', sxpar(hdr0, 'ECHANGL')
            mwrfits, a, 'hires'+string(imlist[i], format='(I04)')+'_'+string(j, format='(I1)')+'.fits', hdr, /create
        endfor

    endfor
end
