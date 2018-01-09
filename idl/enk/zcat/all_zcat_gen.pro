pro all_zcat_gen
    fzcat = full_zcat()
    sczcat = ucsc_zcat()
    match, fzcat.objno, sczcat.objno, w1, w2
    fzcat = fzcat[complement(w1, n_elements(fzcat))]
    all_zcat = [sczcat, fzcat]
    all_zcat = all_zcat[uniq(all_zcat.objno, sort(all_zcat.objno))]
    mwrfits, all_zcat, 'all_zcat.fits'
end
