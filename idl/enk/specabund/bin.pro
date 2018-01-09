function bin, array1, binsize
    binsize = long(binsize)
    n = n_elements(array1)/binsize + 1
    w = (lindgen(binsize, n) + 1) < (n_elements(array1) - 1)
    return, total(array1[w], 1) / binsize
end
