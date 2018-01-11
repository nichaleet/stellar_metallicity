function contains, array, x
    w = where(array eq x, c)
    if c gt 0 then return, 1 else return, 0
end
