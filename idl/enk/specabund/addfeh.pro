pro addfeh, moogify
    ;constant = 0.15855343
    constant = 0.0
    tags = tag_names(moogify)
    if contains(tags, 'FEH') then moogify.feh += constant
    ;if contains(tags, 'MGFE') then moogify.mgfe -= 0.15
    if contains(tags, 'SIFE') then moogify.sife -= constant
    if contains(tags, 'SIFE') then moogify.sife += 0.15
    if contains(tags, 'CAFE') then moogify.cafe -= constant
    if contains(tags, 'TIFE') then moogify.tife -= constant
    if contains(tags, 'ALPHAFE') then moogify.alphafe -= constant
    if contains(tags, 'ALPHAFE2') then moogify.alphafe2 -= constant
    if contains(tags, 'ALPHAFE3') then moogify.alphafe3 -= constant
    if contains(tags, 'ALPHAFE4') then moogify.alphafe4 -= constant
end
