function append_struct, tx, ty
    tags = tag_names(tx)
    s = 'tmp = create_struct("' + tags[0] + '", [tx.' + tags[0] + $
        ', ty.' + tags[0] + '])'
    void = execute(s)
    for i = 1, n_elements(tags) - 1 do begin
        s = 'tmp = create_struct(tmp, "' + tags[i] + '", [tx.' + tags[i] + $
            ', ty.' + tags[i] + '])'
        void = execute(s)
    endfor
    return, tmp
end
