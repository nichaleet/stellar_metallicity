pro out2bin, f, newmoog=newmoog
    moogspec = read_moog_spec(f, newmoog=newmoog)
    lambdafile = file_dirname(f, /mark_directory)+'lambda.bin'
    if ~file_test(lambdafile) then begin
        openw, 2, lambdafile
        writeu, 2, moogspec.lambda
        close, 2
    endif
    binfname = file_dirname(f, /mark_directory)+file_basename(f, '.out2')+'.bin'
    file_delete, binfname, /allow_nonexistent
    openw, 1, binfname
    writeu, 1, moogspec.spec
    close, 1
    file_delete, f
    file_delete, binfname+'.gz', /allow_nonexistent
    spawn, 'gzip --rsyncable '+binfname
end


pro convert_grid1
    f = file_search(getenv('ahome')+'grid/synths/*.out2', count=c)
    for i=0L,c-1 do begin
        out2bin, f[i]
    endfor
end
