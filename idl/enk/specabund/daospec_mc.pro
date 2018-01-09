pro daospec_mc, name, processors=processors
    if ~keyword_set(processors) then processors = 8
    
    ;case name of
    ;    'hd115444': begin
    ;        minvr = -28
    ;        maxvr = -24
    ;    end
    ;    'hd122563': begin
    ;        minvr = -25
    ;        maxvr = -23
    ;    end
    ;    'UMi20103': begin
    ;        minvr = -245
    ;        maxvr = -242
    ;    end
    ;    'Scl1019417': begin
    ;        minvr = 100
    ;        maxvr = 103
    ;    end
    ;    else: message, "I don't know the RV limits for this star."
    ;endcase
    
    nmc = 100
    seed = 726405L
    for i=0L,nmc-1 do begin
        print, 'Trial '+string(i+1, format='(I5)')+' / '+string(nmc, format='(I5)')
        spawn, 'ps -e | grep daospec_hires | wc', numproc
        shouldiwait = fix(strmid(numproc, 0, 7)) ge processors
        while shouldiwait do begin
            spawn, 'ps -e | grep daospec_hires | wc', numproc
            shouldiwait = fix(strmid(numproc, 0, 7)) ge processors
            wait, 1
        endwhile
        split_hires, name, /daospec, mcnum=i+1, seed=seed;, minvr=minvr, maxvr=maxvr
        spawn, 'cp *_0??e.fits '+string(i+1, format='(I05)')+'/.'
        cd, string(i+1, format='(I05)')
        spawn, 'daospec_hires '+name+' &'
        cd, '..'
    endfor
    
    print, 'Waiting for daospec to finish.'
    shouldiwait = 1
    while shouldiwait do begin
        spawn, 'ps -e | grep daospec_hires | wc', numproc
        shouldiwait = fix(strmid(numproc, 0, 7)) ge 1
        wait, 1
    endwhile
    ;for i=0,nmc-1 do begin
    ;    daospec_to_moog, name, /jlc, mctrial=i+1
    ;endfor
    for i=0,nmc-1 do begin
        print, 'Trial '+string(i+1, format='(I5)')+' / '+string(nmc, format='(I5)')
        merge_hires, name, mc=i+1
        cd, string(i+1, format='(I05)')
        spawn, 'rm *_0??e.fits *_0??.fits *_0??C.fits *_0??R.fits *.opt'
        cd, '..'
    endfor

    exit
end
