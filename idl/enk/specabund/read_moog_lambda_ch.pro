function read_moog_lambda_ch, n=n, fullres=fullres, mds=mds
    lambdafile = getenv('chome')+'gridch'+(keyword_set(mds) ? '_mds/' : '/')+(keyword_set(fullres) ? 'synths' : 'bin')+'/lambda.bin'
    lambda = dblarr(20001)
    lambda_in = 0d
    n = 0L
    openr, 1, lambdafile
    while ~eof(1) do begin
        readu, 1, lambda_in
        lambda[n] = lambda_in
        n += 1
    endwhile
    close, 1
    lambda = lambda[0:n-1]
    return, lambda
end
