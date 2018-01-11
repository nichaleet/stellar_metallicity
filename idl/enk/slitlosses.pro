function gaussian2d, x, y
    common com, seeing, slitlength
    return, exp((-x*x - y*y) / (2d * seeing*seeing))
end

function pq_limits, x
    common com, seeing, slitlength
    return, [0, slitlength/2]
end

function slitlosses, slitwidth, seeingin
    common com, seeing, slitlength
    slitlength = 120.0
    a = 2d * sqrt(2d * alog(2d))
    seeing = seeingin / a
    total = int_2d('gaussian2d', [0, 100.], 'pq_limits', 96, /double)
    truncated = int_2d('gaussian2d', [0, slitwidth], 'pq_limits', 96, /double)
    loss = truncated / total
    return, loss
end
