;1 - ordinary least squares, y vs. x
;2 - ordinary least squares, x vs. y
;3 - bisector
;4 - orthogonal regression

function ab96, x, y, dx, dy, type, sigma=sigma, forceint=int, forceslope=forceslope
    if type ne 1 and type ne 2 and type ne 3 and type ne 4 then message, 'Incorrect fit type.'

    if 0 then begin
        mx = weightedmean(x, dx, error=dmx)
        my = weightedmean(y, dy, error=dmy)
        covxy = total((dx*dy)^(-1.))*total((dx*dy)^(-1.)*(x-mx)*(y-my)) / (total((dx*dy)^(-1.))^2. - total((dx*dy)^(-2.)))
        ddx = weightedmean(dx^2., dx)
        ddy = weightedmean(dy^2., dy)
    endif else begin
        mx = mean(x)
        my = mean(y)
        covxy = total((x-mx)*(y-my))/double(n_elements(x)-1)
        dmx = stddev(x)
        dmy = stddev(y)
        ddx = mean(dx^2.)
        ddy = mean(dy^2.)
    endelse

    ;beta1 = total((x-mx)*(y-my))/(total((x-mx)^2.) - total(dx^2.))
    ;beta2 = (total((y-my)^2.) - total(dy^2.))/total((x-mx)*(y-my))

    beta1 = covxy / (dmx^2. - ddx)
    beta2 = (dmy^2. - ddy) / covxy
    case type of
        1: beta = keyword_set(forceslope) ? forceslope : beta1
        2: beta = keyword_set(forceslope) ? forceslope : beta2
        3: beta = keyword_set(forceslope) ? forceslope : (beta1*beta2 - 1.0 + sqrt((1.0+beta1^2.)*(1.0+beta2^2.))) / (beta1 + beta2)
        4: beta = keyword_set(forceslope) ? forceslope : 0.5*((beta2-beta1^(-1.)) + sign(covxy)*sqrt(4.0+(beta2-beta1^(-1.))^2.))
    endcase
    alpha = keyword_set(forceint) ? int : my - beta*mx
    db1 = ((x-mx)*(y-beta1*x-alpha) + beta1*dx^2.) / (dmx^2. - ddx)
    db2 = ((y-my)*(y-beta2*x-alpha) - dy^2.) / covxy
    case type of
        1: db = db1
        2: db = db2
        3: db = ((1.0+beta2^2.)*beta*db1 + (1.0+beta1^2.)*beta*db2)/((beta1+beta2)*sqrt((1.0+beta1^2.)*(1.0+beta2^2.)))
        4: db = ((beta/beta1^2.)*db1+beta*db2)/sqrt(4.0+(beta2-beta1^(-1.))^2.)
    endcase
    da = y - beta*x - mx*db
    dbeta = keyword_set(forceslope) ? 0.0 : (mean((db-mean(db))^2.)/n_elements(db))^0.5
    dalpha = keyword_set(forceint) ? 0.0 : (mean((da-mean(da))^2.)/n_elements(da))^0.5
    sigma = [dalpha, dbeta]
    return, [alpha, beta]
end
