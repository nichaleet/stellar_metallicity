;+
; NAME:
;       weightedmean
;
; PURPOSE:
;       Calculate the mean of a set of measurements downweighted by
;       the standard error in each measurement.
;
; CALLING SEQUENCE:
;       result = weightedmean(x, sd)
;
; INPUT:
;       x       an array of measurement values
;       sd      an array of standard error which must be the same
;               length as x
;
; OPTIONAL OUTPUT:
;       error   standard error of the weighted mean
;-

function weightedmean, x, sd, error=error, meanerror=meanerror
     N = n_elements(x)
     if n_elements(sd) ne N then message, 'SD must be same length as X'
     xmean = total(x/sd^2.) / total(1./sd^2.)
     meanerror = 1. / sqrt(total(1./sd^2.))
     error = sqrt((total(x^2./sd^2.)*total(1./sd^2.) - (total(x/sd^2.))^2.) / ((total(1./sd^2.))^2. - total(1./sd^4.)))
     return, xmean
end
