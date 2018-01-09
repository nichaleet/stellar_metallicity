;--------------------------------------------------------------------
;+
; NAME:
;
; 
; PURPOSE:
;
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
;
; MODIFICATION HISTORY:
;~~~mg
;- 
;--------------------------------------------------------------------

pro ks2d2s, x1,y1,x2,y2,d,probks

    n1 = n_elements(x1)
    n2 = n_elements(x2)

    d1=0
    for j=0L,n1-1 do begin
       quadct,x1[j],y1[j],x1,y1,n1,fa,fb,fc,fd
       quadct,x1[j],y1[j],x2,y2,n2,ga,gb,gc,gd
       d1 = max([d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd)])
    endfor

    d2=0
    for j=0L,n2-1 do begin
       quadct,x2[j],y2[j],x1,y1,n1,fa,fb,fc,fd
       quadct,x2[j],y2[j],x2,y2,n2,ga,gb,gc,gd
       d2 = max([d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd)])
    endfor

    d = 0.5*float(d1+d2)

    sqen = sqrt(float(n1*n2)/float(n1+n2))
    r1 = correlate(x1,y1)
    r2 = correlate(x2,y2)

    rr=sqrt(1.0 - 0.5*(r1^2. + r2^2.))

    d_in  = d*sqen/(1.0 + rr*(0.25-0.75/sqen))

    ;neff=float(n1*n2)/float(n1+n2)

    ;prob_ks, d_in, neff, probks
    ;prob_ks, rr, sqen, probks

    a2 = -2.*d_in^2.
    fac = 2
    probks = 0
    termbf = 0
    eps1 = 1.e-8
    eps2 = 1.e-10
    for j=1,100 do begin
       term = fac*exp(a2*j^2.)
       probks = probks + term
       if abs(term) le eps1*termbf or abs(term) le eps2*probks then return
       fac = -fac
       termbf = abs(term)
    endfor
    probks = 1
    return
 end


