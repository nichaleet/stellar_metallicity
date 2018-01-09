pro scr
;random
;see the original models

  common sps_spec, sps, spsz, spsage
  spsspec = sps_interp(0.0, 5.0)

  for i=10,21 do begin
     for j=8,18 do begin
        z = spsz(i)
        age = spsage(j*10)
        spsstruct = sps_interp(z, age)
        lambda = spsstruct.lambda
        spsspec = spsstruct.spec
        spsspec = spsspec/lambda^2
        good = where(lambda gt 3000 and lambda lt 8000)
        plot,lambda(good),spsspec(good)
        xyouts,0.1,0.1,'z='+strtrim(string(z),2)+ 'age ='+strtrim(string(age),2),/normal
        wait,1
     endfor
  endfor
stop
end
