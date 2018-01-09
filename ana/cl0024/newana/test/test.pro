pro test
readcol,'POTTHOFF.dat',y,x,id
idea = where(id eq 1)
nidea = where(id eq 0)
idea_param =linfit(x(idea),y(idea),sigma=sigmai,yfit=yfiti)
nidea_param = linfit(x(nidea),y(nidea),sigma=sigman,yfit=yfitn)

print, idea_param
print, nidea_param
stop
end
