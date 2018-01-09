pro plotlines,z,color=color
  readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt',linemin, linemax, type,name,format ='D,D,A,A' 
  linewave = 0.5*(linemin+linemax)
  nlines = n_elements(linewave)
  linewave = linewave*(z+1.)
  linename = strsplit(name,'#',/extract)
  goodlines = where(type eq 'a',cgood)
  linewave=linewave(goodlines)
  linename=linename(goodlines)
  for j=0,cgood-1 do if linewave(j) ge !x.crange[0] and linewave(j) lt !x.crange[1] then begin
     ymin = (!y.crange)[0]
     ymax = (!y.crange)[1]*0.5
     if ~keyword_set(color) then cnow = 'red' else cnow = color
     oplot,[linewave[j],linewave[j]],[ymin,ymax],color=160,linestyle=2 
     xyouts,linewave[j],0.5,linename[j],color=160,font=0,charsize=0
  endif
end
