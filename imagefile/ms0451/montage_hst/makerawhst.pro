pro makerawhst
files=file_search('hst/*.fits',count=cfile)
for i=0,cfile-1 do begin
   im = mrdfits(files[i],1,hdr)
   name = strmid(files[i],4)
   writefits,'rawdir/nl_'+name,im,hdr
endfor

end
