pro makerawhst_lowres
files=file_search('/scr2/nichal/workspace2/imagefile/ms0451/hst/*.fits',count=cfile)
for i=0,cfile-1 do begin
   im = mrdfits(files[i],1,hdr)
   name = file_basename(files[i])
   oldsize = size(im,/dimension)
   newx = round(oldsize[0]/4.)
   newy = round(oldsize[1]/4.)
   hrebin,im,hdr,newim,newhdr,newx,newy
   writefits,'rawdir/nl_lowres_'+name,newim,newhdr
endfor

end
