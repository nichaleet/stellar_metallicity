pro find_candidates_clash,clustername
filein = 'clash_'+clustername+'_cat.sav'
dir    ='/scr2/nichal/workspace2/red_arcs/catalogs/CLASH/'
filein = dir+filein
restore,filein
dim    = size(Data,/dimensions)
if dim[0] eq 125 then begin
   f475 = {wl:475,mag:data[43,*],magerr:data[44,*]}
   f555 = {wl:555,mag:data[49,*],magerr:data[50,*]}
   f606 = {wl:606,mag:data[55,*],magerr:data[56,*]}
   f625 = {wl:625,mag:data[61,*],magerr:data[62,*]}
   f775 = {wl:775,mag:data[67,*],magerr:data[68,*]}
   f814 = {wl:814,mag:data[73,*],magerr:data[74,*]}
   f105 = {wl:1050,mag:data[85,*],magerr:data[86,*]}
   f110 = {wl:1100,mag:data[91,*],magerr:data[92,*]}
   f125 = {wl:1250,mag:data[97,*],magerr:data[98,*]}
   f140 = {wl:1400,mag:data[103,*],magerr:data[104,*]}
   f160 = {wl:1600,mag:data[109,*],magerr:data[110,*]}
   mag = [f475,f555,f606,f625,f775,f814,f105,f110,f125,f140,f160]
   i_flux = reform(data[70,*])
   i_mag  = reform(data[67,*])
   g_mag  = reform(data[43,*])
   z = reform(data[123,*])
   bpz = reform(data[115,*])
endif else if dim[0] eq 119 then begin
   f475 = {wl:475,mag:data[43,*],magerr:data[44,*]}
   f606 = {wl:606,mag:data[49,*],magerr:data[50,*]}
   f625 = {wl:625,mag:data[55,*],magerr:data[56,*]}
   f775 = {wl:775,mag:data[61,*],magerr:data[62,*]}
   f814 = {wl:814,mag:data[67,*],magerr:data[68,*]}
   f105 = {wl:1050,mag:data[79,*],magerr:data[80,*]}
   f110 = {wl:1100,mag:data[85,*],magerr:data[86,*]}
   f125 = {wl:1250,mag:data[91,*],magerr:data[92,*]}
   f140 = {wl:1400,mag:data[97,*],magerr:data[98,*]}
   f160 = {wl:1600,mag:data[103,*],magerr:data[104,*]}
   mag = [f475,f606,f625,f775,f814,f105,f110,f125,f140,f160]
   i_Flux = reform(data[64,*])
   i_mag  = reform(data[61,*])
   g_mag  = reform(data[43,*])
   z = reform(data[117,*])
   bpz = reform(data[109,*])
endif else stop,'Stop(41): Number of Columns is unknown.'

;g-i selection
lessgood = where(z gt 0.6 and g_mag-i_mag gt 1. and i_mag lt 21.5 and i_mag ne -99 and z lt 1.2)

;UVJ selection
;loop by redshift interval
;z_list = [0.4,0.5,0.6,0.7,0.8,0.9,1.0]
z_list = [0.6,0.7,0.8,0.9,1.0]
goodarr=[]
for i=0,n_elements(z_list)-1 do begin
   znow = z_list[i]
   findu = min(abs((1.+znow)*365-mag.wl),upos)
   findv = min(abs((1.+znow)*550-mag.wl),vpos)
   findj = min(abs((1.+znow)*1200-mag.wl),jpos)
   umag = mag(upos).mag
   vmag = mag(vpos).mag
   check = where(upos-vpos eq 0,c)
   if c ne 0 then begin
      print, 'check if u and v fall in the same band'
      print, c
   endif
   if findj lt 150. then begin
      jmag = mag(jpos).mag
      good = where(i_mag lt 21.5 and i_mag gt 18.0 and z gt znow and z lt znow+0.1 and umag-vmag gt 1.3 and vmag-jmag lt 1.5 and umag-vmag gt (vmag-jmag)*0.88+0.69 and abs(bpz-z) lt 0.1,c)
   endif else good =  where(i_mag lt 21.5 and i_mag gt 18.0 and z gt znow and z lt znow+0.1 and umag-vmag gt 1.3 and abs(bpz-z) lt 0.1,c)
   if c ne 0 then goodarr=[goodarr,good]
endfor
print, clustername
print, 'Found total of', n_elements(goodarr),' candidates'
goodarr_ori = goodarr
goodarr     = setunion(goodarr_ori,lessgood)
if clustername eq 'macs0429' then goodarr = [goodarr,1107,1311]
findgood    = setintersection(goodarr,goodarr_ori,positions=goodpos)
goodtag     = bytarr(n_elements(goodarr))
goodtag(goodpos) = 1
;save structure
if n_elements(goodarr) gt 0 then begin
   ids = long(data[0,goodarr])
   ra = data[1,goodarr]
   dec = data[2,goodarr]
   if dim[0] eq 125 then begin
      f475 = {wl:475,mag:data[43,goodarr],magerr:data[44,goodarr]}
      f555 = {wl:555,mag:data[49,goodarr],magerr:data[50,goodarr]}
      f606 = {wl:606,mag:data[55,goodarr],magerr:data[56,goodarr]}
      f625 = {wl:625,mag:data[61,goodarr],magerr:data[62,goodarr]}
      f775 = {wl:775,mag:data[67,goodarr],magerr:data[68,goodarr]}
      f814 = {wl:814,mag:data[73,goodarr],magerr:data[74,goodarr]}
      f105 = {wl:1050,mag:data[85,goodarr],magerr:data[86,goodarr]}
      f110 = {wl:1100,mag:data[91,goodarr],magerr:data[92,goodarr]}
      f125 = {wl:1250,mag:data[97,goodarr],magerr:data[98,goodarr]}
      f140 = {wl:1400,mag:data[103,goodarr],magerr:data[104,goodarr]}
      f160 = {wl:1600,mag:data[109,goodarr],magerr:data[110,goodarr]}
      mag = [f475,f555,f606,f625,f775,f814,f105,f110,f125,f140,f160]
      i_flux = reform(data[70,goodarr])
      g_mag  = reform(data[43,goodarr])
      r_mag  = reform(data[61,goodarr])
      i_mag  = reform(data[67,goodarr])
      z = reform(data[123,goodarr])
      bpz = reform(data[115,goodarr])
   endif else if dim[0] eq 119 then begin
      f475 = {wl:475,mag:data[43,goodarr],magerr:data[44,goodarr]}
      f606 = {wl:606,mag:data[49,goodarr],magerr:data[50,goodarr]}
      f625 = {wl:625,mag:data[55,goodarr],magerr:data[56,goodarr]}
      f775 = {wl:775,mag:data[61,goodarr],magerr:data[62,goodarr]}
      f814 = {wl:814,mag:data[67,goodarr],magerr:data[68,goodarr]}
      f105 = {wl:1050,mag:data[79,goodarr],magerr:data[80,goodarr]}
      f110 = {wl:1100,mag:data[85,goodarr],magerr:data[86,goodarr]}
      f125 = {wl:1250,mag:data[91,goodarr],magerr:data[92,goodarr]}
      f140 = {wl:1400,mag:data[97,goodarr],magerr:data[98,goodarr]}
      f160 = {wl:1600,mag:data[103,goodarr],magerr:data[104,goodarr]}
      mag = [f475,f606,f625,f775,f814,f105,f110,f125,f140,f160]
      i_Flux = reform(data[64,goodarr])
      g_mag  = reform(data[43,goodarr])
      r_mag  = reform(data[55,goodarr])
      i_mag  = reform(data[61,goodarr])
      z = reform(data[117,goodarr])
      bpz = reform(data[109,goodarr])
   endif
   struct = {ids:ids,g_mag:g_mag,r_mag:r_mag,i_mag:i_mag,z:z,bpz:bpz,ra:ra,dec:dec,goodtag:goodtag}
      
      print, '     id     i_mag    g-i      z'   
      for i=0,n_elements(ids)-1 do begin
      print, strtrim(string(ids(i)),2),i_mag(i),g_mag(i)-i_mag(i),z(i)
      endfor
      comment_Short = strtrim(string(reform(ids)),2)+','+strtrim(string(long(goodtag)),2)
      comments = 'i='+strtrim(string(i_mag),2)+', z='+strtrim(string(z),2)+',g-i='+strtrim(string(g_mag-i_mag),2) 
     write_ds9_regionfile,ra,dec,comment=comments,filename='/scr2/nichal/workspace2/red_arcs/imaging_data/'+clustername+'/'+clustername+'_candidate.reg',color='blue'
      write_ds9_regionfile,ra,dec,comment=comment_short,filename='/scr2/nichal/workspace2/red_arcs/imaging_data/'+clustername+'/'+clustername+'_candidate_id.reg',color='blue'
      save,struct,mag,filename=dir+clustername+'_candidate.sav'
   endif else print, 'no good candidates in this clusters'
stop
end


