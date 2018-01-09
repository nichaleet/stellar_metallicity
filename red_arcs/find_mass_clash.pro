pro find_mass_clash,cluster,id,mu
;Find mass using Bell et al 2013
mu=float(mu)
file=file_search('/scr2/nichal/workspace2/red_arcs/catalogs/CLASH/'+cluster+'_candidate.sav',count=c)
if c eq 0 then stop,'Stop: no file found'
restore,file
s=struct
m=mag

pos = where(s.ids eq id)
if pos eq -1 then stop
z   = s.z(pos)
bpz = s.bpz(pos)
ra  = s.ra(pos)
dec = s.dec(pos)
mag = reform(m.mag(pos))+2.5*alog10(mu)
magerr = reform(m.magerr(pos))
mag_ivar = 1./magerr^2
wl  = m.wl
badmag = where(mag lt 0 or mag gt 50)
remove,badmag,wl,mag,magerr,mag_ivar

file_kcorrect='/scr2/nichal/workspace2/red_arcs/k_correct/'+cluster+'_'+strtrim(string(id),2)+'.sav'
if ~file_test(file_kcorrect) then begin
;kcorrection
   for i=0,n_elements(wl)-1 do begin
      wlnow = wl(i)
      if wlnow lt 1000 then filter = file_search('/scr2/nichal/kcorrect/data/filters/clash_acs_f'+strtrim(string(wlnow),2)+'w.par',count=c)
      if wlnow gt 1000 then begin
         wlnow=wlnow/10
         filter = file_search('/scr2/nichal/kcorrect/data/filters/clash_wfc3_f'+strtrim(string(wlnow),2)+'w.par',count=c)
      endif
      if c eq 0 then stop
      filter = strmid(filter,1+strpos(filter,'/',/reverse_search))
      if i eq 0 then filterlist = filter else filterlist=[filterlist,filter]
   endfor
   
   kcorrect,mag,mag_ivar,z(0),kcorrect1,chi2=chi2,filterlist=filterlist,/magnitude
   
;griz
   sdss_wl = [477,623,762,913]
   good = intarr(4)
   for fil=0,3 do  good(fil) = where(abs(sdss_wl(fil)-wl) eq min(abs(sdss_wl(fil)-wl)))
   kcorrect,mag(good),mag_ivar(good),z(0),kcorrect2,chi2=chi2,filterlist=['sdss_g0.par', 'sdss_r0.par', 'sdss_i0.par', 'sdss_z0.par'],/magnitude

   print, kcorrect1
   print, kcorrect2
   save,kcorrect1,kcorrect2,good,filename=file_kcorrect
endif else restore,file_kcorrect

;now calculate absolute griz
z=z(0)
dist = lumdist(z)*1.e6        ;parsecs
abs_griz = mag(good)-5.*(alog10(dist)-1.)+kcorrect2
g = abs_griz(0)
r = abs_griz(1)
i = abs_griz(2)
z = abs_griz(3)
bandtags=['g','r','i','z']
color=[g-r,g-i,g-z,r-i,r-z]
colortags=['g-r','g-i','g-z','r-i','r-z']
solar_abs=[5.15,4.67,4.56,4.53] ;griz from Bell2003
LogL  = (solar_abs-abs_griz)/2.5 ;10^logL

a_arr = -0.001*[[499,379,367,106,124],[306,220,215,22,41],[222,152,153,-6,18],[223,175,171,52,41]]
b_arr = 0.001*[[1519,914,698,1982,1067],[1097,661,508,1431,780],[864,518,402,1114,623],[689,421,322,923,463]]

logM_arr = fltarr(5,4)
;loop over bands and colors to find Mass
for cc=0,4 do begin 
   for bb=0,3 do begin
      a=a_arr[cc,bb]
      b=b_arr[cc,bb]
      logM_arr[cc,bb] = logL[bb]+a+b*color(cc) 
   endfor ;loop over band
endfor ;loop over color
print, logm_arr
print, 'median logM is',median(logm_arr)
print, 'median logL is',median(logL)
stop
end
