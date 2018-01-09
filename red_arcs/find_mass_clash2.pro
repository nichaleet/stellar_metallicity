pro find_mass_clash2,cluster,id,mu,redo=redo
;Find mass using Kcorrect program with full HST band
mu=float(mu)
file=file_search('/scr2/nichal/workspace2/red_arcs/catalogs/CLASH/'+cluster+'_candidate.sav',count=c)
if c eq 0 then stop,'Stop: no file found'
restore,file
s=struct
m=mag

pos = where(s.ids eq id)
if pos eq -1 then stop
z   = (s.z(pos))[0]
bpz = s.bpz(pos)
ra  = s.ra(pos)
dec = s.dec(pos)
mag = reform(m.mag(pos))+2.5*alog10(mu)
magerr = reform(m.magerr(pos))
mag_ivar = 1./magerr^2
wl  = m.wl
badmag=where(mag gt 50 or mag lt 0,cbad)
if cbad ne 0 then remove,badmag,mag,magerr,mag_ivar,wl

file_kcorrect='/scr2/nichal/workspace2/red_arcs/k_correct/'+cluster+'_'+strtrim(string(id),2)+'.sav'
if ~file_test(file_kcorrect) or keyword_set(redo) then begin
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
   badmag = where(mag lt 0,cbad);case of -99 and stuff
   if cbad ne 0 then remove,cbad,mag,mag_ivar,filterlist
   kcorrect,mag,mag_ivar,z,kcorrect1,chi2=chi2,filterlist=filterlist,/magnitude,mass=mass,b300=b300

   print, kcorrect1
   mass=0.52*mass ;fix for the hubble constant h=0.72
   save,kcorrect1,mass,b300,filterlist,filename=file_kcorrect
endif else restore,file_kcorrect

;find rest frame g-i
wl = long(wl)
findg = min(abs(wl-475.),posg)
findi = min(abs(wl-762.),posi)
if findg lt 150. and findi lt 150. then g_i=(mag(posg)-kcorrect1(posg))-(mag(posi)-kcorrect1(posi))
;m = M+DM(z)+K(z)
;note: the magnitude is corrected for mu but it's fine since it's a constant.
print, 'M is',mass
print, 'LogM is',alog10(mass)
print, 'rest frame g-i:',g_i,wl([posg,posi])
print, 'b300:',b300
print, 'imag' ,mag(posi)-2.5*alog10(mu)
stop
end
