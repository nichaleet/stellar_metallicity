pro lick_measurement,outname,ewfile,plot=plot
;lick_measurement,'all_cl0024_gauss_nomaxflux_emisssub','all_cl0024_lickindices_gauss_nomaxflux_emisssub.fits',/plot

if ~keyword_set(plot) then plot=0 else plot=1
;output directory
outdir = '/scr2/nichal/workspace2/lickindices/output/results/'
!p.multi= [0,1,3]
!p.font = 0
;reading EW files
ewdir = '/scr2/nichal/workspace2/lickindices/output/ew/'
ewfile = ewdir+ewfile
ew = mrdfits(ewfile,1)
ngals = n_elements(ew)
goodobj = where(ew.good eq 1,cgoodobj)

;get the thomas11b prediction grid at specified resolution
res=0.025
;check the file
filename = '/scr2/nichal/workspace2/lickindices/thomas_fine'+strtrim(string(res,format='(F5.3)'),2)+'.sav'
filethomas = file_search(filename,count=fcount)
if fcount eq 0 then begin
   thomasgrid = thomas_interp(res) 
   save,thomasgrid,filename='/scr2/nichal/workspace2/lickindices/thomas_fine'+strtrim(string(res,format='(F5.3)'),2)+'.sav'
endif else restore,filethomas[0]

ngrids = size(thomasgrid,/dimension)
nage = ngrids(0)
nz = ngrids(1)
nalpha = ngrids(2)
agethomas = thomasgrid.age ; log(gyr)
zthomas = thomasgrid.zh    ;log
alphathomas = thomasgrid.alphafe ;log

tagthomas = tag_names(thomasgrid)
tagmatch  = [where(strmatch(tagthomas,'cn1',/fold_case) eq 1),where(strmatch(tagthomas,'cn2',/fold_case) eq 1),where(strmatch(tagthomas,'ca4227',/fold_case) eq 1),where(strmatch(tagthomas,'g4300',/fold_case) eq 1),where(strmatch(tagthomas,'c24668',/fold_case) eq 1),where(strmatch(tagthomas,'mg1',/fold_case) eq 1),where(strmatch(tagthomas,'mg2',/fold_case) eq 1),where(strmatch(tagthomas,'mgb',/fold_case) eq 1),where(strmatch(tagthomas,'fe4383',/fold_case) eq 1),where(strmatch(tagthomas,'fe4531',/fold_case) eq 1),where(strmatch(tagthomas,'fe5270',/fold_case) eq 1),where(strmatch(tagthomas,'fe5335',/fold_case) eq 1),where(strmatch(tagthomas,'fe5406',/fold_case) eq 1),where(strmatch(tagthomas,'hda',/fold_case) eq 1),where(strmatch(tagthomas,'hga',/fold_case) eq 1),where(strmatch(tagthomas,'hgf',/fold_case) eq 1)]
if n_elements(tagmatch) ne 16 then stop, 'Oh no'

;OUTPUT STRUCTURE
;make output structure which is a copy of EW but add more fields
outarr=[]
for i=0,ngals-1 do begin
   ew_new=ew[i]
   ew_new=create_struct(ew_new,'chisq',fltarr(nage,nz,nalpha)-99.)
   ew_new=create_struct(ew_new,'prob_age',fltarr(nage,2)-99.)
   ew_new=create_struct(ew_new,'prob_ZH',fltarr(nz,2)-99.)
   ew_new=create_struct(ew_new,'prob_alpha',fltarr(nalpha,2)-99.)
   ew_new=create_struct(ew_new,'best_age',fltarr(4)-99.)
   ew_new=create_struct(ew_new,'best_ZH',fltarr(4)-99.)
   ew_new=create_struct(ew_new,'best_alpha',fltarr(4)-99.)
   ew_new=create_struct(ew_new,'keepindex',bytarr(16))
   outarr = [outarr,ew_new]
endfor

;Run loop over each object
for io=0,cgoodobj-1 do begin
   obj = ew(goodobj(io))
   print, 'Now doing object: ',obj.objname
   print, 'Mask: ',obj.mask
   print, 'object number: ',io+1,'/',cgoodobj
   index = [obj.cn1[0],obj.cn2[0],obj.ca4227[0],obj.g4300[0],obj.c24668[0],obj.mg1[0],obj.mg2[0],obj.mgb[0],obj.fe4383[0],obj.fe4531[0],obj.fe5270[0],obj.fe5335[0],obj.fe5406[0],obj.hda[0],obj.hga[0],obj.hgf[0]]
   index_err = [obj.cn1[1],obj.cn2[1],obj.ca4227[1],obj.g4300[1],obj.c24668[1],obj.mg1[1],obj.mg2[1],obj.mgb[1],obj.fe4383[1],obj.fe4531[1],obj.fe5270[1],obj.fe5335[1],obj.fe5406[1],obj.hda[1],obj.hga[1],obj.hgf[1]]
   keepindex= bytarr(16)+1
   chisqcheck=0
   cutcount=0
   cutcountmax = 2
   while chisqcheck eq 0 and cutcount le cutcountmax do begin
      goodind = where(index ne -99 and index_err ne -99 and keepindex eq 1,cgoodind) ;goodind of observed spec
      if cgoodind gt 0 then begin ;only do if there are good index measurements
    
         index = index(goodind)
         index_err = index_err(goodind)
         tagmatchnow = tagmatch(goodind)
         ;;calculate chisq array
         chisqarr = fltarr(nage,nz,nalpha)
         for ii=0,cgoodind-1 do begin ;loop over each index
            residual = (index(ii)-thomasgrid.(tagmatchnow(ii)))^2/(index_err(ii))^2
            chisqarr += residual ;chisquare array of different models
         endfor      
         ;;normalize the probability
         minchisq = min(chisqarr)
         prob = exp(-0.5*(chisqarr-minchisq))      
         ;;integration
         volume = 0d
         for iage=0L,nage-2 do begin
            for iz=0L,nz-2 do begin
               for ialpha=0L,nalpha-2 do begin
                  dage = agethomas[iage+1,iz,ialpha]-agethomas[iage,iz,ialpha]
                  dz = zthomas[iage,iz+1,ialpha]-zthomas[iage,iz,ialpha]
                  dalpha = alphathomas[iage,iz,ialpha+1]-alphathomas[iage,iz,ialpha]
                  volume+= dage*dz*dalpha*prob(iage,iz,ialpha)
               endfor
            endfor
         endfor
         print, 'volume under surface is ',volume
         prob = prob/volume
         ;;finish normalizing the probability
         ;;Get 1D prob for each parameter
         probage = fltarr(nage)
         probz   = fltarr(nz)
         probalpha = fltarr(nalpha)
         for i=0,nage-1 do probage(i) = int_tabulated_2d(zthomas[i,*,*],alphathomas[i,*,*],prob[i,*,*])
         for i=0,nz-1 do probz(i) = int_tabulated_2d(agethomas[*,i,*],alphathomas[*,i,*],prob[*,i,*])
         for i=0,nalpha-1 do probalpha(i) = int_tabulated_2d(agethomas[*,*,i],zthomas[*,*,i],prob[*,*,i])
         ;;Get confidence interval for each parameter
      
         bestage = confidence_interval(agethomas[*,0,0],probage,xtitle='log(age)')
         bestz = confidence_interval(zthomas[0,*,0],probz,xtitle='[Z/H]')
         bestalpha = confidence_interval(alphathomas[0,0,*],probalpha,xtitle='[alpha/Fe]')
         
         ;;check chisquare contribution of each index
         maxprob = max(probage,ageloc)
         maxprob = max(probz,zloc)
         maxprob = max(probalpha,alphaloc)
         bestchiarr = fltarr(16)-99.
         for ii=0,cgoodind-1 do begin ;loop over each index
            bestchiarr(goodind(ii)) = (index(ii)-thomasgrid[ageloc,zloc,alphaloc].(tagmatchnow(ii)))^2/(index_err(ii))^2
         endfor
         ;;if the contribution from the index with the most contribution is greater than 30% then take it out.
         if max(bestchiarr,locmax)/total(bestchiarr) gt 0.3 then begin
            keepindex(locmax)= 0
            cutcount+=1
         endif else chisqcheck=1
         if cutcount eq cutcountmax then chisqcheck=1
         
      endif
   endwhile
   if chisqcheck eq 1 then begin
      print,'Cut ', cutcount, ' indices.'
      if cutcount gt 0 then print,tagthomas(tagmatch(where(keepindex eq 0)))
      ;;plot
      if plot eq 1 then begin
         plot,10^(agethomas[*,0,0]-9.),probage,xtitle='Age(Gyr)',/xlog
         for ip=1,3 do oplot,10^([bestage[ip],bestage[ip]]-9.),!y.crange,color=fsc_color('red'),linestyle=2
         oplot,10^([bestage[0],bestage[0]]-9.),!y.crange,color=fsc_color('green'),linestyle=2
         plot,zthomas[0,*,0],probz,xtitle='[z/H]'
         for ip=1,3 do oplot,[bestz[ip],bestz[ip]],!y.crange,color=fsc_color('red'),linestyle=2
         oplot,[bestz[0],bestz[0]],!y.crange,color=fsc_color('green'),linestyle=2	
         plot,alphathomas[0,0,*],probalpha,xtitle='[alpha/Fe]'
         for ip=1,3 do oplot,[bestalpha[ip],bestalpha[ip]],!y.crange,color=fsc_color('red'),linestyle=2
         oplot,[bestalpha[0],bestalpha[0]],!y.crange,color=fsc_color('green'),linestyle=2
      endif
                                ;print out
      print, 'age(Gyr): ',10^(bestage-9.)
      print, '[Z/H]: ' ,bestz
      print, '[alpha/Fe]: ',bestalpha
      print, ';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
         ;;write to output structure
      outarr(goodobj(io)).chisq = chisqarr
      outarr(goodobj(io)).prob_age = [[reform(agethomas[*,0,0])],[probage]]
      outarr(goodobj(io)).prob_zh = [[reform(zthomas[0,*,0])],[probz]]
      outarr(goodobj(io)).prob_alpha = [[reform(alphathomas[0,0,*])],[probalpha]]      
      outarr(goodobj(io)).best_age = bestage
      outarr(goodobj(io)).best_zh = bestz
      outarr(goodobj(io)).best_alpha = bestalpha  
      outarr(goodobj(io)).keepindex = keepindex
   endif
endfor

outfile = outdir+outname+'_lickmeasurement.fits'
mwrfits,outarr,outfile,/create,/silent
stop
end
