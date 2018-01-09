pro chisq_contribution, lickoutfile
;input: lickoutfile --> output file from lick_measurement.pro
;purpose: calculate the contribution to total chisq from each index
;example: chisq_contribution,'all_cl0024_gauss_nomaxflux_emisssub_lickmeasurement.fits'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

lickfilename = (strsplit(lickoutfile,'.',/extract))[0]

;read input file
inputdir = '/scr2/nichal/workspace2/lickindices/output/results/'
lickstr = mrdfits(inputdir+lickoutfile,1)
goodobj = where(lickstr.good eq 1, cgoodobj)
ngals = n_elements(goodobj)
;output directory
outdir = '/scr2/nichal/workspace2/lickindices/output/chisq_contribution/'

;get the thomas11b prediction grid 
res=0.025
filename = '/scr2/nichal/workspace2/lickindices/thomas_fine'+strtrim(string(res,format='(F5.3)'),2)+'.sav'
filethomas = file_search(filename,count=fcount)
restore,filethomas[0]
agearr = reform(thomasgrid[*,0,0].age)
zarr = reform(thomasgrid[0,*,0].zh)
alphaarr = reform(thomasgrid[0,0,*].alphafe)
tagthomas = tag_names(thomasgrid)
tagmatch  = [where(strmatch(tagthomas,'cn1',/fold_case) eq 1),where(strmatch(tagthomas,'cn2',/fold_case) eq 1),where(strmatch(tagthomas,'ca4227',/fold_case) eq 1),where(strmatch(tagthomas,'g4300',/fold_case) eq 1),where(strmatch(tagthomas,'c24668',/fold_case) eq 1),where(strmatch(tagthomas,'mg1',/fold_case) eq 1),where(strmatch(tagthomas,'mg2',/fold_case) eq 1),where(strmatch(tagthomas,'mgb',/fold_case) eq 1),where(strmatch(tagthomas,'fe4383',/fold_case) eq 1),where(strmatch(tagthomas,'fe4531',/fold_case) eq 1),where(strmatch(tagthomas,'fe5270',/fold_case) eq 1),where(strmatch(tagthomas,'fe5335',/fold_case) eq 1),where(strmatch(tagthomas,'fe5406',/fold_case) eq 1),where(strmatch(tagthomas,'hda',/fold_case) eq 1),where(strmatch(tagthomas,'hga',/fold_case) eq 1),where(strmatch(tagthomas,'hgf',/fold_case) eq 1)]
tagthomasout = tagthomas(tagmatch)
if n_elements(tagmatch) ne 16 then stop, 'Oh no'

;check if the input file assumed [alpha/FE] =0
chisqsize = size(lickstr[0].chisq,/dimensions)
if n_elements(chisqsize) eq 3 then alphafix = 0
if n_elements(chisqsize) eq 2 then alphafix = 1

set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
!p.charsize=1

;begin loop over each galaxy
for i=0,ngals-1 do begin
   obj = lickstr(goodobj(i))   

   ;find where the best values are
   minagediff = min(agearr-obj.best_age[0],ageloc)
   minzdiff = min(zarr-obj.best_zh[0],zloc)
   if alphafix eq 0 then $
      minalphadiff = min(alphaarr-obj.best_alpha[0],alphaloc)
   if alphafix eq 1 then $
      minalphadiff = min(alphaarr-0.,alphaloc)
   
   index = [obj.cn1[0],obj.cn2[0],obj.ca4227[0],obj.g4300[0],obj.c24668[0],obj.mg1[0],obj.mg2[0],obj.mgb[0],obj.fe4383[0],obj.fe4531[0],obj.fe5270[0],obj.fe5335[0],obj.fe5406[0],obj.hda[0],obj.hga[0],obj.hgf[0]]
   index_err = [obj.cn1[1],obj.cn2[1],obj.ca4227[1],obj.g4300[1],obj.c24668[1],obj.mg1[1],obj.mg2[1],obj.mgb[1],obj.fe4383[1],obj.fe4531[1],obj.fe5270[1],obj.fe5335[1],obj.fe5406[1],obj.hda[1],obj.hga[1],obj.hgf[1]]
   goodind = where(index ne -99 and index_err ne -99,cgoodind) ;goodind of observed spec
   chiarr = fltarr(16)-99.
   if cgoodind gt 0 then begin
      index = index(goodind)
      index_err = index_err(goodind)
      tagmatchnow = tagmatch(goodind)
      for ii=0,cgoodind-1 do begin ;loop over each index
         chiarr(goodind(ii)) = (index(ii)-thomasgrid[ageloc,zloc,alphaloc].(tagmatchnow(ii)))^2/(index_err(ii))^2
      endfor   
      psname=outdir+lickfilename+'_chisqcontribution_'+strtrim(string(goodobj(i)),2)+'.eps'
      device, filename = psname,xsize = 30,ysize = 10, $
              xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      vsym,4,/fill,rot=45
      xval = findgen(16)
      plot,xval,chiarr,min_value=min(chiarr(goodind))-5.,xtickv=xval,xticks=15,xtickname=tagthomasout,ytitle='chi square',psym=8,xrange=[-0.5,15.5],xstyle=1
      device,/close
   endif
   
   ;plot
   
endfor
 
stop
end
