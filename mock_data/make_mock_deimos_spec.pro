pro make_mock_deimos_spec,dirname,notell=notell,notput=notput,alpha=alpha
  ;initial set up
  redshift = 0.55
  z = -0.1
  age = 3.0
  vdisp = 250.  ;FWHM!!
 
  sample1 = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/nicha_reduced/smm1a/spec1d.smm1a.010.N31733.fits.gz',1,hdrblue)
  sample2 = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/nicha_reduced/smm1a/spec1d.smm1a.010.N31733.fits.gz',2,hdrred)
  lambdablue = sample1.lambda
  lambdared  = sample2.lambda
  ;dlambda
  restore, '/scr2/nichal/workspace2/sps_fit/data/smm1a/specres_poly.sav'
  dlamblue = poly(lambdablue/1000.-7.8, specres_poly)/2.35
  dlamred  = poly(lambdared/1000.-7.8, specres_poly)/2.35
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;make SPS spectrum
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  spsstruct = sps_interp(z, age)
  spslambda = spsstruct.lambda
  spsspec = spsstruct.spec
  w = where(spslambda gt 0 and spslambda lt 100000, c)
  if c lt 25 then message, 'Not enough pixels.'
  spslambda = spslambda[w]
  spsspec = spsspec[w]
  clight = 299792.458
  spsspec = spsspec*clight/spslambda^2      ;change fnu(Lsun/Hz) to flambda 
  spsmedian = median(spsspec)
  spsspec = spsspec/spsmedian      ;normalize to around 1
  ;smooth to vdisp and deimos resolution
  spsspec = smooth_gauss_wrapper(spslambda, spsspec, spslambda, vdisp/clight/2.35*spslambda)
  spsspecblue = smooth_gauss_wrapper(spslambda*(1.+redshift), spsspec, lambdablue, dlamblue)
  spsspecred = smooth_gauss_wrapper(spslambda*(1.+redshift), spsspec, lambdared, dlamred)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;add ALPHA enhancement at Mg region
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if keyword_Set(alpha) then begin
     z_alpha = alpha
     spsstruct_a = sps_interp(z_alpha,age)
     spslambda_a = spsstruct_a.lambda
     spsspec_a   = spsstruct_a.spec
     
     w = where(spslambda_a gt 0 and spslambda_a lt 100000, c)
     if c lt 25 then message, 'Not enough pixels.'
     spslambda_a = spslambda_a[w]
     spsspec_a = spsspec_a[w]
     spsspec_a = spsspec_a*clight/spslambda_a^2 ;change fnu(Lsun/Hz) to flambda 
     spsspec_a = spsspec_a/spsmedian    ;normalize to the same level as the main
     ;;smooth to vdisp and deimos resolution
     spsspec_a = smooth_gauss_wrapper(spslambda_a, spsspec_a, spslambda_a, vdisp/clight/2.35*spslambda_a)
     spsspecblue_a = smooth_gauss_wrapper(spslambda_a*(1.+redshift), spsspec_a, lambdablue, dlamblue)
     spsspecred_a = smooth_gauss_wrapper(spslambda_a*(1.+redshift), spsspec_a, lambdared, dlamred)
     
     ;;Continuum normalize both spectra
     mask = bytarr(n_Elements(lambdared))
     lambdarest = lambdared/(1.+redshift)
     readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indstart,indend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,format='I,D,D,D,D,D,D,I,A',comment='#',/silent
     use_indices = ['Mg_b','Mg_2','Mg_1']
     for i=0,n_elements(use_indices)-1 do begin
        indnow = where(indname eq use_indices(i),cindnow)
        if cindnow eq 0 then stop
        w = where(lambdarest gt indstart(indnow[0]) and lambdarest lt indend(indnow[0]),cw)
        if cw gt 0 then mask[w]=1
     endfor
     won = where(mask eq 0,con)
     
     bkpt = slatec_splinefit(lambdared[won], spsspecred[won], coeff, bkspace=150, upper=3, lower=3, /silent)
     if bkpt[0] eq -1 then return
     cont_main = slatec_bvalu(lambdared, bkpt, coeff)
     
     bkpt = slatec_splinefit(lambdared[won], spsspecred_a[won], coeff, bkspace=150, upper=3, lower=3, /silent)
     cont_a = slatec_bvalu(lambdared,bkpt,coeff)

     spsspecred_contdiv = spsspecred/cont_main
     spsspecred_a_contdiv = spsspecred_a/cont_a

     ;;Insert the alpha enhanced spectrum to the main spectrum in the alpha sensitive region
     ;;The range of Mg1,Mg2,Mgb are from 4926 to 5333A which is only in the red spectra
     newspec = spsspecred_contdiv
     won = where(mask eq 1, con)
     if con gt 0 then newspec[won] = spsspecred_a_contdiv[won] else stop
     
     set_plot,'x'
     !p.multi=[0,1,2]
     plot,lambdarest,spsspecred_contdiv,xrange=[4895,5366]
     oplot,lambdarest,spsspecred_a_contdiv,color=fsc_color('red')
     oplot,lambdarest[won],newspec[won],color=fsc_color('green')
     al_legend,['main','alpha','final'],color=fsc_color(['white','red','green']),psym=15
     
     plot,lambdarest,spsspecred
     oplot,lambdarest,cont_main
     oplot,lambdarest,spsspecred_a,color=fsc_color('red')
     oplot,lambdarest,cont_a,color=fsc_color('red')
     oplot,lambdarest[won],newspec[won]*cont_main[won],color=fsc_color('green')
     ;;apply back the continuum
     spsspecred = newspec*cont_main
  endif
 ; stop
  ;;;;;;;;;get telluric features;;;;;;;;;;;;;;;;;;;
  tell = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent)

  airmass = sxpar(hdrblue, 'AIRMASS')
  aratio = airmass/tell[0].airmass
  telllambda = tell[0].lambda
  tellspec = (tell[0].spec)^aratio
  tellivar = ((tell)[0].ivar)*(((tell)[0].spec)/(aratio*tellspec))^2.
  ivarmissing = 10d10
  w = where(tellivar ge 1d8, c)
  if c gt 0 then tellivar[w] = ivarmissing
  f = where(finite(tellspec) and tellspec gt 0 and finite(tellivar) and tellivar gt 0 and tellivar lt 1d8)
  telllambda = telllambda[f]
  tellspec = tellspec[f]
  tellivarfull = tellivar
  tellivar = tellivar[f]

  tellspecblue= interpolate(tellspec, findex(telllambda, lambdablue), missing=1.)
  tellivarblue = interpolate(tellivarfull, findex(telllambda, lambdablue), missing=ivarmissing)
  tellspecred= interpolate(tellspec, findex(telllambda, lambdared), missing=1.)
  tellivarred = interpolate(tellivarfull, findex(telllambda, lambdared), missing=ivarmissing)

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;generate different signal to noise sets of data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;targetSN = 10.^(findgen(25)/30.+0.7)   
  targetSN = 10.^(findgen(20)/18.+0.4)
  npixblue = n_elements(lambdablue)
  npixred = n_elements(lambdared)
;  for i=0,19 do begin
  for i=0,n_elements(targetSN)-1 do begin
     ;add noise
     sn = targetSN[i]
     ;sn = 8.
     sigmasqblue = spsspecblue^2/sn^2
     sigmasqred  = spsspecred^2/sn^2
     specblue = spsspecblue+randomn(seed,npixblue)*sqrt(sigmasqblue)
     specred  = spsspecred+randomn(seed,npixred)*sqrt(sigmasqred)
     ivarblue = 1./sigmasqblue
     ivarred  = 1./sigmasqred
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;add telluric
     if ~keyword_set(notell) then begin
        specblue_tell = specblue*tellspecblue
        specred_tell  = specred*tellspecred
        ivarblue = 1./(specblue^2*(sigmasqblue/specblue^2+1./(tellivarblue*tellspecblue^2)))
        ivarred  = 1./(specred^2*(sigmasqred/specred^2+1./(tellivarred*tellspecred^2)))
     endif else begin
        specblue_tell = specblue
        specred_tell = specred
     endelse
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;;add DEIMOS instrument throughput
     if ~keyword_set(notput) then begin
        readcol,'thr_g900_70_gg455.asc',wltput,tput
        tputred  = interpol(tput,wltput,lambdared)
        tputblue = interpol(tput,wltput,lambdablue)
        tputred = tputred/median(tputred)
        tputblue = tputblue/median(tputblue)
        specblue_tell = specblue_tell*tputblue
        specred_tell = specred_tell*tputred
        ivarblue = ivarblue/tputblue^2
        ivarred = ivarred/tputred^2
     endif
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     bluestruct = {spec:specblue_tell,lambda:lambdablue,ivar:ivarblue,skyspec:sample1.skyspec,FeH:z,age:age,vdisp:vdisp,redshift:redshift,sn:sn}
     redstruct = {spec:specred_tell,lambda:lambdared,ivar:ivarred,skyspec:sample2.skyspec,FeH:z,age:age,vdisp:vdisp,redshift:redshift,sn:sn}
     if i lt 10 then prefix='/scr2/nichal/workspace2/mock_data/'+dirname+'/spec1d.'+dirname+'.00'
     if i ge 10 then prefix='/scr2/nichal/workspace2/mock_data/'+dirname+'/spec1d.'+dirname+'.0'
     if i gt 99 then prefix='/scr2/nichal/workspace2/mock_data/'+dirname+'/spec1d.'+dirname+'.'
     if i eq 0 then begin
        testdir = file_search('/scr2/nichal/workspace2/mock_data/'+dirname,/mark_directory,count=ctestdir)
	if ctestdir eq 0 then file_mkdir,'/scr2/nichal/workspace2/mock_data/'+dirname 
     endif
     mwrfits,bluestruct,prefix+strtrim(string(fix(i)),2)+'.mockspec.fits',hdrblue,/create
     mwrfits,redstruct,prefix+strtrim(string(fix(i)),2)+'.mockspec.fits',hdrred,/silent
  endfor
  if keyword_Set(alpha) then alphasave=alpha else alphasave=0
  params_save = {redshift:redshift,age:age,z:z,vdisp:vdisp,notell:keyword_set(notell),alpha:alphasave,notput:keyword_set(notput),SN:targetSN}
  save,params_save,filename='/scr2/nichal/workspace2/mock_data/'+dirname+'/variables.sav'
end
