;+
; NAME:
;   deimos_traceall
;
; PURPOSE:
;   find traces of slitlets in a chip
;
; CALLING SEQUENCE:
;   deimos_traceall, chipno, flatimage, flativar, $
;         slitcoords, nslitlets, xpos1, xpos2, ypos1, ypos2, indblu, $
;         ybin=ybin, plot=plot, pspath=pspath  ,badmatch=badmatch, $
;         model=model
;
; INPUTS:
;    chipno -- which of the 8 chips (1-8)
;    flatimage -- spectral flat for this mask
;    flativar  -- invvar of flat
;    slitcoords -- structure detailing bluslit information
;    
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   ybin -- (default 8) binning of raw spectrum for traceset
;   plot  -- set to generate test plots on screen
;   pspath -- if set, path for ps files 
;   model -- optical model structure
;
; OUTPUTS:
;   nslitlets -- how many slitlets on this chip
;   xpos1 -- xposition of start edge of traces
;   xpos2 -- xposition of end edge of traces
;   ypos1 -- yposition of start edge of traces
;   ypos2 -- yposition of end edge of traces
;   synth1 -- array[nslit] =1 if xpos1 was synthesized from neighbors
;   synth2 -- array[nslit] =1 if xpos2 was synthesized from neighbors
;   indblu -- index of structure detailing slitlets, from bluprint FITS header
;   badmatch - set to 1 on output if the slit tracing was problematic
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   split off from deimos_spslit because this code will be called from both
;    quick-look and final spec2d analysis
;
; REVISION HISTORY:
;
;       Wed Feb 20 17:52:02 2002, Douglas Finkbeiner (dfink)
;             Split of from 2dtest, from 2001-Aug-22
;       2002-apr-13  MD - split from deimos_spslit 
;       2002-apr-25  DPF - cleaned up a bit more
;       2002-oct-14  DPF - removed flipsign and other deadwood
;
;----------------------------------------------------------------------
pro deimos_traceall, chipno, flatimage, flativar, $
          slitcoords, nslitlets, xpos1, xpos2, ypos1, ypos2, $
          synth1, synth2, indblu, $
          ybin=ybin, plot=plot, pspath=pspath, badmatch=badmatch  , $
                     model=model
  
; -------- Set defaults:
  if NOT keyword_set(ybin) then ybin = 8
  if NOT keyword_set(flativar) then begin 
     print, 'DEIMOS_TRACEALL:  You should really set flativar!'
  endif else ivar = flativar NE 0  ; treat as mask (1=good) for now


; determine gap between slitlets for smoothing later
  diffs=model.xb-shift(model.xt,1)
  gooddiff=where(lindgen(n_elements(diffs)) NE 0 AND diffs ne 0,ct)
  if ct gt 2 then gap=median(diffs[gooddiff],/even) else gap=5.
  gap = gap > 3
  gap = gap < 45

; -------- Find some traces - ncoeff=4 means third order fit 
  deimos_trace_crude, flatimage, tset1, tset2, flatbin, ybin=ybin, ncoeff=4, $
    invvar=ivar,gap=gap

; -------- reject traces with strange (don't match neighbors) coefficients
  tset1dim = size(tset1.coeff, /dimens)                                                                           ;ENK
  tset2dim = size(tset2.coeff, /dimens)                                                                           ;ENK
  if n_elements(tset1dim) eq 1 or n_elements(tset2dim) eq 1 then ntset = 1 else ntset=tset1dim[1] < tset2dim[1]   ;ENK
  ;ntset=(size(tset1.coeff, /dimens))[1] < (size(tset2.coeff, /dimens))[1]

  if ntset gt 5 then begin
      clean_tset, tset1
      clean_tset, tset2
  endif

; -------- compute tracesets for top and bottom of slits
  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2

;; -------- ENK: what if xpos1 > xpos2?
;  cbrokenslit = 1
;  halfy1 = floor((size(xpos1, /dimens))[0] / 2.)
;  halfy2 = floor((size(xpos2, /dimens))[0] / 2.)
;  while cbrokenslit gt 0 do begin
;      wbrokenslit = where(xpos1[halfy1,*] ge xpos2[halfy2,*], cbrokenslit)
;      if cbrokenslit gt 0 then begin
;          do1 = n_elements(xpos1[halfy1,*]) le n_elements(xpos2[halfy2,*])
;          do2 = n_elements(xpos1[halfy1,*]) ge n_elements(xpos2[halfy2,*])
;          if do1 then begin
;              startfix = min(wbrokenslit)
;              fix2 = complement([startfix], n_elements(xpos2[halfy2,*]))
;              xpos2 = xpos2[*,fix2]
;              tset2 = {func:tset2.func, xmin:tset2.xmin, xmax:tset2.xmax, coeff:tset2.coeff[*,fix2]}
;          endif
;          if do2 then begin
;              endfix = max(wbrokenslit)
;              fix1 = complement([endfix], n_elements(xpos1[halfy1,*]))
;              xpos1 = xpos1[*,fix1]
;              tset1 = {func:tset1.func, xmin:tset1.xmin, xmax:tset1.xmax, coeff:tset1.coeff[*,fix1]}
;          endif
;      endif
;  endwhile
;
;; -------- ENK: remove duplicate traces matched to the same slit
;  mintracewidth = 2  ;pixels
;
;  xpos1dim = (size(xpos1, /dimens))                                                  ;ENK
;  if n_elements(xpos1dim) eq 1 then n1 = 1 else n1 = xpos1dim[1]                     ;ENK
;  ;n1 = (size(xpos1, /dimens))[1]
;  if n1 gt 1 then begin
;      goodtrace = bytarr(n1)+1B
;      for i=1,n1-1 do if xpos1[halfy1,i] lt xpos1[halfy1,i-1]+mintracewidth then goodtrace[i-1]=0
;      wgoodtrace = where(goodtrace)
;      tset1 = {func:tset1.func, xmin:tset1.xmin, xmax:tset1.xmax, coeff:tset1.coeff[*,wgoodtrace]}
;      xpos1=xpos1[*,wgoodtrace]
;      ypos1=ypos1[*,wgoodtrace]
;  endif
;
;  xpos2dim = (size(xpos2, /dimens))                                                  ;ENK
;  if n_elements(xpos2dim) eq 1 then n2 = 1 else n2 = xpos2dim[1]                     ;ENK
;  ;n2 = (size(xpos2, /dimens))[1]
;  if n2 gt 1 then begin
;      goodtrace = bytarr(n2)+1B
;      for i=1,n2-1 do if xpos2[halfy2,i] lt xpos2[halfy2,i-1]+mintracewidth then goodtrace[i-1]=0
;      wgoodtrace = where(goodtrace)
;      tset2 = {func:tset2.func, xmin:tset2.xmin, xmax:tset2.xmax, coeff:tset2.coeff[*,wgoodtrace]}
;      xpos2=xpos2[*,wgoodtrace]
;      ypos2=ypos2[*,wgoodtrace]
;  endif
;
;  xpos1dim = (size(xpos1, /dimens))                                                  ;ENK
;  if n_elements(xpos1dim) eq 1 then nslitlets = 1 else nslitlets = xpos1dim[1]       ;ENK
;  ;nslitlets = (size(xpos1, /dimens))[1]
;
;; -------- ENK: remove overlapping slits
;  i = 0
;  repeat begin
;      if xpos2[halfy2,i] gt xpos1[halfy1,i+1] then begin
;          w1 = complement(i, nslitlets, c)
;          tset1 = {func:tset1.func, xmin:tset1.xmin, xmax:tset1.xmax, coeff:tset1.coeff[*,w1]}
;          xpos1=xpos1[*,w1]
;          ypos1=ypos1[*,w1]
;          w2 = complement(nslitlets-1, nslitlets, c)
;          tset2 = {func:tset2.func, xmin:tset2.xmin, xmax:tset2.xmax, coeff:tset2.coeff[*,w2]}
;          xpos2=xpos2[*,w2]
;          ypos2=ypos2[*,w2]
;          nslitlets--
;      endif
;      i++
;  endrep until i ge nslitlets-1


; -------- debugging plots
  ny = 4096
  yind = lindgen(ny/ybin)*ybin
  if keyword_set(plot) then begin 
     atv, flatbin, min=0, max=60000
     atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=4
     atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=6
  endif 

; -------- postscript plots (for paper)
  if nslitlets gt 1 and keyword_set(pspath) then begin 
     dfpsplot, pspath+'flat1.ps', /sq, /color, bits=8
     imps = flatbin[0:ny/ybin-1, *]
     display,bytscl(imps,min=-100,max=65000), xtit='spatial [pix]', ytit=$
       'lambda [pix]', chars=1.5, xmargin=[7, 2]
     ;nline1 = (size(xpos1))[2]
     xpos1dim = size(xpos1)                                                             ;ENK
     nline1 = xpos1dim[0]                                                               ;ENK
     ;nline2 = (size(xpos2))[2]
     xpos2dim = size(xpos2)                                                             ;ENK
     nline2 = xpos2dim[0]                                                               ;ENK
     if nline1 eq 1 then begin                                                          ;ENK
         for i=0, nline1-1 do oplot, xpos1[yind], ypos1[yind]/ybin,color=4              ;ENK
     endif else begin                                                                   ;ENK
         for i=0, nline1-1 do oplot, xpos1[yind, i], ypos1[yind, i]/ybin,color=4
     endelse                                                                            ;ENK
     if nline2 eq 1 then begin                                                          ;ENK
         for i=0, nline2-1 do oplot, xpos2[yind], ypos2[yind]/ybin,color=6              ;ENK
     endif else begin                                                                   ;ENK
         for i=0, nline2-1 do oplot, xpos2[yind, i], ypos2[yind, i]/ybin,color=6
     endelse                                                                            ;ENK
     dfpsclose
  endif 

; -------- use information from slitcoords (from blueprint) to match
;           found slit positions with slit numbers.             
;          This overwrites tset1, tset2 to include anything we missed. 

  deimos_slit_match, chipno, tset1, tset2, indblu, $
    synth1, synth2, slitcoords=slitcoords, badmatch=badmatch, model=model

; indblu is index in bluslit (or slitcoords) for each slit found

; -------- compute tracesets for top and bottom of slits (again)
  ;xpos1old = xpos1;;
  ;xpos2old = xpos2;;
  delvarx, ypos1
  delvarx, ypos2
  traceset2xy, tset1, ypos1, xpos1
  traceset2xy, tset2, ypos2, xpos2
  ;stop                                                                         ;;
  xpos1dim = (size(xpos1, /dimens))                                                  ;ENK
  if n_elements(xpos1dim) eq 1 then nslitlets = 1 else nslitlets = xpos1dim[1]       ;ENK
  ;nslitlets = (size(xpos1, /dimens))[1]

  if keyword_set(plot) then begin 
     yind = lindgen(ny/ybin)*ybin
     atvplot, xpos1[yind, *], ypos1[yind, *]/ybin, ps=3,color=5
     atvplot, xpos2[yind, *], ypos2[yind, *]/ybin, ps=3,color=7
     ws1 = where(synth1, nws1)
     if nws1 NE 0 then atvplot, xpos1[*, ws1], ypos1[*, ws1]/ybin, ps=3,color=1
     ws2 = where(synth2, nws2)
     if nws2 NE 0 then atvplot, xpos2[*, ws2], ypos2[*, ws2]/ybin, ps=3,color=2
     abandy = 300-(slitcoords[indblu].ymm*.87)
     for i=0, nslitlets-1 do begin 
        ii = i*12+100
        yband = abandy[i]*ybin
        atvplot, [xpos1[yband, i], xpos2[yband, i]], $
          [1, 1]*abandy[i], thick=5, color=3
        atvxyouts, 0.5*(xpos1[yband, i]+xpos2[yband, i]), abandy[i]+8, $
          string(indblu[i], format='(I3.3)'), align=0.5, charsize=1.8, $
          color=2
     endfor
  endif 

  return
end
