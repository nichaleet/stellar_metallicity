;+
; NAME:
;       enk_crvlength
;
; PURPOSE:
;       This function computes the length of a curve with a tabular
;       representation, y(i) = F(x(i)).  The function need not be
;       monotonic, one-to-one, or onto.
;
; CATEGORY:
;       Numerical Analysis
;
; CALLING SEQUENCE:
;       Result = enk_crvlength(X, Y)
;
; INPUTS:
;       X:    An N-element vector (N >= 3) of type float or double. These 
;             values need not be specified in ascending order.
;
;       Y:    An N-element vector of type float or double.
;
; RESTRICTIONS:
;       Data that is highly oscillatory requires a sufficient number
;       of samples for an accurate curve length computation.
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, March 1996
;       modified to accept non-monotonic X by ENK 2006-11-03
;-

FUNCTION enk_crvlength, X, Y, print=print

  ON_ERROR, 2

  TypeX = SIZE(X)
  TypeY = SIZE(Y)

  ;Check Y data type.
  if TypeY[TypeY[0]+1] ne 4 and TypeY[TypeY[0]+1] ne 5 then $
    MESSAGE, "Y values must be float or double."

  ;Check length.
  ;if TypeX[TypeX[0]+2] lt 3 then $
  ;  MESSAGE, "X and Y arrays must contain 3 or more elements."
  if TypeX[TypeX[0]+2] ne TypeY[TypeY[0]+2] then $
    MESSAGE, "X and Y arrays must have the same number of elements."

  ;Check duplicate values.
  ;if TypeX[TypeX[0]+2] ne N_ELEMENTS(UNIQ(X)) then $
  ;  MESSAGE, "X array contains duplicate points."

  length = 0d
  step = -1d * ts_diff(X, 1)
  inc = where(step gt 0, cinc)
  dec = where(step lt 0, cdec)
  vert = where(step eq 0, cvert)
  cvert = cvert - 1

  if cinc eq 1 then length = length + sqrt((x[inc+1]-x[inc])^2d + (y[inc+1]-y[inc])^2d)
  if ~finite(length) then stop
  if cinc gt 1 then begin
      incstep = fix(-1*ts_diff(inc, 1))
      winc = where(incstep ne 1, cincstep)
      k = 0L
      for i=0,cincstep-1 do begin
          range = [inc[k:winc[i]], inc[winc[i]]+1]
          n = n_elements(range)
          if n eq 1 then message, 'Unknown error.'
          if n eq 2 then length = length + sqrt((x[range[1]]-x[range[0]])^2d + (y[range[1]]-y[range[0]])^2d)
          if n gt 2 then begin
              dydx = (shift(y[range],-1) - shift(y[range],1)) / (shift(x[range],-1) - shift(x[range],1) + 0.0d)
              dydx[0] = (-3.0d*y[range[0]] + 4.0d*y[range[1]] - y[range[2]]) / (x[range[2]] - x[range[0]])
              dydx[n-1] = (3.0d*y[range[n-1]] - 4.0d*y[range[n-2]] + y[range[n-3]]) / (x[range[n-1]] - x[range[n-3]])
              length = length + trap_int(x[range], sqrt(1d + dydx^2d))
          endif
          k = winc[i] + 1
      endfor
  endif
  if ~finite(length) then stop
  
  if cdec eq 1 then length = length + sqrt((x[dec+1]-x[dec])^2d + (y[dec+1]-y[dec])^2d)
  if ~finite(length) then stop
  if cdec gt 1 then begin
      decstep = fix(-1*ts_diff(dec, 1))
      wdec = where(decstep ne 1, cdecstep)
      k = 0L
      for i=0,cdecstep-1 do begin
          range = reverse([dec[k:wdec[i]], dec[wdec[i]]+1])
          n = n_elements(range)
          if n eq 1 then message, 'Unknown error.'
          if n eq 2 then length = length + sqrt((x[range[1]]-x[range[0]])^2d + (y[range[1]]-y[range[0]])^2d)
          if n gt 2 then begin
              dydx = (shift(y[range],-1) - shift(y[range],1)) / (shift(x[range],-1) - shift(x[range],1) + 0.0)
              dydx[0] = (-3.0*y[range[0]] + 4.0*y[range[1]] - y[range[2]]) / (x[range[2]] - x[range[0]])
              dydx[n-1] = (3.0*y[range[n-1]] - 4.0*y[range[n-2]] + y[range[n-3]]) / (x[range[n-1]] - x[range[n-3]])
              length = length + trap_int(x[range], sqrt(1d + dydx^2d))
          endif
          k = wdec[i] + 1
      endfor
  endif
  if ~finite(length) then stop

  if cvert eq 1 then length = length + abs(y[vert[1]] - y[vert[0]])
  if ~finite(length) then stop
  if cvert gt 1 then begin
      vert = vert[0:cvert-1]
      vertstep = fix(-1*ts_diff(vert, 1))
      wvert = where(vertstep ne 1, cvertstep)
      k = 0L
      for i=0,cvertstep-1 do begin
          range = vert[k:wvert[i]]
          n = n_elements(range)
          if n eq 1 then message, 'Unknown error.'
          for j=0,n-2 do length = length + abs(y[range[j+1]] - y[range[j]])
          k = wvert[i] + 1
      endfor
  endif
  if ~finite(length) then stop

  return, length
END
