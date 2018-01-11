;+
; NAME: deimos_pix2arcsec
;
; PURPOSE: Convert pixel position on a DEIMOS CCD to offset in arcsec
;      from the center of the mosaic.
;
; CATEGORY: DEIMOS astrometry
;
; CALLING SEQUENCE:
;      deimos_pix2arcsec, ccd, xpix, ypix, xarcsec, yarcsec
;
;
; INPUTS:
;      ccd: CCD number, 1 through 8
;      xpix: x pixel, which can be fractional
;      ypix: y pixel, which can be fractional
;
; OPTIONAL INPUTS:
;      none
;
; KEYWORD PARAMETERS:
;      none
;
; OUTPUTS:
;      xarcsec: offset in arcsec along the x axis from the center of
;               the mosaic
;      yarcsec: offset in arcsec along the y axis from the center of
;               the mosaic
;
; OPTIONAL OUTPUTS:
;      none
;
; MODIFICATION HISTORY:
;      2008jun12: created - ENK
;-


pro deimos_pix2arcsec, ccd, xpix, ypix, xarcsec, yarcsec
    xnom = [-3203.8, -1067.6, 1067.6, 3203.8, -3203.8, -1067.6, 1067.6, 3203.8]
    ynom = [-2056., -2056., -2056., -2056., 2056., 2056., 2056., 2056.]
    delx = [-20.05, -12.64, 0.0, -1.34, -19.02, -9.65, 1.88, 4.81]
    dely = [14.12, 7.25, 0.0, -19.92, 16.46, 8.95, 1.02, -24.01]
    theta = [-0.082, -0.030, 0.0, -0.1206, 0.136, -0.060, -0.019, -0.082]
    theta *= !DTOR

    platescale = 0.1185 ;arcsec/pix

    ccdtype = size(ccd, /type)
    if ccdtype lt 1 or (ccdtype gt 3 and ccdtype lt 12) or ccdtype gt 15 then message, 'CCD must be an integer.'
    if ccd lt 1 or ccd gt 8 then message, 'CCD must be between 1 and 8.'
    ccd -= 1

    ;if xpix lt 0.0 or xpix gt 2048.0 then message, 'xpix is out of range.'
    ;if ypix lt 0.0 or ypix gt 4096.0 then message, 'ypix is out of range.'

    ;if ccd le 4 then xpix = 2048.0 - xpix
    ;if ccd mod 2 eq 0 then ypix = 4096.0 - ypix

    xpix -= 1024.0
    ypix -= 2048.0

    x = xpix*cos(theta[ccd]) - ypix*sin(theta[ccd]) + xnom[ccd] + delx[ccd]
    y = xpix*sin(theta[ccd]) + ypix*cos(theta[ccd]) + ynom[ccd] + dely[ccd]

    xarcsec = x*platescale
    yarcsec = y*platescale

    return
end
