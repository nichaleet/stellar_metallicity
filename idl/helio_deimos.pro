;+
; NAME:
;   lco_helio
;
; PURPOSE:
;   Compute correction term to add to velocities to convert to heliocentric.
;
; CALLING SEQUENCE:
;   vcorr = lco_helio( ra, dec, [ epoch, jd=, tai=, $
;    longitude=, latitude=, altitude= ] )
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   epoch          - Epoch of observation for RA, DEC; default to 2000.
;
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.  Note this should probably be
;                    type DOUBLE.
;   tai            - Number of seconds since Nov 17 1858; either JD or TAI
;                    must be specified.  Note this should probably either
;                    be type DOUBLE or LONG64.
;   longitude      - Longitude of observatory;
;                    default to LCO
;   latitute       - Latitude of observatory; default to LCO
;   altitude       - Altitude of observatory; default to LCO
;                    Keck
;
; OUTPUTS:
;   vcorr          - Velocity correction term, in km/s, to add to measured
;                    radial velocity to convert it to the heliocentric frame.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   baryvel
;   ct2lst
;
; REVISION HISTORY:
;   09-May-2000  Written by S. Burles & D. Schlegel
;   30-Aug-2002  Revised by JXP
;   05-Feb-2005  Revised by MG
;-
;------------------------------------------------------------------------------

function helio_deimos, ra, dec, epoch, jd=jd, tai=tai, $
 longitude=longitude, latitude=latitude, altitude=altitude 

   if (NOT keyword_set(epoch)) then epoch = 2000.0

   ; Default to location of KECK
   if (NOT keyword_set(longitude)) then longitude = 360. -155.478
   if (NOT keyword_set(latitude)) then latitude = 19.8283
   if (NOT keyword_set(altitude)) then altitude = 4160.

   vcor = heliocentric(ra,dec,2000, longitude=longitude, latitude=latitude, $
                       altitude=altitude,jd=jd)

   return, vcor
end

