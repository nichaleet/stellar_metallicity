function coelho_interp, feh, age, afe
  ;return response function at given z, age, and alpha
  common coelho_spec, c07, c07feh, c07age, c07afe

;;Get the whole set of spectra;;;;;;;;;;
  if (size(c07))[1] eq 0 then begin
     f = file_search('/scr2/nichal/workspace2/coelho/coelho07/coelho+07*.fits')
     c07 = coelho_read_spec(f)
  endif

  c07feh = c07[uniq(c07.feh,sort(c07.feh))].feh
  c07afe = c07[uniq(c07.afe,sort(c07.afe))].afe
  c07age = c07[uniq(c07.agegyr,sort(c07.agegyr))].agegyr
  
  nc07feh = n_elements(c07feh) 
  nc07age = n_elements(c07age)
  nc07afe = n_elements(c07afe)

  ;rearrange into a structure of 3 dimensions [feh,age,afe]
  str  = replicate(c07[0],nc07feh,nc07age,nc07afe)
  for i=0,nc07feh-1 do for j=0,nc07age-1 do for k=0,nc07afe-1 do begin
     wc07 = where(c07.feh eq c07feh[i] and c07.agegyr eq c07age[j] and c07.afe eq c07afe[k],cwc07)
     if cwc07 ne 1 then message,'you have problem with unique values'
     str[i,j,k] = c07[wc07]
  endfor
  c07 = str
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;Get an interpolated spectrum of the specified params

  ;if feh lt min(c07feh) or feh gt max(c07feh) then $
  ;   message, '[Fe/H] is off grid.',/continue
  ;if age lt min(c07age) or age gt max(c07age) then $
  ;   message, 'Age is off grid.',/continue
  ;if afe lt min(c07afe) or afe gt max(c07afe) then $
  ;   message, '[alpha/Fe] is off grid.',/continue
   
  wfeh = value_locate(c07feh,feh)
  wage = value_locate(c07age,age)
  wafe = value_locate(c07afe,afe)

  if wfeh eq -1 then ifeh = [0,1] else $
     if wfeh eq nc07feh-1 then ifeh =[nc07feh-2,nc07feh-1] else $
        ifeh = [wfeh,wfeh+1]
  
  if wage eq -1 then iage = [0,1] else $
     if wage eq nc07age-1 then iage =[nc07age-2,nc07age-1] else $
        iage = [wage,wage+1]

  if wafe eq -1 then iafe = [0,1] else $
     if wafe eq nc07afe-1 then iafe =[nc07afe-2,nc07afe-1] else $
        iafe = [wafe,wafe+1]

  ;trilinear interpolation
  dfeh = (feh-c07feh[ifeh[0]])/(c07feh[ifeh[1]]-c07feh[ifeh[0]])
  dage = (age-c07age[iage[0]])/(c07age[iage[1]]-c07age[iage[0]])
  dafe = (afe-c07afe[iafe[0]])/(c07afe[iafe[1]]-c07afe[iafe[0]])
  
 
  c00 = c07[ifeh[0],iage[0],iafe[0]].spec*(1.-dfeh)+c07[ifeh[1],iage[0],iafe[0]].spec*dfeh
  ;plot,c07[ifeh[0],iage[0],iafe[0]].spec  
  ;oplot,c07[ifeh[1],iage[0],iafe[0]].spec,color=fsc_color('green')
  ;oplot,c00,color=fsc_color('red')
  c01 = c07[ifeh[0],iage[0],iafe[1]].spec*(1.-dfeh)+c07[ifeh[1],iage[0],iafe[1]].spec*dfeh
  c10 = c07[ifeh[0],iage[1],iafe[0]].spec*(1.-dfeh)+c07[ifeh[1],iage[1],iafe[0]].spec*dfeh
  c11 = c07[ifeh[0],iage[1],iafe[1]].spec*(1.-dfeh)+c07[ifeh[1],iage[1],iafe[1]].spec*dfeh

  c0 = c00*(1.-dage)+c10*dage
  c1 = c01*(1.-dage)+c11*dage

  c07spec = c0*(1.-dafe)+c1*dafe

  return, {lambda:c07[ifeh[0],iage[0],iafe[1]].lambda,spec:c07spec}
end
