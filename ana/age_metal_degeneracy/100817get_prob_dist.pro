pro get_prob_dist,filein=filein,fileout=fileout
;input is cube of chisq must be in dimension of [feh,age,galaxy]
;convert chisq cube into probability distributions (normalized)
    chisqcube = mrdfits(filein,0)
    fehgrid = mrdfits(filein,1)
    agegrid = mrdfits(filein,2)
    objname = mrdfits(filein,3)
    feharr = fehgrid[*,0]
    agearr = reform(agegrid[0,*]) 
    cube_dimen = size(chisqcube,/dimensions)
    ngals = cube_dimen(2)
    
    deltafeh = fehgrid[1,0]-fehgrid[0,0]
    deltaage = agegrid[0,1]-agegrid[0,0]
    probfeh_age_str = {objname:'name',feh:feharr,probfeh:fltarr(cube_dimen(0)),age:agearr,probage:fltarr(cube_dimen(1)),probfehage:fltarr(cube_dimen(0),cube_dimen(1))}
    probfeh_age_all = replicate(probfeh_age_str,ngals)

    for gg=0,ngals-1 do begin
       Lgrid = -0.5*(chisqcube[*,*,gg]-min(chisqcube[*,*,gg]))
       probgrid = exp(double(Lgrid))
       volume = 0
       for ii=0,cube_dimen(0)-1 do begin
          for jj=0,cube_dimen(1)-1 do begin
             volume = volume+deltafeh*deltaage*probgrid[ii,jj]
          endfor
       endfor
       probgrid = probgrid/volume
       probfeh_age_all(gg).probfehage = probgrid

       probfeh = fltarr(cube_dimen(0))
       probage = fltarr(cube_dimen(1))
       for ii=0,cube_dimen(0)-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
       for jj=0,cube_dimen(1)-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
       probfeh = probfeh/int_tabulated(feharr,probfeh)
       probage = probage/int_tabulated(agearr,probage)
       probfeh_age_all(gg).objname = objname.objname(gg)
       probfeh_age_all(gg).probfeh = probfeh
       probfeh_age_all(gg).probage = probage
    endfor
    mwrfits,probfeh_age_all,fileout,/create,/silent
end
