pro find_id_from_pos,ra_pos,dec_pos,data,raind=raind,decind=decind,idind=idind,id,distance
id = indgen(n_elements(ra_pos))
distance = findgen(n_elements(ra_pos))
for ii=0,n_elements(ra_pos)-1 do begin
   ra = ra_pos(ii)
   dec = dec_pos(ii)
   diff_ra = data(raind,*)-ra
   diff_dec = data(decind,*)-dec
   diff_dis = sqrt(diff_ra^2+diff_dec^2)
   position = where(diff_dis eq min(diff_dis))
   if position(0) eq -1 then stop
   id(ii) = data(idind,position(0))
   distance(ii) = diff_dis(0)
endfor
end
