pro run_spsfit2fast_gallazzi1
  for i=0,19 do begin
     ni=10*i
     nf=(10*(i+1))-1
     spsfit2fast,'gallazzi1',nbegin=ni,nfinal=nf
  endfor
end
