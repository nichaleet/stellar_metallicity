pro runmanyiter_sn_mock_deimos_spec

dirname='mockv3_240717_'+strtrim(string(indgen(20)),2)
for i=0,19 do make_mock_deimos_spec,dirname(i)

end
