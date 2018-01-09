pro oplotlinenames
       linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
      linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', ' ', ' ', 'H8', 'CaK', ' ', 'CaH', ' ', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', ' ', 'Mgb', ' ', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
      linecolors = ['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue']
      ;mark and label lines
      n = n_elements(linewaves)
      for j=0,n-1 do begin
         if (linewaves)[j] le !X.CRANGE[0] or (linewaves)[j] ge !X.CRANGE[1] then continue
         oplot, [(linewaves)[j], (linewaves)[j]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[j])
         xyouts, (linewaves)[j]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[j], orientation=90, alignment=1, color=fsc_color((linecolors)[j]),charsize=0.8
      endfor


end
