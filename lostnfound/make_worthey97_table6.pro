pro make_worthey97_table6
table = []
;age = 1.5
Fe=[-0.225,0.000,0.250,0.500]
HDA = [5.208,3.864,1.771,0.421] 
HGA = [3.850,2.120,-0.619,-2.303]
for i=0,3 do begin
   entry =create_struct('age',1.5,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 2.
Fe=[-0.225,0.000,0.250,0.500]
HDA = [3.747,1.959,-0.030,-1.376] 
HGA = [2.014,-0.340,-3.068,-4.625]
for i=0,3 do begin
   entry =create_struct('age',2.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 3.
Fe=[-0.225,0.000,0.250,0.500]
HDA = [1.995,0.314,-1.192,-2.481]
HGA = [-0.247,-2.639,-4.561,-5.937] 
for i=0,3 do begin
   entry =create_struct('age',3.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 5.
Fe=[-0.225,0.000,0.250,0.500]
HDA = [0.287,-0.632,-2.335,-3.691]
HGA = [-2.790,-3.953,-5.943,-7.318] 
for i=0,3 do begin
   entry =create_struct('age',5.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 8.
Fe=[-2.00,-1.500,-1.000,-0.500,-0.250,0.000,0.250,0.500]
HDA = [6.374,5.034,3.980,0.396,-0.498,-1.751,-3.135,-4.325]
HGA = [5.109,3.421,2.162,-2.828,-3.896,-5.342,-6.854,-7.960]
for i=0,7 do begin
   entry =create_struct('age',8.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 12.
Fe=[-2.00,-1.500,-1.000,-0.500,-0.250,0.000,0.250,0.500]
HDA = [5.205,4.056,3.422,-0.303,-1.428,-2.609,-4.093,-5.234]
HGA = [3.525,2.049,1.144,-3.842,-5.074,-6.323,-7.802,-8.758]
for i=0,7 do begin
   entry =create_struct('age',12.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

;age = 17.
Fe=[-2.00,-1.500,-1.000,-0.500,-0.250,0.000,0.250,0.500]
HDA = [4.443,3.581,2.978,-0.974,-2.116,-3.354,-4.733,-5.704]
HGA = [2.509,1.254,0.426,-4.762,-5.893,-7.091,-8.409,-9.162]
for i=0,7 do begin
   entry =create_struct('age',17.0,'FeH',Fe(i),'HDA',HDA(i),'HGA',HGA(i))
   table = [table,entry]
endfor

worthey_table = table
save,worthey_table,filename='wortheytable.sav'
end
