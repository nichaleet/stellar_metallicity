pro make_cat_reg,txt_cat,name,lineskip=lineskip
if ~keyword_set(lineskip) then lineskip=146
data=read_ascii(txt_cat,data_start=lineskip)
data=data.field001
save,data,filename='clash_'+name+'_cat.sav'
ids=string(long(data[0,*]))
ra = data[1,*]
dec= data[2,*]
g=data[43,*]
i=data[67,*]
dim=size(data,/dimensions)
z = data[dim(0)-2,*]
good = where(z gt 0.6 and i lt 21)

;make string
g_i=strtrim(string(g-i),2)
z = strtrim(string(data[dim(0)-2,*]))
i = strtrim(string(i),2)
comments=ids+','+g_i+','+z+','+i
write_ds9_regionfile,ra(good),dec(good),comment=comments(good),filename='all_'+name+'.reg',color='blue'
stop
end
