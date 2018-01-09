pro mkreg

data=read_ascii('2mass_sources.txt')
data = data.field01
ra   =data(6,*)
dec  =data(7,*)
jmag =data(9,*)
good = where(jmag gt 10 and jmag lt 15)
comment = strtrim(string(jmag),2)
write_ds9_regionfile,ra(good),dec(good),comment=comment(good),filename='2mass.reg'

end
