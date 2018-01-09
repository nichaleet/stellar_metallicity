pro gallazzi

mass = [8.91,9.11,9.31,9.51,9.72,9.91,10.11,10.30,10.51,10.72,10.91,11.11,11.30,11.51,11.72,11.91]
logz = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13,0.14,0.15]
minz = [-1.11,-1.07,-1.10,-1.03,-0.97,-0.90,-0.80,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04,-0.03,-0.03]
maxz = [-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28,0.29,0.30]
age  = [9.06,9.09,9.11,90.17,9.23,9032,9046,9.61,9.73,9.82,9.87,9.90,9.92,9.94,9.95,9.96]
minage = [8.80,8.81,8.85,8.89,8.94,9.00,9.09,9.23,9.34,9.48,9.60,9.67,9.72,9.75,9.76,9.77]
maxage = [9.46,9.48,9.44,9.49,9.57,9.71,9.85,9.93,9.98,10.03,10.06,10.08,10.09,10.11,10.12,10.12]

str = {logmass:mass,z:logz,zmin:minz,zmax:maxz,age:age,agemin:minage,agemax:maxage}
mwrfits,str,'gallazzi_data.fits',/create

end