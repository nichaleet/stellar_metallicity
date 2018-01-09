pro makedeimosreg_sps, mask, color=color, zrange=zrange
  if ~keyword_set(color) then color='green'
  if ~keyword_set(zrange) then zrange=[0,2]

  file = file_search('/scr2/nichal/workspace2/sps_fit/data/'+mask+ '/sps_fit.fits.gz', count=c)
  if c eq 0 then stop,'no file found'

  str = mrdfits(file[0],1)
  
  good = where(str.zspec gt zrange[0] and str.zspec lt zrange[1],cgood)
  ras = str(good).ra
  decs= str(good).dec
  imag = str(good).i
  objnames = mask+',i='+strtrim(string(imag,format='(f4.1)'))+',sn='+strtrim(string(str(good).sn,format='(f4.1)'))
  objcolor = strarr(cgood)

  objcolor(where(imag lt 19.5)) = 'purple'
  objcolor(where(imag ge 19.5 and imag lt 20.0)) = 'cyan'
  objcolor(where(imag ge 20.0 and imag lt 20.5)) = 'blue'
  objcolor(where(imag ge 20.5 and imag lt 21.0)) = 'green'
  objcolor(where(imag ge 21.0 and imag lt 21.5)) = 'yellow'
  objcolor(where(imag ge 21.5 and imag lt 22.0)) = 'orange'
  objcolor(where(imag ge 22.0 and imag lt 22.5)) = 'red'
  objcolor(where(imag ge 22.5)) = 'black'
  write_ds9_regionfile,ras,decs,comment=objnames,filename=mask+'.reg',symbol='diamond',color=color,objcolor=objcolor
end
