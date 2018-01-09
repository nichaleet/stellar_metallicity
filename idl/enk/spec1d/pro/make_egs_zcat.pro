

pro make_egs_zcat



  zcat = mrdfits(concat_dir(getenv('DEEP2PRODUCTS'), 'zcat.latest.fits'), 1, hdr, /silent)

  date = strcompress(zcat.date, /rem)

  zstruct=structindex()

  whegs = where(strmid(strcompress(zcat.objname, /rem), 0, 1) eq '1' $
                OR strmid(strcompress(zcat.objname, /rem), 0, 1) eq '5' $
                OR strmid(strcompress(zcat.objname, /rem), 0, 1) eq '6' $
                OR strmid(strcompress(zcat.objname, /rem), 0, 1) eq '7', nobj)

  zcat = zcat[whegs]

 whegs = where(strmid(strcompress(zstruct.objname, /rem), 0, 1) eq '1' $
                OR strmid(strcompress(zstruct.objname, /rem), 0, 1) eq '5' $
                OR strmid(strcompress(zstruct.objname, /rem), 0, 1) eq '6' $
                OR strmid(strcompress(zstruct.objname, /rem), 0, 1) eq '7', nobj)


  zstruct = zstruct[whegs]

  writeegsselection,zcat

  zcat=translate_objno(zcat)
  zstruct=translate_objno(zstruct)

  tmp = {OBJNO:0L, RA:0D, DEC:0D, MAGB:0.0, MAGR:0.0, MAGI:0.0, PGAL:0.0, $
         SFD_EBV:0.0, CLASS:'', SUBCLASS:'', OBJNAME:'', MASKNAME:'', $
         SLITNAME:'', DATE:'', MJD:0.0, z:0.0, zhelio:0.0, z_err:0.0, $
         rchi2:0.0, dof:0.0, VDISP:0.0, VDISP_ERR:0.0, ZQUALITY:0, $
         COMMENT:' '}

  zz = replicate(tmp, nobj)
  struct_assign, zcat, zz

  mwrfits, zz, 'zcat.egs.fits', /silent, /create
  spawn, 'gzip -1vf zcat.egs.fits'

  uzz = uniqzcat(zz)
  mwrfits, uzz, 'zcat.egs.uniq.fits', /silent, /create
  spawn, 'gzip -1vf zcat.egs.uniq.fits'
  
  mwrfits,zstruct,'zstructs.fits',/silent

; make ascii version.
  file = 'zcat.egs.dat'
  openw, unit, file, /get_lun
  printf, unit, '# DEEP2 EGS Redshift Catalog'
  printf, unit, '# OBJNO '
  printf, unit, '# RA'
  printf, unit, '# DEC'
  printf, unit, '# MAGB'
  printf, unit, '# MAGR'
  printf, unit, '# MAGI'
  printf, unit, '# PGAL'
  printf, unit, '# SFD_EBV'
  printf, unit, '# CLASS'
  printf, unit, '# SUBCLASS'
  printf, unit, '# OBJNAME'
  printf, unit, '# MASKNAME'
  printf, unit, '# SLITNAME'
  printf, unit, '# DATE'
  printf, unit, '# MJD'  
  printf, unit, '# Z'
  printf, unit, '# ZHELIO'
  printf, unit, '# Z_ERR'
  printf, unit, '# RCHI2'
  printf, unit, '# DOF'
  printf, unit, '# VDISP'
  printf, unit, '# VDISP_ERR'
  printf, unit, '# ZQUALITY'
  printf, unit, '# COMMENT'

  
  sp = ' '
  frmt = '(I8,A1,D,A1,D,A1,' + $
    'F11,A1,F11,A1,F11,A1,' + $
    'F9,A1,F9,A1,A10,A1,' + $
    'A10,A1,A8,A1,A4,A1,' + $
    'A3,A1,A12,A1,F13,A1,F11,A1,' + $
    'F11,A1,F11,A1,F11,A1,I8,A1,' + $
    'F11,A1,F11,A1,I2,A1,' + $
    'A20)'
  for ii=0,nobj-1 do printf, unit, zz[ii].objno, sp, zz[ii].ra, sp, $
    zz[ii].dec, sp, zz[ii].magb, sp, zz[ii].magr, sp, zz[ii].magi, sp, $
    zz[ii].pgal, sp, zz[ii].sfd_ebv, sp, zz[ii].class, sp, $
    zz[ii].subclass, sp, zz[ii].objname, sp, zz[ii].maskname, sp, $
    zz[ii].slitname, sp, zz[ii].date, sp, zz[ii].mjd, sp, zz[ii].z, sp, $
    zz[ii].zhelio, sp, zz[ii].z_err, sp, zz[ii].rchi2, sp, zz[ii].dof, sp, $
    zz[ii].vdisp, sp, zz[ii].vdisp_err, sp, zz[ii].zquality, sp, $
    zz[ii].comment;, format=frmt

  close, unit
  free_lun, unit

; make ascii version.
  file = 'zcat.egs.uniq.dat'
  openw, unit, file, /get_lun
  printf, unit, '# DEEP2 EGS Redshift Catalog'
  printf, unit, '# OBJNO '
  printf, unit, '# RA'
  printf, unit, '# DEC'
  printf, unit, '# MAGB'
  printf, unit, '# MAGR'
  printf, unit, '# MAGI'
  printf, unit, '# PGAL'
  printf, unit, '# SFD_EBV'
  printf, unit, '# CLASS'
  printf, unit, '# SUBCLASS'
  printf, unit, '# OBJNAME'
  printf, unit, '# MASKNAME'
  printf, unit, '# SLITNAME'
  printf, unit, '# DATE'
  printf, unit, '# MJD'  
  printf, unit, '# Z'
  printf, unit, '# ZHELIO'
  printf, unit, '# Z_ERR'
  printf, unit, '# RCHI2'
  printf, unit, '# DOF'
  printf, unit, '# VDISP'
  printf, unit, '# VDISP_ERR'
  printf, unit, '# ZQUALITY'
  printf, unit, '# COMMENT'

  
  sp = ' '
  frmt = '(I8,A1,D,A1,D,A1,' + $
    'F11,A1,F11,A1,F11,A1,' + $
    'F9,A1,F9,A1,A10,A1,' + $
    'A10,A1,A8,A1,A4,A1,' + $
    'A3,A1,A12,A1,F13,A1,F11,A1,' + $
    'F11,A1,F11,A1,F11,A1,I8,A1,' + $
    'F11,A1,F11,A1,I2,A1,' + $
    'A20)'

  zz = uzz
  nobj = n_elements(uzz)
  for ii=0,nobj-1 do printf, unit, zz[ii].objno, sp, zz[ii].ra, sp, $
    zz[ii].dec, sp, zz[ii].magb, sp, zz[ii].magr, sp, zz[ii].magi, sp, $
    zz[ii].pgal, sp, zz[ii].sfd_ebv, sp, zz[ii].class, sp, $
    zz[ii].subclass, sp, zz[ii].objname, sp, zz[ii].maskname, sp, $
    zz[ii].slitname, sp, zz[ii].date, sp, zz[ii].mjd, sp, zz[ii].z, sp, $
    zz[ii].zhelio, sp, zz[ii].z_err, sp, zz[ii].rchi2, sp, zz[ii].dof, sp, $
    zz[ii].vdisp, sp, zz[ii].vdisp_err, sp, zz[ii].zquality, sp, $
    zz[ii].comment;, format=frmt

  close, unit
  free_lun, unit






end
