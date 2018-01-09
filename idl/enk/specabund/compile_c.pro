pro compile_c
    file_delete, getenv('ENK_IDL')+'specabund/smooth_gauss'+strtrim(!VERSION.MEMORY_BITS, 2)+'.so', /allow_nonexistent
    make_dll, 'smooth_gauss', 'smooth_gauss'+strtrim(!VERSION.MEMORY_BITS, 2), ['smooth_gauss', 'smooth_gauss_noivar'], compile_directory=getenv('ENK_IDL')+'specabund'

    file_delete, getenv('ENK_IDL')+'specabund/smooth_lorentz'+strtrim(!VERSION.MEMORY_BITS, 2)+'.so', /allow_nonexistent
    make_dll, 'smooth_lorentz', 'smooth_lorentz'+strtrim(!VERSION.MEMORY_BITS, 2), ['smooth_lorentz', 'smooth_lorentz_noivar'], compile_directory=getenv('ENK_IDL')+'specabund'
end
