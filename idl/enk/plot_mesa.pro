pro plot_mesa, filename, xname, yname, plotname, logx=logx, logy=logy, unlogx=unlogx, unlogy=unlogy, reversex=reversex, reversey=reversey, xtitle=xtitle, ytitle=ytitle, noopen=noopen, noclose=noclose
    varnames = ' '
    openr, lun, filename, /get_lun
    skip_lun, lun, 5, /lines
    readf, lun, varnames, format='(A3000)'
    varnames = strsplit(varnames, /extract)
    nvars = n_elements(varnames)
    
    match, strlowcase(varnames), strtrim(strlowcase(xname), 2), xcolnum, w1
    match, strlowcase(varnames), strtrim(strlowcase(yname), 2), ycolnum, w2
    if w1 eq -1 then message, 'xname not found.'
    if w2 eq -1 then message, 'yname not found.'

    x = dblarr(100000)
    y = dblarr(100000)
    strin = ' '
    i = 0L

    while ~eof(lun) do begin
        readf, lun, strin, format='(A3000)'
        vars = strsplit(strin, /extract)
        x[i] = vars[xcolnum]
        y[i] = vars[ycolnum]
        i++
    endwhile
    x = x[0:i-1]
    y = y[0:i-1]

    close, lun
    free_lun, lun

    ;case strlowcase(xname) of
    ;    'model_number': x = model_number
    ;    'model': x = model_number
    ;    'modeln': x = model_number
    ;    'star_age': x = star_age
    ;    'age': x = star_age
    ;    'log_age': x = alog10(star_age)
    ;    'logage': x = alog10(star_age)
    ;    'star_mass': x = star_mass
    ;    'mass': x = star_mass
    ;    'm': x = star_mass
    ;    'log_dt': x = log_dt
    ;    'logdt': x = log_dt
    ;    'dt': x = 10.0^(log_dt)
    ;    'num_zones': x = num_zones
    ;    'nzones': x = num_zones
    ;    'zones': x = num_zones
    ;    'conv_mx1_top': x = conv_mx1_top
    ;    'convmx1top': x = conv_mx1_top
    ;    'conv_mx1_bot': x = conv_mx1_bot
    ;    'convmx1bot': x = conv_mx1_bot
    ;    'conv_mx2_top': x = conv_mx2_top
    ;    'convmx2top': x = conv_mx2_top
    ;    'conv_mx2_bot': x = conv_mx2_bot
    ;    'convmx2bot': x = conv_mx2_bot
    ;    'mx1_top': x = mx1_top
    ;    'mx1top': x = mx1_top
    ;    'mx1_bot': x = mx1_bot
    ;    'mx1bot': x = mx1_bot
    ;    'mx2_top': x = mx2_top
    ;    'mx2top': x = mx2_top
    ;    'mx2_bot': x = mx2_bot
    ;    'mx2bot': x = mx2_bot
    ;    'epsnuc_M_1': x = epsnuc_M_1
    ;    'epsnucM1': x = epsnuc_M_1
    ;    'epsnuc_M_2': x = epsnuc_M_2
    ;    'epsnucM2': x = epsnuc_M_2
    ;    'epsnuc_M_3': x = epsnuc_M_3
    ;    'epsnucM3': x = epsnuc_M_3
    ;    'epsnuc_M_4': x = epsnuc_M_4
    ;    'epsnucM4': x = epsnuc_M_4
    ;    'epsnuc_M_5': x = epsnuc_M_5
    ;    'epsnucM5': x = epsnuc_M_5
    ;    'epsnuc_M_6': x = epsnuc_M_6
    ;    'epsnucM6': x = epsnuc_M_6
    ;    'epsnuc_M_7': x = epsnuc_M_7
    ;    'epsnucM7': x = epsnuc_M_7
    ;    'epsnuc_M_8': x = epsnuc_M_8
    ;    'epsnucM8': x = epsnuc_M_8
    ;    'h1_boundary_mass': x = h1_boundary_mass
    ;    'h1boundarymass': x = h1_boundary_mass
    ;    'h1m': x = h1_boundary_mass
    ;    'he4_boundary_mass': x = he4_boundary_mass
    ;    'he4boundarymass': x = he4_boundary_mass
    ;    'he4m': x = he4_boundary_mass
    ;    'kh_timescale': x = kh_timescale
    ;    'tkh': x = kh_timescale
    ;    'log_lh': x = log_LH
    ;    'loglh': x = log_LH
    ;    'lh': x = 10.0^(log_LH)
    ;    'log_lhe': x = log_LHe
    ;    'loglhe': x = log_LHe
    ;    'lhe': x = 10.0^(log_LHe)
    ;    'log_l': x = log_L
    ;    'logl': x = log_L
    ;    'l': x = 10.0^(log_L)
    ;    'log_teff': x = log_Teff
    ;    'log_t_eff': x = log_Teff
    ;    'logteff': x = log_Teff
    ;    't_eff': x = 10.0^(log_Teff)
    ;    'teff': x = 10.0^(log_Teff)
    ;    'log_r': x = log_R
    ;    'logr': x = log_R
    ;    'r': x = 10.0^(log_R)
    ;    'log_g': x = log_g
    ;    'logg': x = log_g
    ;    'g': x = 10.0^(log_g)
    ;    'log_center_t': x = log_center_T
    ;    'logcentert': x = log_center_T
    ;    'logt_c': x = log_center_T
    ;    'logtc': x = log_center_T
    ;    'centert': x = 10.0^(log_center_T)
    ;    't_c': x = 10.0^(log_center_T)
    ;    'tc': x = 10.0^(log_center_T)
    ;    'log_center_rho': x = log_center_Rho
    ;    'logcenterrho': x = log_center_Rho
    ;    'logrho_c': x = log_center_Rho
    ;    'logrhoc': x = log_center_Rho
    ;    'centerrho': x = 10.0^(log_center_Rho)
    ;    'rho_c': x = 10.0^(log_center_Rho)
    ;    'rhoc': x = 10.0^(log_center_Rho)
    ;    'log_center_p': x = log_center_P
    ;    'logcenterp': x = log_center_P
    ;    'logp_c': x = log_center_P
    ;    'logpc': x = log_center_P
    ;    'centerp': x = 10.0^(log_center_P)
    ;    'p_c': x = 10.0^(log_center_P)
    ;    'pc': x = 10.0^(log_center_P)
    ;    'log_average_h1': x = log_average_h1
    ;    'logaverageh1': x = log_average_h1
    ;    'log<h1>': x = log_average_h1
    ;    'log <h1>': x = log_average_h1
    ;    'average_h1': x = 10.0^(log_average_h1)
    ;    'averageh1': x = 10.0^(log_average_h1)
    ;    '<h1>': x = 10.0^(log_average_h1)
    ;    'log_average_h2': x = log_average_h2
    ;    'logaverageh2': x = log_average_h2
    ;    'log<h2>': x = log_average_h2
    ;    'log <h2>': x = log_average_h2
    ;    'average_h2': x = 10.0^(log_average_h2)
    ;    'averageh2': x = 10.0^(log_average_h2)
    ;    '<h2>': x = 10.0^(log_average_h2)
    ;    'log_average_he3': x = log_average_he3
    ;    'logaveragehe3': x = log_average_he3
    ;    'log<he3>': x = log_average_he3
    ;    'log <he3>': x = log_average_he3
    ;    'average_he3': x = 10.0^(log_average_he3)
    ;    'averagehe3': x = 10.0^(log_average_he3)
    ;    '<he3>': x = 10.0^(log_average_he3)
    ;    'log_average_he4': x = log_average_he4
    ;    'logaveragehe4': x = log_average_he4
    ;    'log<he4>': x = log_average_he4
    ;    'log <he4>': x = log_average_he4
    ;    'average_he4': x = 10.0^(log_average_he4)
    ;    'averagehe4': x = 10.0^(log_average_he4)
    ;    '<he4>': x = 10.0^(log_average_he4)
    ;    'log_average_li7': x = log_average_li7
    ;    'logaverageli7': x = log_average_li7
    ;    'log<li7>': x = log_average_li7
    ;    'log <li7>': x = log_average_li7
    ;    'average_li7': x = 10.0^(log_average_li7)
    ;    'averageli7': x = 10.0^(log_average_li7)
    ;    '<li7>': x = 10.0^(log_average_li7)
    ;    'num_retries': x = num_retries
    ;    'nretries': x = num_retries
    ;    'retries': x = num_retries
    ;    'num_backups': x = num_backups
    ;    'nbackups': x = num_backups
    ;    'backups': x = num_backups
    ;    else: message, 'xname unknown.'
    ;endcase
    ;case strlowcase(yname) of
    ;    'model_number': y = model_number
    ;    'model': y = model_number
    ;    'modeln': y = model_number
    ;    'star_age': y = star_age
    ;    'age': y = star_age
    ;    'star_mass': y = star_mass
    ;    'mass': y = star_mass
    ;    'm': y = star_mass
    ;    'log_dt': y = log_dt
    ;    'logdt': y = log_dt
    ;    'dt': y = 10.0^(log_dt)
    ;    'num_zones': y = num_zones
    ;    'nzones': y = num_zones
    ;    'zones': y = num_zones
    ;    'conv_mx1_top': y = conv_mx1_top
    ;    'convmx1top': y = conv_mx1_top
    ;    'conv_mx1_bot': y = conv_mx1_bot
    ;    'convmx1bot': y = conv_mx1_bot
    ;    'conv_mx2_top': y = conv_mx2_top
    ;    'convmx2top': y = conv_mx2_top
    ;    'conv_mx2_bot': y = conv_mx2_bot
    ;    'convmx2bot': y = conv_mx2_bot
    ;    'mx1_top': y = mx1_top
    ;    'mx1top': y = mx1_top
    ;    'mx1_bot': y = mx1_bot
    ;    'mx1bot': y = mx1_bot
    ;    'mx2_top': y = mx2_top
    ;    'mx2top': y = mx2_top
    ;    'mx2_bot': y = mx2_bot
    ;    'mx2bot': y = mx2_bot
    ;    'epsnuc_M_1': y = epsnuc_M_1
    ;    'epsnucM1': y = epsnuc_M_1
    ;    'epsnuc_M_2': y = epsnuc_M_2
    ;    'epsnucM2': y = epsnuc_M_2
    ;    'epsnuc_M_3': y = epsnuc_M_3
    ;    'epsnucM3': y = epsnuc_M_3
    ;    'epsnuc_M_4': y = epsnuc_M_4
    ;    'epsnucM4': y = epsnuc_M_4
    ;    'epsnuc_M_5': y = epsnuc_M_5
    ;    'epsnucM5': y = epsnuc_M_5
    ;    'epsnuc_M_6': y = epsnuc_M_6
    ;    'epsnucM6': y = epsnuc_M_6
    ;    'epsnuc_M_7': y = epsnuc_M_7
    ;    'epsnucM7': y = epsnuc_M_7
    ;    'epsnuc_M_8': y = epsnuc_M_8
    ;    'epsnucM8': y = epsnuc_M_8
    ;    'h1_boundary_mass': y = h1_boundary_mass
    ;    'h1boundarymass': y = h1_boundary_mass
    ;    'h1m': y = h1_boundary_mass
    ;    'he4_boundary_mass': y = he4_boundary_mass
    ;    'he4boundarymass': y = he4_boundary_mass
    ;    'he4m': y = he4_boundary_mass
    ;    'kh_timescale': y = kh_timescale
    ;    'tkh': y = kh_timescale
    ;    'log_lh': y = log_LH
    ;    'loglh': y = log_LH
    ;    'lh': y = 10.0^(log_LH)
    ;    'log_lhe': y = log_LHe
    ;    'loglhe': y = log_LHe
    ;    'lhe': y = 10.0^(log_LHe)
    ;    'log_l': y = log_L
    ;    'logl': y = log_L
    ;    'l': y = 10.0^(log_L)
    ;    'log_teff': y = log_Teff
    ;    'log_t_eff': y = log_Teff
    ;    'logteff': y = log_Teff
    ;    't_eff': y = 10.0^(log_Teff)
    ;    'teff': y = 10.0^(log_Teff)
    ;    'log_r': y = log_R
    ;    'logr': y = log_R
    ;    'r': y = 10.0^(log_R)
    ;    'log_g': y = log_g
    ;    'logg': y = log_g
    ;    'g': y = 10.0^(log_g)
    ;    'log_center_t': y = log_center_T
    ;    'logcentert': y = log_center_T
    ;    'logt_c': y = log_center_T
    ;    'logtc': y = log_center_T
    ;    'centert': y = 10.0^(log_center_T)
    ;    't_c': y = 10.0^(log_center_T)
    ;    'tc': y = 10.0^(log_center_T)
    ;    'log_center_rho': y = log_center_Rho
    ;    'logcenterrho': y = log_center_Rho
    ;    'logrho_c': y = log_center_Rho
    ;    'logrhoc': y = log_center_Rho
    ;    'centerrho': y = 10.0^(log_center_Rho)
    ;    'rho_c': y = 10.0^(log_center_Rho)
    ;    'rhoc': y = 10.0^(log_center_Rho)
    ;    'log_center_p': y = log_center_P
    ;    'logcenterp': y = log_center_P
    ;    'logp_c': y = log_center_P
    ;    'logpc': y = log_center_P
    ;    'centerp': y = 10.0^(log_center_P)
    ;    'p_c': y = 10.0^(log_center_P)
    ;    'pc': y = 10.0^(log_center_P)
    ;    'log_average_h1': y = log_average_h1
    ;    'logaverageh1': y = log_average_h1
    ;    'log<h1>': y = log_average_h1
    ;    'log <h1>': y = log_average_h1
    ;    'average_h1': y = 10.0^(log_average_h1)
    ;    'averageh1': y = 10.0^(log_average_h1)
    ;    '<h1>': y = 10.0^(log_average_h1)
    ;    'log_average_h2': y = log_average_h2
    ;    'logaverageh2': y = log_average_h2
    ;    'log<h2>': y = log_average_h2
    ;    'log <h2>': y = log_average_h2
    ;    'average_h2': y = 10.0^(log_average_h2)
    ;    'averageh2': y = 10.0^(log_average_h2)
    ;    '<h2>': y = 10.0^(log_average_h2)
    ;    'log_average_he3': y = log_average_he3
    ;    'logaveragehe3': y = log_average_he3
    ;    'log<he3>': y = log_average_he3
    ;    'log <he3>': y = log_average_he3
    ;    'average_he3': y = 10.0^(log_average_he3)
    ;    'averagehe3': y = 10.0^(log_average_he3)
    ;    '<he3>': y = 10.0^(log_average_he3)
    ;    'log_average_he4': y = log_average_he4
    ;    'logaveragehe4': y = log_average_he4
    ;    'log<he4>': y = log_average_he4
    ;    'log <he4>': y = log_average_he4
    ;    'average_he4': y = 10.0^(log_average_he4)
    ;    'averagehe4': y = 10.0^(log_average_he4)
    ;    '<he4>': y = 10.0^(log_average_he4)
    ;    'log_average_li7': y = log_average_li7
    ;    'logaverageli7': y = log_average_li7
    ;    'log<li7>': y = log_average_li7
    ;    'log <li7>': y = log_average_li7
    ;    'average_li7': y = 10.0^(log_average_li7)
    ;    'averageli7': y = 10.0^(log_average_li7)
    ;    '<li7>': y = 10.0^(log_average_li7)
    ;    'num_retries': y = num_retries
    ;    'nretries': y = num_retries
    ;    'retries': y = num_retries
    ;    'num_backups': y = num_backups
    ;    'nbackups': y = num_backups
    ;    'backups': y = num_backups
    ;    else: message, 'yname unknown.'
    ;endcase

    if ~keyword_set(plotname) then plotname = 'mesa.eps'

    if keyword_set(unlogx) then x = 10.0^x
    if keyword_set(unlogy) then y = 10.0^y
    if keyword_set(logx) then x = alog10(x)
    if keyword_set(logy) then y = alog10(y/0.000954265748)
    xrange = keyword_set(reversex) ? [max(x), min(x)] : [min(x), min(x)]
    yrange = keyword_set(reversey) ? [max(y), min(y)] : [min(y), min(y)]
    if ~keyword_set(xtitle) then xtitle = (keyword_set(logx) ? 'log ' : '')+(keyword_set(unlogx) ? '10^' : '')+xname
    if ~keyword_set(ytitle) then ytitle = (keyword_set(logy) ? 'log ' : '')+(keyword_set(unlogy) ? '10^' : '')+yname
    xtitle = '!7'+xtitle
    ytitle = '!7'+ytitle

    if ~keyword_set(noopen) then begin
        setplot
        device, filename=plotname, xsize=7, ysize=7, /inches, /color, /encapsulated
        plot, x, y, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, /nodata
    endif
    oplot, x, y, color=fsc_color(keyword_set(noopen) ? 'red' : 'black')
    if ~keyword_set(noclose) then begin
        device, /close
        resetplot
    endif
    stop
end
