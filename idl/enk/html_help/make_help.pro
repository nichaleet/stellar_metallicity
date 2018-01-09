pro make_help
    pros1 = file_search('/usr/local/rsi/idl/lib/', '*.pro')
    pros2 = file_search('/home/ekirby/idl/', '*.pro')
    pros = [pros1, pros2]
    mk_html_help, pros, 'idl_help.html', title="IDL Help"
    pros = file_search('/usr/local/rsi/idl/lib/', '*.pro')
    mk_html_help, pros, 'idl_help_rsi.html', title="IDL Help (RSI)"
    pros = file_search('/home/ekirby/idl/enk/', '*.pro')
    mk_html_help, pros, 'idl_help_enk.html', title="IDL Help (ENK)"
    pros = file_search('/home/ekirby/idl/npk/', '*.pro')
    mk_html_help, pros, 'idl_help_npk.html', title="IDL Help (NPK)"
    pros = file_search('/home/ekirby/idl/idlutils/', '*.pro')
    mk_html_help, pros, 'idl_help_goddard.html', title="IDL Help (Goddard)"
    pros1 = file_search('/home/ekirby/idl/deep/', '*.pro')
    pros2 = file_search('/home/ekirby/idl/deepsci/', '*.pro')
    pros = [pros1, pros2]
    mk_html_help, pros, 'idl_help_deep.html', title="IDL Help (DEEP)"    
end
