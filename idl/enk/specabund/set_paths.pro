pro set_paths, mask
    common paths, feh_p05, directory, vrlim, photfile, atmfile, epsfile
    case strlowcase(mask) of
        'n288': begin
            directory = getenv('M31')+'gc/n288/'
            vrlim = [-40, 0]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n288_abund.eps'
            feh_p05 = -1.39
        end
        '288l1': begin
            directory = getenv('CALTECH')+'moogify/288l1/'
            vrlim = [-80, -25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '288l1_abund.eps'
            feh_p05 = -1.39
        end
        '288l2': begin
            directory = getenv('CALTECH')+'moogify/288l2/'
            vrlim = [-80, -25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '288l2_abund.eps'
            feh_p05 = -1.39
        end
        '288l3': begin
            directory = getenv('CALTECH')+'moogify/288l3/'
            vrlim = [-80, -25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '288l3_abund.eps'
            feh_p05 = -1.39
        end
        '288l4': begin
            directory = getenv('CALTECH')+'moogify/288l4/'
            vrlim = [-80, -25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '288l4_abund.eps'
            feh_p05 = -1.39
        end
        '288l5': begin
            directory = getenv('CALTECH')+'moogify/288l5/'
            vrlim = [-80, -25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '288l5_abund.eps'
            feh_p05 = -1.39
        end
        'n1904': begin
            directory = getenv('M31')+'data/gc/n1904/'
            vrlim = [210, 250]
            photfile = getenv('M31')+'gc/n1904/n1904_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n1904/n1904_atm.txt'
            epsfile = 'n1904_abund.eps'
            feh_p05 = -1.42
        end
        '1904l1': begin
            directory = getenv('CALTECH')+'moogify/1904l1/'
            vrlim = [165, 210]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '1904l1_abund.eps'
            feh_p05 = -1.42
        end
        '1904l2': begin
            directory = getenv('CALTECH')+'moogify/1904l2/'
            vrlim = [165, 210]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '1904l2_abund.eps'
            feh_p05 = -1.42
        end
        '1904l3': begin
            directory = getenv('CALTECH')+'moogify/1904l3/'
            vrlim = [165, 210]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '1904l3_abund.eps'
            feh_p05 = -1.42
        end
        '1904l4': begin
            directory = getenv('CALTECH')+'moogify/1904l4/'
            vrlim = [165, 210]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '1904l4_abund.eps'
            feh_p05 = -1.42
        end
        'n2419': begin
            directory = getenv('M31')+'data/gc/n2419/'
            vrlim = [-20, 10]
            photfile = getenv('M31')+'gc/n2419/n2419_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n2419/n2419_atm.txt'
            epsfile = 'n2419_abund.eps'
            feh_p05 = -2.32
        end
        'n2419b': begin
            directory = getenv('CALTECH')+'moogify/n2419b/'
            vrlim = [-10, 25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n2419b_abund.eps'
            feh_p05 = -2.32
        end
        'n2419c': begin
            directory = getenv('CALTECH')+'moogify/n2419c/'
            vrlim = [-75, -40]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n2419c_abund.eps'
            feh_p05 = -2.32
        end
        'n2419c_test5': begin
            directory = getenv('CALTECH')+'moogify/n2419c_test5/'
            vrlim = [-80, -40]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n2419c_abund.eps'
            feh_p05 = -2.32
        end
        'n2419c_test6': begin
            directory = getenv('CALTECH')+'moogify/n2419c_test6/'
            vrlim = [-80, -40]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n2419c_abund.eps'
            feh_p05 = -2.32
        end
        'n4590a': begin
            directory = getenv('CALTECH')+'moogify/n4590a/'
            vrlim = [-90, -50]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n4590_abund.eps'
            feh_p05 = -2.23
        end
        'n4590b': begin
            directory = getenv('CALTECH')+'moogify/n4590b/'
            vrlim = [-90, -50]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n4590_abund.eps'
            feh_p05 = -2.23
        end
        '4590l1': begin
            directory = getenv('CALTECH')+'moogify/4590l1/'
            vrlim = [-90, -65]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n4590_abund.eps'
            feh_p05 = -2.23
        end
        'n5272c': begin
            directory = getenv('CALTECH')+'moogify/n5272c/'
            vrlim = [-150, -115]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n5272_abund.eps'
            feh_p05 = -1.50
        end
        'm3_iv-101': begin
            directory = getenv('CALTECH')+'moogify/m3_iv-101/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'm3_iv-101_abund.eps'
            feh_p05 = -1.50
        end
        'n5634a': begin
            directory = getenv('CALTECH')+'moogify/n5634a/'
            vrlim = [-60, -15]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n5634_abund.eps'
            feh_p05 = -1.88
        end
        'n5634b': begin
            directory = getenv('CALTECH')+'moogify/n5634b/'
            vrlim = [-20, 10]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n5634_abund.eps'
            feh_p05 = -1.88
        end
        '5897a': begin
            directory = getenv('CALTECH')+'moogify/5897a/'
            vrlim = [100, 160]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '5897a_abund.eps'
            feh_p05 = -1.90
        end
        '5897l1': begin
            directory = getenv('CALTECH')+'moogify/5897l1/'
            vrlim = [110, 155]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '5897l1_abund.eps'
            feh_p05 = -1.90
        end
        '5897l2': begin
            directory = getenv('CALTECH')+'moogify/5897l2/'
            vrlim = [95, 135]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '5897l2_abund.eps'
            feh_p05 = -1.90
        end
        '5897l3': begin
            directory = getenv('CALTECH')+'moogify/5897l3/'
            vrlim = [100, 140]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '5897l3_abund.eps'
            feh_p05 = -1.90
        end
        '5897l4': begin
            directory = getenv('CALTECH')+'moogify/5897l4/'
            vrlim = [110, 130]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '5897l4_abund.eps'
            feh_p05 = -1.90
        end
        'n6205': begin
            directory = getenv('M31')+'data/gc/n6205/'
            vrlim = [-250, -210]
            photfile = getenv('M31')+'gc/n6205/n6205_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n6205/n6205_atm.txt'
            epsfile = 'n6205_abund.eps'
            feh_p05 = -1.57
        end
        '6229a': begin
            directory = getenv('CALTECH')+'moogify/6229a/'
            vrlim = [-150, -105]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6229a_abund.eps'
            feh_p05 = -1.47
        end
        'n6341a': begin
            directory = getenv('CALTECH')+'moogify/n6341a/'
            vrlim = [-140, -100]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n6341_abund.eps'
            feh_p05 = -2.34
        end
        'n6341b': begin
            directory = getenv('CALTECH')+'moogify/n6341b/'
            vrlim = [-140, -100]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n6341_abund.eps'
            feh_p05 = -2.34
        end 
        '6341l1': begin
            directory = getenv('CALTECH')+'moogify/6341l1/'
            vrlim = [-140, -90]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n6341_abund.eps'
            feh_p05 = -2.34
        end
        '6341l2': begin
            directory = getenv('CALTECH')+'moogify/6341l2/'
            vrlim = [-140, -90]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n6341_abund.eps'
            feh_p05 = -2.34
        end
       'n6656b': begin
            directory = getenv('CALTECH')+'moogify/n6656b/'
            vrlim = [-140, -105]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n6656b_abund.eps'
            feh_p05 = -1.56
        end
        '6779l1': begin
            directory = getenv('CALTECH')+'moogify/6779l1/'
            vrlim = [-155, -115]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6779l1_abund.eps'
            feh_p05 = -999.0
        end
        '6779l2': begin
            directory = getenv('CALTECH')+'moogify/6779l2/'
            vrlim = [-135, -100]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6779l2_abund.eps'
            feh_p05 = -999.0
        end
        '6779l3': begin
            directory = getenv('CALTECH')+'moogify/6779l3/'
            vrlim = [-140, -110]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6779l3_abund.eps'
            feh_p05 = -999.0
        end
       'n6838': begin
            directory = getenv('M31')+'data/gc/n6838/'
            vrlim = [-10, 10]
            photfile = getenv('M31')+'gc/n6838/n6838_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n6838/n6838_atm.txt'
            epsfile = 'n6838_abund.eps'
            feh_p05 = -0.76
        end
       '6864ab': begin
            directory = getenv('CALTECH')+'moogify/6864aB/'
            vrlim = [-200, -160]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6864aB_abund.eps'
            feh_p05 = -1.29
        end
       '6864l1': begin
            directory = getenv('CALTECH')+'moogify/6864l1/'
            vrlim = [-215, -190]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6864aB_abund.eps'
            feh_p05 = -1.29
        end
       '6864l2': begin
            directory = getenv('CALTECH')+'moogify/6864l2/'
            vrlim = [-215, -190]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '6864aB_abund.eps'
            feh_p05 = -1.29
        end
        'n7006': begin
            directory = getenv('M31')+'data/gc/n7006/'
            vrlim = [-360, -350]
            photfile = getenv('M31')+'gc/n7006/n7006_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n7006/n7006_atm.txt'
            epsfile = 'n7006_abund.eps'
            feh_p05 = -1.55
        end
        '7006a': begin
            directory = getenv('CALTECH')+'moogify/7006a/'
            vrlim = [-430, -390]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7006a_abund.eps'
            feh_p05 = -1.55
        end
        'n7078': begin
            directory = getenv('M31')+'data/gc/n7078/'
            vrlim = [-100, -60]
            photfile = getenv('M31')+'gc/n7078/n7078_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n7078/n7078_atm.txt'
            epsfile = 'n7078_abund.eps'
            feh_p05 = -2.38
        end
        'n7078d': begin
            directory = getenv('CALTECH')+'moogify/n7078d/'
            vrlim = [-100, -60]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7078d_abund.eps'
            feh_p05 = -2.38
        end
        'n7078e': begin
            directory = getenv('CALTECH')+'moogify/n7078e/'
            vrlim = [-100, -60]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7078e_abund.eps'
            feh_p05 = -2.38
        end
        '7078l1b': begin
            directory = getenv('CALTECH')+'moogify/7078l1B/'
            vrlim = [-130, -70]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7078l1B_abund.eps'
            feh_p05 = -2.38
        end
        '7078l2b': begin
            directory = getenv('CALTECH')+'moogify/7078l2B/'
            vrlim = [-130, -70]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7078l2B_abund.eps'
            feh_p05 = -2.38
        end
        '7078l3b': begin
            directory = getenv('CALTECH')+'moogify/7078l3B/'
            vrlim = [-130, -70]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7078l3B_abund.eps'
            feh_p05 = -2.38
        end
        '7078l4b': begin
            directory = getenv('CALTECH')+'moogify/7078l4B/'
            vrlim = [-130, -70]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7078l4B_abund.eps'
            feh_p05 = -2.38
        end
        '7078l5b': begin
            directory = getenv('CALTECH')+'moogify/7078l5B/'
            vrlim = [-130, -70]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7078l5B_abund.eps'
            feh_p05 = -2.38
        end
        'n7089b': begin
            directory = getenv('CALTECH')+'moogify/n7089b/'
            vrlim = [0, 40]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089b_abund.eps'
            feh_p05 = -1.65
        end
        '7089c': begin
            directory = getenv('CALTECH')+'moogify/7089c/'
            vrlim = [-55, -5]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7089c_abund.eps'
            feh_p05 = -1.65
        end
        '7089l1': begin
            directory = getenv('CALTECH')+'moogify/7089l1/'
            vrlim = [-55, -15]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089bl1_abund.eps'
            feh_p05 = -1.65
        end
        '7089l2': begin
            directory = getenv('CALTECH')+'moogify/7089l2/'
            vrlim = [-50, -10]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089bl2_abund.eps'
            feh_p05 = -1.65
        end
        '7089l3': begin
            directory = getenv('CALTECH')+'moogify/7089l3/'
            vrlim = [-60, 0]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089l3_abund.eps'
            feh_p05 = -1.65
        end
        '7089l4': begin
            directory = getenv('CALTECH')+'moogify/7089l4/'
            vrlim = [-25, 25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089l4_abund.eps'
            feh_p05 = -1.65
        end
        '7089l5': begin
            directory = getenv('CALTECH')+'moogify/7089l5/'
            vrlim = [-25, 25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089l5_abund.eps'
            feh_p05 = -1.65
        end
        '7089l6': begin
            directory = getenv('CALTECH')+'moogify/7089l6/'
            vrlim = [-25, 25]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089l6_abund.eps'
            feh_p05 = -1.65
        end
        '7089m1': begin
            directory = getenv('CALTECH')+'moogify/7089m1/'
            vrlim = [-50, -10]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089m1_abund.eps'
            feh_p05 = -1.65
        end
        '7089m2': begin
            directory = getenv('CALTECH')+'moogify/7089m2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7089m2_abund.eps'
            feh_p05 = -1.65
        end
        'n7099': begin
            directory = getenv('CALTECH')+'moogify/n7099/'
            vrlim = [-165, -140]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'n7099_abund.eps'
            feh_p05 = -2.27
        end
        '7099l1': begin
            directory = getenv('CALTECH')+'moogify/7099l1/'
            vrlim = [-205, -145]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l1_abund.eps'
            feh_p05 = -2.27
        end
        '7099l2': begin
            directory = getenv('CALTECH')+'moogify/7099l2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l2_abund.eps'
            feh_p05 = -2.27
        end
        '7099l3': begin
            directory = getenv('CALTECH')+'moogify/7099l3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l3_abund.eps'
            feh_p05 = -2.27
        end
        '7099l4': begin
            directory = getenv('CALTECH')+'moogify/7099l4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l4_abund.eps'
            feh_p05 = -2.27
        end
        '7099l5': begin
            directory = getenv('CALTECH')+'moogify/7099l5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l5_abund.eps'
            feh_p05 = -2.27
        end
        '7099l6': begin
            directory = getenv('CALTECH')+'moogify/7099l6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l6_abund.eps'
            feh_p05 = -2.27
        end
        '7099l7': begin
            directory = getenv('CALTECH')+'moogify/7099l7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = '7099l7_abund.eps'
            feh_p05 = -2.27
        end
        'n7492': begin
            directory = getenv('M31')+'data/gc/n7492/'
            vrlim = [-160, -140]
            photfile = getenv('M31')+'gc/n7492/n7492_hires_phot.txt'
            atmfile = getenv('M31')+'gc/n7492/n7492_atm.txt'
            epsfile = 'n7492_abund.eps'
            feh_p05 = -1.85
        end
        'pal2': begin
            directory = getenv('M31')+'gc/pal2/'
            vrlim = [-150, -110]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal2_abund.eps'
            feh_p05 = -1.30       ;Harris 1996
        end
        'pal13': begin
            directory = getenv('CALTECH')+'moogify/pal13/'
            vrlim = [30, 50]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal13-3': begin
            directory = getenv('CALTECH')+'gc/pal13_geha/Pal13-3/'
            vrlim = [-999, 999]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal13-5': begin
            directory = getenv('CALTECH')+'gc/pal13_geha/Pal13-5/'
            vrlim = [-999, 999]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal13-6': begin
            directory = getenv('CALTECH')+'gc/pal13_geha/Pal13-6/'
            vrlim = [-999, 999]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal13-7': begin
            directory = getenv('CALTECH')+'gc/pal13_geha/Pal13-7/'
            vrlim = [-999, 999]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal13-8': begin
            directory = getenv('CALTECH')+'gc/pal13_geha/Pal13-8/'
            vrlim = [-999, 999]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal13_abund.eps'
            feh_p05 = -1.74       ;Harris 1996
        end
        'pal14a': begin
            directory = getenv('CALTECH')+'moogify/pal14a/'
            vrlim = [80, 120]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'pal14a_abund.eps'
            feh_p05 = -1.62       ;Harris 1996
        end
        'eri_1': begin
            directory = getenv('CALTECH')+'moogify/Eri_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'eri_abund.eps'
            feh_p05 = -1.46       ;Harris 1996
        end
        'eri_2': begin
            directory = getenv('CALTECH')+'moogify/Eri_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            atmfile = 'none'
            epsfile = 'eri_abund.eps'
            feh_p05 = -1.46       ;Harris 1996
        end
        'ng1904r': begin
            directory = getenv('M31')+'gc2/ng1904/'
            vrlim = [210, 240]
            photfile = getenv('M31')+'/gc2/n1904/n1904_phot.fits'
            epsfile = 'ng1904r_abund.eps'
            feh_p05 = -1.42
        end
        'ng5024r': begin
            directory = getenv('M31')+'gc2/ng5024/'
            vrlim = [-95, -60]
            photfile = getenv('M31')+'/gc2/n5024/n5024_phot.fits'
            epsfile = 'ng5024r_abund.eps'
            ;feh_p05 = -1.84         ;Martell, Smith, & Briley 2008
            feh_p05 = -999
        end
        'ng5053r': begin
            directory = getenv('M31')+'gc2/ng5053/'
            vrlim = [10, 45]
            photfile = getenv('M31')+'/gc2/n5053/n5053_phot.fits'
            epsfile = 'ng5053r_abund.eps'
            ;feh_p05 = -2.29         ;Harris 1996
            feh_p05 = -999
        end
        'ng5904r': begin
            directory = getenv('M31')+'gc2/ng5904/'
            vrlim = [5, 45]
            photfile = getenv('M31')+'/gc2/n5904/n5904_phot.fits'
            epsfile = 'ng5904r_abund.eps'
            feh_p05 = -1.30
        end
        'ng5904_byteff': begin
            directory = getenv('CALTECH')+'moogify/ng5904_byteff/'
            vrlim = [-1000, 1000]
            photfile = getenv('M31')+'/gc2/n5904/n5904_byteff_phot.fits'
            epsfile = 'ng5904_byteff_abund.eps'
            feh_p05 = -1.30
        end
        'ng5904_bycolor': begin
            directory = getenv('CALTECH')+'moogify/ng5904_bycolor/'
            vrlim = [-1000, 1000]
            photfile = getenv('M31')+'/gc2/n5904/n5904_bycolor_phot.fits'
            epsfile = 'ng5904_bycolor_abund.eps'
            feh_p05 = -1.30
        end
        'ter5a': begin
            directory = getenv('CALTECH')+'moogify/ter5a/'
            vrlim = [-1000, 1000]
            photfile = getenv('CALTECH')+'/gc/ter5/ter5_phot.fits'
            epsfile = 'ter5a_abund.eps'
            feh_p05 = 0.00
        end
        'ter5b': begin
            directory = getenv('CALTECH')+'moogify/ter5b/'
            vrlim = [-1000, 1000]
            photfile = getenv('CALTECH')+'/gc/ter5/ter5_phot.fits'
            epsfile = 'ter5b_abund.eps'
            feh_p05 = 0.00
        end
        'ter5c': begin
            directory = getenv('CALTECH')+'moogify/ter5c/'
            vrlim = [-1000, 1000]
            photfile = getenv('CALTECH')+'/gc/ter5/ter5_phot.fits'
            epsfile = 'ter5b_abund.eps'
            feh_p05 = 0.00
        end
        'leoi': begin
            directory = getenv('M31')+'leoi/deimos/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/deimos/leoi_phot.txt'
            epsfile = 'leoi_abund.eps'
            feh_p05 = -999
        end
        'lin1_1': begin
            directory = getenv('M31')+'leoi/LIN1_1/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN1_1/LIN1_1_phot.fits'
            epsfile = 'LIN1_3_abund.eps'
            feh_p05 = -999
        end
        'lin1_2': begin
            directory = getenv('M31')+'leoi/LIN1_2/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN1_2/LIN1_2_phot.fits'
            epsfile = 'LIN1_4_abund.eps'
            feh_p05 = -999
        end
        'lin1_3': begin
            directory = getenv('M31')+'leoi/LIN1_3/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN1_3/LIN1_3_phot.fits'
            epsfile = 'LIN1_3_abund.eps'
            feh_p05 = -999
        end
        'lin1_4': begin
            directory = getenv('M31')+'leoi/LIN1_4/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN1_4/LIN1_4_phot.fits'
            epsfile = 'LIN1_4_abund.eps'
            feh_p05 = -999
        end
        'lin2_1': begin
            directory = getenv('M31')+'leoi/LIN2_1/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN2_1/LIN2_1_phot.fits'
            epsfile = 'LIN2_3_abund.eps'
            feh_p05 = -999
        end
        'lin2_2': begin
            directory = getenv('M31')+'leoi/LIN2_2/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN2_2/LIN2_2_phot.fits'
            epsfile = 'LIN2_4_abund.eps'
            feh_p05 = -999
        end
        'lin2_3': begin
            directory = getenv('M31')+'leoi/LIN2_3/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN2_3/LIN2_3_phot.fits'
            epsfile = 'LIN2_3_abund.eps'
            feh_p05 = -999
        end
        'lin2_4': begin
            directory = getenv('M31')+'leoi/LIN2_4/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN2_4/LIN2_4_phot.fits'
            epsfile = 'LIN2_4_abund.eps'
            feh_p05 = -999
        end
        'lin3_1': begin
            directory = getenv('M31')+'leoi/LIN3_1/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN3_1/LIN3_1_phot.fits'
            epsfile = 'LIN3_3_abund.eps'
            feh_p05 = -999
        end
        'lin3_2': begin
            directory = getenv('M31')+'leoi/LIN3_2/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN3_2/LIN3_2_phot.fits'
            epsfile = 'LIN3_4_abund.eps'
            feh_p05 = -999
        end
        'lin3_3': begin
            directory = getenv('M31')+'leoi/LIN3_3/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN3_3/LIN3_3_phot.fits'
            epsfile = 'LIN3_3_abund.eps'
            feh_p05 = -999
        end
        'lin3_4': begin
            directory = getenv('M31')+'leoi/LIN3_4/'
            vrlim = [225, 320]
            photfile = getenv('M31')+'leoi/LIN3_4/LIN3_4_phot.fits'
            epsfile = 'LIN3_4_abund.eps'
            feh_p05 = -999
        end
        'l2a': begin
            directory = getenv('M31')+'dsph/leoii/L2A/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'l2b': begin
            directory = getenv('M31')+'dsph/leoii/L2B/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'l2c': begin
            directory = getenv('M31')+'dsph/leoii/L2C/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'l2d': begin
            directory = getenv('M31')+'dsph/leoii/L2D/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'l2e': begin
            directory = getenv('M31')+'dsph/leoii/L2E/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'l2f': begin
            directory = getenv('M31')+'dsph/leoii/L2F/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoii_abund.eps'
            feh_p05 = -999
        end
        'cb': begin
            directory = getenv('M31')+'udwarf/CB/'
            ;vrlim = [-40, 120]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'CB_abund.eps'
            feh_p05 = -999
        end
        'combi1': begin
            directory = getenv('ahome')+'deimos/Combi1/'
            vrlim = [-1000, 1000]
        end
        'combi2': begin
            directory = getenv('ahome')+'deimos/Combi2/'
            vrlim = [-1000, 1000]
        end
        'cvni': begin
            directory = getenv('M31')+'udwarf/CVnI/'
            ;vrlim = [-20, 70]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'CVnI_abund.eps'
            feh_p05 = -999
        end
        'cvibi1': begin
            directory = getenv('ahome')+'deimos/CVIbi1/'
            vrlim = [-1000, 1000]
        end
        'cvibi2': begin
            directory = getenv('ahome')+'deimos/CVIbi2/'
            vrlim = [-1000, 1000]
        end
        'cvnii': begin
            directory = getenv('M31')+'udwarf/CVnII/'
            ;vrlim = [-160, -120]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'CVnII_abund.eps'
            feh_p05 = -999
        end
        'cvnii-3': begin
            directory = getenv('CALTECH')+'moogify/CVnII-3/'
            ;vrlim = [-160, -120]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'CVnII-3_abund.eps'
            feh_p05 = -999
        end
        'herc': begin
            directory = getenv('M31')+'udwarf/Herc/'
            ;vrlim = [0, 50]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Herc_abund.eps'
            feh_p05 = -999
        end
        'leoiv': begin
            directory = getenv('M31')+'udwarf/LeoIV/'
            ;vrlim = [90, 150]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'LeoIV_abund.eps'
            feh_p05 = -999
        end
        'leo4-2': begin
            directory = getenv('CALTECH')+'moogify/leo4-2/'
            ;vrlim = [90, 150]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leo4-2_abund.eps'
            feh_p05 = -999
        end
        'leot': begin
            directory = getenv('M31')+'udwarf/LeoT/'
            ;vrlim = [0, 100]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'LeoT_abund.eps'
            feh_p05 = -999
        end
        'umai': begin
            directory = getenv('M31')+'udwarf/UMaI/'
            ;vrlim = [-70, -20]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'UMaI_abund.eps'
            feh_p05 = -999
        end
        'umaii': begin
            directory = getenv('M31')+'udwarf/UMaII/'
            ;vrlim = [-140, -70]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'UMaII_abund.eps'
            feh_p05 = -999
        end
        'uma2-4': begin
            directory = getenv('ahome')+'deimos/UMa2-4/'
            vrlim = [-1000, 1000]
        end
        'uma2-5': begin
            directory = getenv('ahome')+'deimos/UMa2-5/'
            vrlim = [-1000, 1000]
        end
        'uma2-7': begin
            directory = getenv('ahome')+'deimos/UMa2-7/'
            vrlim = [-1000, 1000]
        end
        'booi-5': begin
            directory = getenv('CALTECH')+'moogify/BooI-5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'BooI-5_abund.eps'
            feh_p05 = -999
        end
        'boo2_1': begin
            directory = getenv('CALTECH')+'moogify/Boo2_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'boo2_2': begin
            directory = getenv('CALTECH')+'moogify/Boo2_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'boo2_3': begin
            directory = getenv('CALTECH')+'moogify/Boo2_3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'boo2_4': begin
            directory = getenv('CALTECH')+'moogify/Boo2_4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'boo2-5': begin
            directory = getenv('CALTECH')+'moogify/Boo2-5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'booiic': begin
            directory = getenv('CALTECH')+'moogify/booiic/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'boo2_abund.eps'
            feh_p05 = -999
        end
        'seg2_1': begin
            directory = getenv('CALTECH')+'moogify/seg2_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_1_abund.eps'
            feh_p05 = -999
        end
        'seg2_2': begin
            directory = getenv('CALTECH')+'moogify/seg2_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_2_abund.eps'
            feh_p05 = -999
        end
        'seg2_3': begin
            directory = getenv('CALTECH')+'moogify/seg2_3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_3_abund.eps'
            feh_p05 = -999
        end
        'seg2_4': begin
            directory = getenv('CALTECH')+'moogify/seg2_4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_4_abund.eps'
            feh_p05 = -999
        end
        'seg2_5': begin
            directory = getenv('CALTECH')+'moogify/seg2_5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_5_abund.eps'
            feh_p05 = -999
        end
        'seg2_6': begin
            directory = getenv('CALTECH')+'moogify/seg2_6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_6_abund.eps'
            feh_p05 = -999
        end
        'seg2_7': begin
            directory = getenv('CALTECH')+'moogify/seg2_7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_7_abund.eps'
            feh_p05 = -999
        end
        'seg2_8': begin
            directory = getenv('CALTECH')+'moogify/seg2_8/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_8_abund.eps'
            feh_p05 = -999
        end
        'seg2_9': begin
            directory = getenv('CALTECH')+'moogify/seg2_9/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_9_abund.eps'
            feh_p05 = -999
        end
        'seg2_stack': begin
            directory = getenv('CALTECH')+'moogify/seg2_stack/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'seg2_stack_abund.eps'
            feh_p05 = -999
        end
        'ari-2': begin
            directory = getenv('CALTECH')+'moogify/Ari-2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'ari-2_abund.eps'
            feh_p05 = -999
        end
        'segue1': begin
            directory = getenv('M31')+'udwarf/SEGUE1/'
            vrlim = [150, 250]
            photfile = 'none'
            epsfile = 'SEGUE1_abund.eps'
            feh_p05 = -999
        end
        'segue1_2': begin
            directory = getenv('M31')+'udwarf/segue1_2/'
            vrlim = [150, 250]
            photfile = 'none'
            epsfile = 'segue1_2_abund.eps'
            feh_p05 = -999
        end
        'segue1_2009': begin
            directory = getenv('ahome')+'deimos/segue1_2009/'
            vrlim = [0, 300]
            photfile = 'none'
            epsfile = 'segue1_2009_abund.eps'
            feh_p05 = -999
        end
        'segue1_2010': begin
            directory = getenv('ahome')+'deimos/segue1_2010/'
            vrlim = [0, 300]
            photfile = 'none'
            epsfile = 'segue1_2010_abund.eps'
            feh_p05 = -999
        end
        'segue1-1': begin
            directory = getenv('CALTECH')+'moogify/Segue1-1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Segue1-1_abund.eps'
            feh_p05 = -999
        end
        'segue1-2': begin
            directory = getenv('CALTECH')+'moogify/Segue1-2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Segue1-2_abund.eps'
            feh_p05 = -999
        end
        'segue1-3': begin
            directory = getenv('CALTECH')+'moogify/Segue1-3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Segue1-3_abund.eps'
            feh_p05 = -999
        end
        'segue1-c': begin
            directory = getenv('CALTECH')+'moogify/Segue1-C/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Segue1-C_abund.eps'
            feh_p05 = -999
        end
        'segtide1': begin
            directory = getenv('CALTECH')+'moogify/segtide1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide1_abund.eps'
            feh_p05 = -999
        end
        'segtide2': begin
            directory = getenv('CALTECH')+'moogify/segtide2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide2_abund.eps'
            feh_p05 = -999
        end
        'segtide3': begin
            directory = getenv('CALTECH')+'moogify/segtide3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide3_abund.eps'
            feh_p05 = -999
        end
        'segtide4': begin
            directory = getenv('CALTECH')+'moogify/segtide4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide4_abund.eps'
            feh_p05 = -999
        end
        'segtide5': begin
            directory = getenv('CALTECH')+'moogify/segtide5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide5_abund.eps'
            feh_p05 = -999
        end
        'segtide6': begin
            directory = getenv('CALTECH')+'moogify/segtide6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide6_abund.eps'
            feh_p05 = -999
        end
        'segtide7': begin
            directory = getenv('CALTECH')+'moogify/segtide7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'segtide7_abund.eps'
            feh_p05 = -999
        end
        'w1': begin
            directory = getenv('M31')+'udwarf/W1/'
            vrlim = [-150, 70]
            photfile = 'none'
            epsfile = 'W1_abund.eps'
            feh_p05 = -999
        end
        'w1_1': begin
            directory = getenv('chome')+'deimos/W1/W1_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'W1_abund.eps'
            feh_p05 = -999
        end
        'w1_2': begin
            directory = getenv('chome')+'deimos/W1/W1_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'W1_abund.eps'
            feh_p05 = -999
        end
        'w1_3': begin
            directory = getenv('chome')+'deimos/W1/W1_3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'W1_abund.eps'
            feh_p05 = -999
        end
        'w1_2010': begin
            directory = getenv('M31')+'udwarf/W1_2010/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'W1_2010_abund.eps'
            feh_p05 = -999
        end
        'hyaii': begin
            directory = getenv('CALTECH')+'moogify/HyaII/'
            ;vrlim = [-70, -20]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'HyaII_abund.eps'
            feh_p05 = -999
        end
        'pscii': begin
            directory = getenv('CALTECH')+'moogify/PscII/'
            ;vrlim = [-70, -20]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'PscII_abund.eps'
            feh_p05 = -999
        end
        'crti': begin
            directory = getenv('CALTECH')+'moogify/CrtI/'
            ;vrlim = [-70, -20]
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'CrtI_abund.eps'
            feh_p05 = -999
        end
        'halo': begin
            directory = getenv('M31_DATA')+'halo/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'halo_abund.eps'
            feh_p05 = -999
        end
        'halo2': begin
            directory = getenv('CALTECH')+'moogify/halo2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'halo2_abund.eps'
            feh_p05 = -999
        end
        'd1_2': begin
            directory = getenv('chome')+'deimos/d1_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'd1_2_abund.eps'
            feh_p05 = -999
        end
        'd2_1': begin
            directory = getenv('M31_DATA')+'and2/d2_1/'
            vrlim = [-250, -170]
            photfile = 'none'
            epsfile = 'd2_1_abund.eps'
            feh_p05 = -999
        end
        'd2_2': begin
            directory = getenv('M31_DATA')+'and2/d2_2/'
            vrlim = [-260, -170]
            photfile = 'none'
            epsfile = 'd2_2_abund.eps'
            feh_p05 = -999
        end
        'd2_3': begin
            directory = getenv('M31_DATA')+'and2/d2_3/'
            vrlim = [-250, -170]
            photfile = 'none'
            epsfile = 'd2_3_abund.eps'
            feh_p05 = -999
        end
        'd2_4': begin
            directory = getenv('M31_DATA')+'and2/d2_4/'
            vrlim = [-210, -160]
            photfile = 'none'
            epsfile = 'd2_4_abund.eps'
            feh_p05 = -999
        end
        'd2_5': begin
            directory = getenv('M31_DATA')+'and2/d2_5/'
            vrlim = [-230, -170]
            photfile = 'none'
            epsfile = 'd2_5_abund.eps'
            feh_p05 = -999
        end
        'd2_7': begin
            directory = getenv('M31_DATA')+'and2/d2_7/'
            vrlim = [-250, -170]
            photfile = 'none'
            epsfile = 'd2_7_abund.eps'
            feh_p05 = -999
        end
        'd2_8': begin
            directory = getenv('M31_DATA')+'and2/d2_8/'
            vrlim = [-250, -170]
            photfile = 'none'
            epsfile = 'd2_8_abund.eps'
            feh_p05 = -999
        end
        'd2_9': begin
            directory = getenv('M31_DATA')+'and2/d2_9/'
            vrlim = [-250, -150]
            photfile = 'none'
            epsfile = 'd2_9_abund.eps'
            feh_p05 = -999
        end
        'd2_10': begin
            directory = getenv('M31_DATA')+'and2/d2_10/'
            vrlim = [-250, -170]
            photfile = 'none'
            epsfile = 'd2_10_abund.eps'
            feh_p05 = -999
        end
        'for1b': begin
            directory = getenv('M31')+'dsph/for/for1B/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'for1B_abund.eps'
            feh_p05 = -999
        end
        'for3b': begin
            directory = getenv('M31')+'dsph/for/for3B/'
            vrlim = [-10, 85]
            photfile = 'none'
            epsfile = 'for3B_abund.eps'
            feh_p05 = -999
        end
        'for4b': begin
            directory = getenv('M31')+'dsph/for/for4B/'
            vrlim = [50, 120]
            photfile = 'none'
            epsfile = 'for4B_abund.eps'
            feh_p05 = -999
        end
        'for6': begin
            directory = getenv('M31')+'dsph/for/for6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'for6_abund.eps'
            feh_p05 = -999
        end
        'for7': begin
            directory = getenv('M31')+'dsph/for/for7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'for7_abund.eps'
            feh_p05 = -999
        end
        'scl1': begin
            directory = getenv('M31')+'dsph/scl/scl1/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'scl1_abund.eps'
            feh_p05 = -999
        end
        'scl2': begin
            directory = getenv('M31')+'dsph/scl/scl2/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'scl2_abund.eps'
            feh_p05 = -999
        end
        'scl3': begin
            directory = getenv('M31')+'dsph/scl/scl3/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'scl3_abund.eps'
            feh_p05 = -999
        end
        'scl5': begin
            directory = getenv('M31')+'dsph/scl/scl5/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'scl5_abund.eps'
            feh_p05 = -999
        end
        'scl6': begin
            directory = getenv('M31')+'dsph/scl/scl6/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'scl6_abund.eps'
            feh_p05 = -999
        end
        'sclhrs': begin
            directory = getenv('M31')+'dsph/scl/'
            vrlim = [50, 130]
            photfile = 'none'
            epsfile = 'sclhrs_abund.eps'
            feh_p05 = -999
        end
        'umi1': begin
            directory = getenv('M31')+'dsph/umi/umi1/'
            vrlim = [-275, -175]
            photfile = 'none'
            epsfile = 'umi1_abund.eps'
            feh_p05 = -999
        end
        'umi2': begin
            directory = getenv('M31')+'dsph/umi/umi2/'
            vrlim = [-275, -175]
            photfile = 'none'
            epsfile = 'umi2_abund.eps'
            feh_p05 = -999
        end
        'umi3': begin
            directory = getenv('M31')+'dsph/umi/umi3/'
            vrlim = [-275, -175]
            photfile = 'none'
            epsfile = 'umi3_abund.eps'
            feh_p05 = -999
        end
        'umi6': begin
            directory = getenv('M31')+'dsph/umi/umi6/'
            vrlim = [-275, -175]
            photfile = 'none'
            epsfile = 'umi6_abund.eps'
            feh_p05 = -999
        end
        'umima1': begin
            directory = getenv('CALTECH')+'moogify/umima1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umima1_abund.eps'
            feh_p05 = -999
        end
        'umima2': begin
            directory = getenv('CALTECH')+'moogify/umima2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umima2_abund.eps'
            feh_p05 = -999
        end
        'umima3': begin
            directory = getenv('CALTECH')+'moogify/umima3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umima3_abund.eps'
            feh_p05 = -999
        end
        'umimi1': begin
            directory = getenv('CALTECH')+'moogify/umimi1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umimi1_abund.eps'
            feh_p05 = -999
        end
        'umimi2': begin
            directory = getenv('CALTECH')+'moogify/umimi2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umimi2_abund.eps'
            feh_p05 = -999
        end
        'umimi3': begin
            directory = getenv('CALTECH')+'moogify/umimi3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umimi3_abund.eps'
            feh_p05 = -999
        end
        'umix1': begin
            directory = getenv('CALTECH')+'moogify/umix1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umix1_abund.eps'
            feh_p05 = -999
        end
        'umix2': begin
            directory = getenv('CALTECH')+'moogify/umix2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umix2_abund.eps'
            feh_p05 = -999
        end
        'umix4': begin
            directory = getenv('CALTECH')+'moogify/umix4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umix4_abund.eps'
            feh_p05 = -999
        end
        'uss-1': begin
            directory = getenv('CALTECH')+'moogify/uss-1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-1_abund.eps'
            feh_p05 = -999
        end
        'uss-2': begin
            directory = getenv('CALTECH')+'moogify/uss-2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-2_abund.eps'
            feh_p05 = -999
        end
        'uss-3': begin
            directory = getenv('CALTECH')+'moogify/uss-3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-3_abund.eps'
            feh_p05 = -999
        end
        'uss-4': begin
            directory = getenv('CALTECH')+'moogify/uss-4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-4_abund.eps'
            feh_p05 = -999
        end
        'uss-5': begin
            directory = getenv('CALTECH')+'moogify/uss-5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-5_abund.eps'
            feh_p05 = -999
        end
        'uss-6': begin
            directory = getenv('CALTECH')+'moogify/uss-6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-6_abund.eps'
            feh_p05 = -999
        end
        'uss-7': begin
            directory = getenv('CALTECH')+'moogify/uss-7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-7_abund.eps'
            feh_p05 = -999
        end
        'uss-8': begin
            directory = getenv('CALTECH')+'moogify/uss-8/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-8_abund.eps'
            feh_p05 = -999
        end
        'uss-9c': begin
            directory = getenv('CALTECH')+'moogify/uss-9c/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-9c_abund.eps'
            feh_p05 = -999
        end
        'uss-10': begin
            directory = getenv('CALTECH')+'moogify/uss-10/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-10_abund.eps'
            feh_p05 = -999
        end
        'uss-11': begin
            directory = getenv('CALTECH')+'moogify/uss-11/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-11_abund.eps'
            feh_p05 = -999
        end
        'uss-12': begin
            directory = getenv('CALTECH')+'moogify/uss-12/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'uss-12_abund.eps'
            feh_p05 = -999
        end
        'umi_jsimon': begin
            directory = getenv('CALTECH')+'moogify/umi_jsimon/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'umi_jsimon_abund.eps'
            feh_p05 = -999
        end
        'sex1': begin
            directory = getenv('M31')+'dsph/sex/sex1/'
            vrlim = [190, 300]
            photfile = 'none'
            epsfile = 'sex1_abund.eps'
            feh_p05 = -999
        end
        'sex2': begin
            directory = getenv('M31')+'dsph/sex/sex2/'
            vrlim = [190, 300]
            photfile = 'none'
            epsfile = 'sex2_abund.eps'
            feh_p05 = -999
        end
        'sex3': begin
            directory = getenv('M31')+'dsph/sex/sex3/'
            vrlim = [190, 300]
            photfile = 'none'
            epsfile = 'sex3_abund.eps'
            feh_p05 = -999
        end
        'sex4': begin
            directory = getenv('M31')+'dsph/sex/sex4/'
            vrlim = [190, 300]
            photfile = 'none'
            epsfile = 'sex4_abund.eps'
            feh_p05 = -999
        end
        'sex6': begin
            directory = getenv('M31')+'dsph/sex/sex6/'
            vrlim = [190, 300]
            photfile = 'none'
            epsfile = 'sex6_abund.eps'
            feh_p05 = -999
        end
        'sexmi1': begin
            directory = getenv('CALTECH')+'moogify/sexmi1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sexmi1_abund.eps'
            feh_p05 = -999
        end
        'sexmi2': begin
            directory = getenv('CALTECH')+'moogify/sexmi2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sexmi2_abund.eps'
            feh_p05 = -999
        end
        'sexmi3': begin
            directory = getenv('CALTECH')+'moogify/sexmi3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sexmi3_abund.eps'
            feh_p05 = -999
        end
        'sexmi4': begin
            directory = getenv('CALTECH')+'moogify/sexmi4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sexmi4_abund.eps'
            feh_p05 = -999
        end
        'sexmi5': begin
            directory = getenv('CALTECH')+'moogify/sexmi5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sexmi5_abund.eps'
            feh_p05 = -999
        end
        'dra1': begin
            directory = getenv('M31')+'dsph/dra/dra1/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra1_abund.eps'
            feh_p05 = -999
        end
        'dra2': begin
            directory = getenv('M31')+'dsph/dra/dra2/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra2_abund.eps'
            feh_p05 = -999
        end
        'dra3': begin
            directory = getenv('M31')+'dsph/dra/dra3/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra3_abund.eps'
            feh_p05 = -999
        end
        'dra4': begin
            directory = getenv('M31')+'dsph/dra/dra4/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra4_abund.eps'
            feh_p05 = -999
        end
        'dra5': begin
            directory = getenv('M31')+'dsph/dra/dra5/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra5_abund.eps'
            feh_p05 = -999
        end
        'dra7': begin
            directory = getenv('M31')+'dsph/dra/dra7/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra7_abund.eps'
            feh_p05 = -999
        end
        'dra8': begin
            directory = getenv('M31')+'dsph/dra/dra8/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra8_abund.eps'
            feh_p05 = -999
        end
        'dra9': begin
            directory = getenv('M31')+'dsph/dra/dra9/'
            vrlim = [-350, -250]
            photfile = 'none'
            epsfile = 'dra9_abund.eps'
            feh_p05 = -999
        end
        'n6822a': begin
            directory = getenv('CALTECH')+'moogify/n6822a/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'n6822a_abund.eps'
            feh_p05 = -999
        end
        'n6822b': begin
            directory = getenv('CALTECH')+'moogify/n6822b/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'n6822b_abund.eps'
            feh_p05 = -999
        end
        'i1613a': begin
            directory = getenv('CALTECH')+'moogify/i1613a/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'i1613a_abund.eps'
            feh_p05 = -999
        end        
        'pega': begin
            directory = getenv('CALTECH')+'moogify/pega/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'pega_abund.eps'
            feh_p05 = -999
        end
        'aqra': begin
            directory = getenv('CALTECH')+'moogify/aqra/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'aqra_abund.eps'
            feh_p05 = -999
        end
        'ceta': begin
            directory = getenv('CALTECH')+'moogify/ceta/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'ceta_abund.eps'
            feh_p05 = -999
        end
        'cetb': begin
            directory = getenv('CALTECH')+'moogify/cetb/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'cetb_abund.eps'
            feh_p05 = -999
        end
        'aqrd': begin
            directory = getenv('CALTECH')+'moogify/aqrd/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'aqra_abund.eps'
            feh_p05 = -999
        end
        'leoaaw': begin
            directory = getenv('CALTECH')+'moogify/leoaaW/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoaaW_abund.eps'
            feh_p05 = -999
        end
        'leoac': begin
            directory = getenv('CALTECH')+'moogify/leoac/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoac_abund.eps'
            feh_p05 = -999
        end
        'leoa_rizzi': begin
            directory = getenv('CALTECH')+'moogify/leoa_rizzi/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'leoa_rizzi_abund.eps'
            feh_p05 = -999
        end
        'sagdia': begin
            directory = getenv('CALTECH')+'moogify/sagdia/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sagdia_abund.eps'
            feh_p05 = -999
        end
        'sagdib': begin
            directory = getenv('CALTECH')+'moogify/sagdib/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sagdib_abund.eps'
            feh_p05 = -999
        end
        'vv124a': begin
            directory = getenv('CALTECH')+'moogify/vv124a/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'vv124a_abund.eps'
            feh_p05 = -999
        end        
        'vv124b': begin
            directory = getenv('CALTECH')+'moogify/vv124b/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'vv124b_abund.eps'
            feh_p05 = -999
        end        
        'vv124ab': begin
            directory = getenv('CALTECH')+'moogify/vv124ab/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'vv124ab_abund.eps'
            feh_p05 = -999
        end        
        'callisto': begin
            directory = getenv('M31')+'callisto/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'callisto_abund.eps'
            feh_p05 = -999
        end
        'arcturus': begin
            directory = getenv('M31')+'arcturus/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'arcturus_abund.eps'
            feh_p05 = -999
        end
        '208bosb': begin
            directory = getenv('CALTECH')+'moogify/208BoSB/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = '208bosb_abund.eps'
            feh_p05 = -999
        end
        'echos1_2': begin
            directory = getenv('CALTECH')+'echos/echos1_2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'echos1_2_abund.eps'
            feh_p05 = -999
        end
        'echos2_3': begin
            directory = getenv('CALTECH')+'echos/echos2_3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'echos2_3_abund.eps'
            feh_p05 = -999
        end
        'echos3_1': begin
            directory = getenv('CALTECH')+'echos/echos3_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'echos3_1_abund.eps'
            feh_p05 = -1.62
        end
        'echos7_1': begin
            directory = getenv('CALTECH')+'echos/echos7_1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'echos7_1_abund.eps'
            feh_p05 = -999
        end
        'alphatest': begin
            directory = getenv('M31')+'specabund/alphatest/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'alphatest_abund.eps'
            feh_p05 = -999
        end
        'v368her': begin
            directory = getenv('chome')+'deimos/2010aug11/v368Her1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'alphatest_abund.eps'
            feh_p05 = -999
        end
        'sesar_her3': begin
            directory = getenv('CALTECH')+'moogify/sesar/Her3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sesar_Her3_abund.eps'
            feh_p05 = -999
        end
        'sesar_her4': begin
            directory = getenv('CALTECH')+'moogify/sesar/Her4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sesar_Her4_abund.eps'
            feh_p05 = -999
        end
        'sesar_herclump': begin
            directory = getenv('CALTECH')+'moogify/sesar/HerClump/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sesar_HerClump_abund.eps'
            feh_p05 = -999
        end
        'sesar_aqr': begin
            directory = getenv('CALTECH')+'moogify/sesar/Aqr/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sesar_Aqr_abund.eps'
            feh_p05 = -999
        end
        'sesar_lvmslitb': begin
            directory = getenv('CALTECH')+'moogify/sesar/LVMslitB/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'sesar_LVMslitB_abund.eps'
            feh_p05 = -999
        end
        'hilke': begin
            directory = getenv('CALTECH')+'hilke/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'alphatest_abund.eps'
            feh_p05 = -999
        end
        'cnc': begin
            directory = getenv('CALTECH')+'moogify/Cnc/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'Cnc_abund.eps'
            feh_p05 = -999
        end
        'vss1': begin
            directory = getenv('CALTECH')+'moogify/VSS1/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS1_abund.eps'
            feh_p05 = -999
        end
        'vss2': begin
            directory = getenv('CALTECH')+'moogify/VSS2/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS2_abund.eps'
            feh_p05 = -999
        end
        'vss3': begin
            directory = getenv('CALTECH')+'moogify/VSS3/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS3_abund.eps'
            feh_p05 = -999
        end
        'vss4': begin
            directory = getenv('CALTECH')+'moogify/VSS4/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS4_abund.eps'
            feh_p05 = -999
        end
        'vss5': begin
            directory = getenv('CALTECH')+'moogify/VSS5/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS5_abund.eps'
            feh_p05 = -999
        end
        'vss6': begin
            directory = getenv('CALTECH')+'moogify/VSS6/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS6_abund.eps'
            feh_p05 = -999
        end
        'vss7': begin
            directory = getenv('CALTECH')+'moogify/VSS7/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS7_abund.eps'
            feh_p05 = -999
        end
        'vss8': begin
            directory = getenv('CALTECH')+'moogify/VSS8/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS8_abund.eps'
            feh_p05 = -999
        end
        'vss9': begin
            directory = getenv('CALTECH')+'moogify/VSS9/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'VSS9_abund.eps'
            feh_p05 = -999
        end
        'scatsc_lr': begin
            directory = getenv('UCI')+'ucsc_workshop/sc_at_sc_lr/'
            vrlim = [-1000, 1000]
            photfile = 'none'
            epsfile = 'scatsc_lr_abund.eps'
            feh_p05 = -999
        end
    endcase
end
