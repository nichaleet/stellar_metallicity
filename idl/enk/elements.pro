function elements
    name = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', $
            'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', $
            'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', $
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', $
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', $
            'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', $
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', $
            'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', $
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', $
            'Pa', 'U',  'Np', 'Pu', 'Am'];, 'Cm', 'Bk', 'Cf', 'Es', 'Fm', $
           ;'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', $
           ;'Rg']

    solar =  [12.00,10.99,3.31, 1.42, 2.88, 8.56, 8.05, 8.93, 4.56, 8.09, $    
              6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.50, 6.56, 5.12, 6.36, $
              3.10, 4.99, 4.00, 5.67, 5.39, 7.52, 4.92, 6.25, 4.21, 4.60, $
              2.88, 3.41, 2.37, 3.35, 2.63, 3.23, 2.60, 2.90, 2.24, 2.60, $
              1.42, 1.92, 0.00, 1.84, 1.12, 1.69, 1.24, 1.86, 0.82, 2.00, $
              1.04, 2.24, 1.51, 2.23, 1.12, 2.13, 1.22, 1.55, 0.71, 1.50, $
              0.00, 1.00, 0.51, 1.12, 0.33, 1.10, 0.50, 0.93, 0.13, 1.08, $    
              0.12, 0.88, 0.13, 0.68, 0.27, 1.45, 1.35, 1.80, 0.83, 1.09, $
              0.82, 1.85, 0.71, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12, $   
              0.00, 0.00, 0.00, 0.00, 0.00]

    xam = [1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, $
           22.99, 24.31, 26.98, 28.08, 30.97, 32.06, 35.45, 39.95, 39.10, 40.08, $
           44.96, 47.90, 50.94, 52.00, 54.94, 55.85, 58.93, 58.71, 63.55, 65.37, $     
           69.72, 72.59, 74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22, $
           92.91, 95.94, 98.91, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, $     
           121.8, 127.6, 126.9, 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, $
           145.0, 150.4, 152.0, 157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, $     
           175.0, 178.5, 181.0, 183.9, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, $
           204.4, 207.2, 209.0, 210.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, $
           231.0, 238.0, 237.0, 244.0, 243.0]                               

    xchi1 = [13.60, 24.59, 5.392, 9.323, 8.298, 11.26, 14.53, 13.62, 17.42, 21.56, $
             5.139, 7.646, 5.986, 8.152, 10.49, 10.36, 12.97, 15.76, 4.340, 6.113, $
             6.562, 6.828, 6.746, 6.766, 7.434, 7.902, 7.881, 7.640, 7.726, 9.394, $
             5.999, 7.899, 9.789, 9.752, 11.81, 14.00, 4.177, 5.695, 6.217, 6.634, $     
             6.579, 7.092,  7.28, 7.361, 7.459, 8.337, 7.576, 8.994, 5.786, 7.344, $     
             8.608, 9.010, 10.45, 12.13, 3.894, 5.212, 5.577, 5.539, 5.473, 5.525, $     
             5.582, 5.644, 5.670, 6.150, 5.864, 5.939, 6.022, 6.108, 6.184, 6.254, $     
             5.426, 6.825, 7.550, 7.864, 7.834, 8.348, 8.967, 8.959, 9.226, 10.44, $     
             6.108, 7.417, 7.286, 8.417,   9.0, 10.75, 4.073, 5.278,  5.17, 6.307, $     
              5.89, 6.194, 6.266, 6.026, 5.974]                                   

    xchi2 = [50.00, 54.42, 75.64, 18.21, 25.16, 24.38, 29.60, 35.12, 34.97, 40.96, $
             47.29, 15.04, 18.83, 16.35, 19.77, 23.34, 23.81, 27.63, 31.63, 11.87, $
             12.80, 13.58, 14.66, 16.50, 15.64, 16.19, 17.08, 18.17, 20.29, 17.96, $
             20.52, 15.94, 18.59, 21.16, 21.81, 24.36, 27.29, 11.03, 12.22, 13.13, $
             14.32, 16.16, 15.26, 16.76, 18.08, 19.43, 21.49, 16.91, 18.87, 14.63, $
             16.53,  18.6, 19.13, 20.98, 23.16, 10.00,  11.1,  10.8,  10.6,  10.7, $
              10.9,  11.1, 11.24,  12.1,  11.5,  11.7,  11.8,  11.9,  12.1, 12.18, $
              13.9,  14.9,  16.2,  17.7,  16.6,  17.0,   30., 18.56,  20.5, 18.76, $
             20.43, 15.03, 16.70,  19.0,   30.,   30.,  22.0, 10.15, 11.75,  11.9, $
              12.0,  11.9,   30.,   30.,   30.]

    xchi3 = [90.00, 90.00, 122.5, 153.9, 37.93, 47.89, 47.45, 54.94, 62.71, 63.46, $
             71.62, 80.14, 28.45, 33.49, 30.20, 34.83, 39.61, 40.91, 45.81, 50.91, $
             24.76, 27.49, 29.31,  31.0, 33.67, 30.65,  33.5,  35.3, 36.84, 39.72, $
             30.73,  34.2,  28.4, 30.82,  35.9, 36.95,  39.2, 42.88, 20.53,  23.1, $
              25.0,  27.2,  29.5,  28.5,  31.1,  32.9,  34.8, 37.47,  28.0, 30.50, $
             25.32, 27.96,  33.0,  31.0,  33.4,  35.8, 19.18, 20.20, 21.62,  22.1, $
              22.3,  23.4,  24.9,  20.6,  21.9,  22.8,  22.8,  22.7,  23.7, 25.05, $
             20.96,  23.3,   60.,   60.,   60.,   60.,   60.,   60.,   34.,  34.2, $
             29.85, 31.94, 25.56,   60.,   60.,   60.,   60.,   60.,   20.,  18.3, $
               60.,   20.,   60.,   60.,  60.]

    n = n_elements(name)
    elem = {name:' ', atomic:0L, solar:0d, xam:0d, xchi1:0d, xchi2:0d, xchi3:0d}
    elem = replicate(elem, n)
    elem.name = name
    elem.atomic = lindgen(n)+1
    elem.solar = solar
    elem.xam = xam
    elem.xchi1 = xchi1
    elem.xchi2 = xchi2
    elem.xchi3 = xchi3
    return, elem
end
