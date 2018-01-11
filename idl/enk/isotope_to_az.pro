function isotope_to_az, isotope
    case isotope of
        'p': begin
            mass = 1
            Z = 1
        end
        'd': begin
            mass = 2
            Z = 1
        end
        else: begin
            isostring = strsplit(isotope, '^', /extract)
            mass = fix(isostring[0])
            case strtrim(isostring[1], 2) of
                'He': Z = 2
                'Li': Z = 3
                'Be': Z = 4
                'B': Z = 5
                'C': Z = 6
                'N': Z = 7
                'O': Z = 8
                'F': Z = 9
                'Ne': Z = 10
                'Na': Z = 11
                'Mg': Z = 12
                'Al': Z = 13
                'Si': Z = 14
                'P': Z = 15
                'S': Z = 16
                'Cl': Z = 17
                'Ar': Z = 18
                'K': Z = 19
                'Ca': Z = 20
                'Sc': Z = 21
                'Ti': Z = 22
                'V': Z = 23
                'Cr': Z = 24
                'Mn': Z = 25
                'Fe': Z = 26
                'Co': Z = 27
                'Ni': Z = 28
                'Cu': Z = 29
                'Zn': Z = 30
                'Ga': Z = 31
                'Ge': Z = 32
            endcase
        end
    endcase
    return, [mass, Z]
end
