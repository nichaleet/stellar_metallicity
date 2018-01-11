function ReadSpec, file
    bs = mrdfits(file,1, /silent)
    rs = mrdfits(file,2, /silent)
    ss = {ss, lambda:[bs.lambda, rs.lambda], spec:[bs.spec, rs.spec], $
          ivar:[bs.ivar, rs.ivar]}
    return, ss
end
