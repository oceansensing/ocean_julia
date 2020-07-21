module EOS
using GSW

export rho_s, rho_c

function rho_s(temp::Float64, salt::Float64, pres::Float64, lon= -75.0, lat=37.7)
    saltA = GSW.gsw_sa_from_sp(salt, pres, lon, lat)
    ctmp = GSW.gsw_ct_from_t(saltA, temp, pres)
    rho = GSW.gsw_rho(saltA, ctmp, pres)
    return rho
end

function rho_c(temp::Float64, cond::Float64, pres::Float64, lon= -75.0, lat=37.7)
    salt = GSW.gsw_sp_from_c(cond, temp, pres)
    saltA = GSW.gsw_sa_from_sp(salt, pres, lon, lat)
    ctmp = GSW.gsw_ct_from_t(saltA, temp, pres)
    rho = GSW.gsw_rho(saltA, ctmp, pres)
    return rho
end

end