using LsqFit

##
m(t,p) = p[1]*exp.(p[2]*t)

const Earth = begin
    R = 6378.137e3
    M = 5.9724e24
    g0 = G*M/R^2
    data = DataFrame(CSV.File(joinpath(@__DIR__,"data","earth_data.csv")))
    alt = data[!,"Altitude"]
    density = data[!,"Density"]
    p0 = [0.5,-0.001]
    fit = curve_fit(m, alt, density, p0)
    const EARTH_DENSITY_PARAMS = fit.param
    fρ = h -> m(h, EARTH_DENSITY_PARAMS)
    Planet(
        R,
        M,
        G*M,
        g0,
        fρ,
        data
) end

const Kerbin = begin
    R = 600.000e3
    M = 5.2915158e22
    g0 = G*M/R^2
    data = DataFrame(CSV.File(joinpath(@__DIR__,"data","kerbin_data.csv")))
    alt = data[!,"Altitude(m)"]
    density = data[!,"Density(kg/m^3)"]
    p0 = [0.5,-0.001]
    fit = curve_fit(m, alt, density, p0)
    const KERBIN_DENSITY_PARAMS = fit.param
    fρ = h -> m(h, KERBIN_DENSITY_PARAMS)
    Planet(
        R,
        M,
        G*M,
        g0,
        fρ,
        data
) end
