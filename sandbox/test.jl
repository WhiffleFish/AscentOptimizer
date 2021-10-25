using CSV, DataFrames
using Plots

df_e = DataFrame(CSV.File("earth_data.csv", delim=','))
alt_e = df_e[!,"Altitude"]
density_e = df_e[!,"Density"]


df_k = DataFrame(CSV.File("kerbin_data.csv"))
alt_k = df_k[!,"Altitude(m)"]
density_k = df_k[!,"Density(kg/m^3)"]



using LsqFit
m(t,p) = p[1]*exp.(p[2]*t)
p0 = [0.5,-0.001]
fit = curve_fit(m, alt_e, density_e, p0)

plot(alt_e, density_e)
Plots.plot!(alt_e, m(alt_e,fit.param))
