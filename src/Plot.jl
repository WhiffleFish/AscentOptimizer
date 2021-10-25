using CairoMakie
using LaTeXStrings
using Plots
import Plots.plot

function Plots.plot(params::OptParams)
    Δt = params.Δt
    N = params.N
    
    model = params.model
    x = model[:x]
    y = model[:y]
    Vx = model[:Vx]
    Vy = model[:Vy]
    θ = model[:θ]
    T = model[:T]

    ts = 0:Δt:Δt*(N-1)
    length(values.(x))

    fig = Figure()
    ax1 = Axis(fig[1,1:2],title="Ascent Profile", xlabel=L"X [m]", ylabel=L"Altitude [$m$]")
    lines!(ax1, value.(x), value.(y))
    ax2 = Axis(fig[2,1], xlabel=L"Time [$s$]", ylabel=L"V_x")
    lines!(ax2, ts, value.(Vx))
    ax3 = Axis(fig[2,2], xlabel=L"Time [$s$]", ylabel=L"V_y")
    lines!(ax3, ts, value.(Vy))
    ax4 = Axis(fig[3,1], xlabel = L"Time [$s$]", ylabel = L"Thrust [$N$]")
    lines!(ax4, ts, value.(T))
    ax5 = Axis(fig[3,2], xlabel = L"Time [$s$]", ylabel = L"$\theta$ [deg]")
    lines!(ax5, ts, rad2deg.(value.(θ)))

    return fig
end
