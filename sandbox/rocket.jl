using Ipopt
using JuMP

Δt = 1.0
N = 200
M = 100
g = 9.81

function orbit_speed(r)
    sqrt(9.81*r)
end

target_y = 10_000
ϵ_y = 10
target_Vx = 500# orbit_speed(target_y)
ϵ_Vx = 10
target_Vy = 0.0
ϵ_Vy = 5
Tmax = 5_000
max_rot = deg2rad(20)
model = JuMP.Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 10_000)

ρ(h) = 0.01*1.25*exp(-0.00014270586685126693*h)

function Fd(h, Vx, Vy)
    ρ(h)*(Vx^2 + Vy^2)
end

function Fd(h, V)
    ρ(h)*V^2
end

JuMP.@variables model begin
    0 ≤ T[1:N] ≤ Tmax
    -max_rot ≤ θ_dot[1:N] ≤ max_rot
    -π/2 ≤ θ[1:N] ≤ π/2
    x[1:N]
    Vx[1:N]
    0 ≤ y[1:N]
    Vy[1:N]
end

JuMP.@constraint(model, x[1] == 0.0)
JuMP.@constraint(model, Vx[1] == 0.0)
JuMP.@constraint(model, target_Vx - ϵ_Vx ≤ Vx[N] ≤ target_Vx + ϵ_Vx)
JuMP.@constraint(model, y[1] == 0.0)
JuMP.@constraint(model, target_y - ϵ_y ≤ y[N] ≤ target_y + ϵ_y)
JuMP.@constraint(model, Vy[1] == 0.0)
JuMP.@constraint(model, target_Vy - ϵ_Vy ≤ Vy[N] ≤ target_Vy + ϵ_Vy)
JuMP.@constraint(model, θ[1] == π/2)

# Dynamic Constraints
JuMP.@NLconstraint(model,
    [i=2:N], x[i] == x[i-1] + Vx[i-1]Δt + 0.5*(T[i-1] - Fd(y[i-1], Vx[i-1], Vy[i-1]))*cos(θ[i-1])*Δt^2
)

JuMP.@NLconstraint(model,
    [i=2:N], Vx[i] == Vx[i-1] + ((T[i-1] - Fd(y[i-1], Vx[i-1], Vy[i-1]))*cos(θ[i-1])/M)*Δt
)

JuMP.@NLconstraint(model,
    [i=2:N], y[i] == y[i-1] + Vy[i-1]Δt + (1/(2*M))*((T[i-1] - Fd(y[i-1], Vx[i-1], Vy[i-1]))*sin(θ[i-1]) - M*g)*Δt^2
)

JuMP.@NLconstraint(model,
    [i=2:N], Vy[i] == Vy[i-1] + ((T[i-1] - Fd(y[i-1], Vx[i-1], Vy[i-1]))*sin(θ[i-1])/M - g)*Δt
)

JuMP.@constraint(model,
    [i=2:N], θ[i] == θ[i-1] + θ_dot[i-1]*Δt
)

JuMP.@objective(model, Min,
    sum(T.^2)
)

optimize!(model)

##
using Plots

using CairoMakie

plot(value.(x), value.(y), label= "" )
plot(value.(Vx), label="")
plot(value.(Vy), label="")
plot(value.(T), label="")
plot(value.(θ) .|> rad2deg, label="")


##
using LaTeXStrings
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

fig

save("AscentProfile.svg", fig)

##
f = Figure()

ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label",
    title = "Title")

f
