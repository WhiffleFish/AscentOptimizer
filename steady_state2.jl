using Ipopt
using JuMP

Δt = 0.1
N = 30
M = 100
g = 9.81

## Initial Conditions
x0 = 0
y0 = 500
Vx0 = 300
Vy0 = 0
θ0 = 0

## Constraints
target_y = 500
ϵ_y = 10
target_Vx = 300
ϵ_Vx = 10
target_Vy = 0.0
ϵ_Vy = 5
Tmax = 5_000
max_rot = deg2rad(20)/Δt
model = JuMP.Model(Ipopt.Optimizer)

# ρ(h) = 1/(h^2+1)

ρ(h)=0.0

function Fd(h, Vx, Vy)
    return 0.0
end

function Fdx(h, Vx, Vy)
    0.5*ρ(h)*sqrt(Vx^2 + Vy^2)*Vx
end

function Fdy(h, Vx, Vy)
    0.5*ρ(h)*sqrt(Vx^2 + Vy^2)*Vy
end

function Fx(h, T, θ, Vx, Vy)
    T*cos(θ) - Fdx(h,Vx,Vy)
end

function Fy(h, T, θ, Vx, Vy)
    T*sin(θ) - Fdy(h,Vx,Vy) - M*g
end

JuMP.@variables model begin
    0 ≤ T[1:N] ≤ Tmax
    -max_rot ≤ θ_dot[1:N] ≤ max_rot
    0 ≤ θ[1:N] ≤ π/2
    x[1:N]
    Vx[1:N]
    0 ≤ y[1:N]
    Vy[1:N]
end

JuMP.@constraint(model, x[1] == x0)
JuMP.@constraint(model, Vx[1] == Vx0)
JuMP.@constraint(model, target_Vx - ϵ_Vx ≤ Vx[N] ≤ target_Vx + ϵ_Vx)
JuMP.@constraint(model, y[1] == y0)
JuMP.@constraint(model, target_y - ϵ_y ≤ y[N] ≤ target_y + ϵ_y)
JuMP.@constraint(model, Vy[1] == Vy0)
JuMP.@constraint(model, target_Vy - ϵ_Vy ≤ Vy[N] ≤ target_Vy + ϵ_Vy)
# JuMP.@constraint(model, θ[1] == θ0)

# Dynamic Constraints
JuMP.@NLconstraint(model,
    [i=2:N], x[i] == x[i-1] + Vx[i-1]Δt + 0.5*(Fy(y[i-1], T[i-1],θ[i-1],Vx[i-1],Vy[i-1])/M)*Δt^2
)

JuMP.@NLconstraint(model,
    [i=2:N], Vx[i] == Vx[i-1] + Fx(y[i-1], T[i-1],θ[i-1],Vx[i-1],Vy[i-1])*Δt/M
)

JuMP.@NLconstraint(model,
    [i=2:N], y[i] == y[i-1] + Vy[i-1]Δt + 0.5*(Fy(y[i-1], T[i-1],θ[i-1],Vx[i-1],Vy[i-1])/M)*Δt^2
)

JuMP.@NLconstraint(model,
    [i=2:N], Vy[i] == Vy[i-1] + Fy(y[i-1], T[i-1],θ[i-1],Vx[i-1],Vy[i-1])*Δt/M
)

JuMP.@constraint(model,
    [i=2:N], θ[i] == θ[i-1] + θ_dot[i-1]*Δt
)

JuMP.@objective(model, Min,
    sum(((Vx ./ target_Vx) .- 1)).^2 + sum(((y ./ target_y) .- 1)).^2 + sum(((T ./ Tmax) .- 1)).^2 )

optimize!(model)

using Plots
plot(value.(x), value.(y), label= "" )
plot(value.(Vx), label="")
plot(value.(Vy), label="")
plot(value.(T), label="")
plot(value.(θ) .|> rad2deg, label="")
