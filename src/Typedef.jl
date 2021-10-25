using JuMP
const G = 6.67408e-11

struct Planet{Fρ<:Function}
    R::Float64
    M::Float64
    μ::Float64
    g0::Float64
    ρ::Fρ
    data::DataFrame
end

struct Rocket
    M::Float64 # kg
    Cd::Float64 # unitless
    Tmax::Float64 # N
    max_rot::Float64 # rad/s
end

Rocket() = Rocket(100, 1.0, 500, deg2rad(20))

mutable struct Targets
    y::Float64
    ϵ_y::Float64
    Vx::Float64
    ϵ_Vx::Float64
    Vy::Float64
    ϵ_Vy::Float64
end

function Targets(p::Planet, h::Float64)
    r = h + p.R
    Vx = √(p.μ/r)
    return Targets(
        h,
        100.0,
        Vx,
        10.0,
        0.0,
        10.0
    )
end

mutable struct OptParams
    model::JuMP.Model
    targets::Targets
    planet::Planet
    rocket::Rocket
    Δt::Float64
    N::Int
end
