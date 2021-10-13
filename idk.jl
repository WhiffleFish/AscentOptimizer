const G = 6.67408e-11

struct Planet
    R::Float64
    M::Float64
    μ::Float64
    data
end
Earth() = Planet(6378.137e3,5.9724e24,G*5.9724e24)

struct Rocket
    M::Float64
    Cd::Float64
    Tmax::Float64
    max_rot::Float64
end
