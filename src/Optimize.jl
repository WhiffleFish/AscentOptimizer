import JuMP.optimize!

#=
function ρ(planet::Planet, h)
    planet.ρ(h)
end

function ∇ρ(ρ, Cd, V) end

function Fd(ρ, Cd, V)
    0.5*ρ*V^2*Cd
end

function ∇Fd(ρ, Cd, V) end
=#

function dx end

function dy end

function dVx end

function dVy end

function Fd(ρ, Vx, Vy)
    ρ*(Vx^2 + Vy^2)
end

function Fd(ρ, V)
    ρ*V^2
end

function instantiate!(params::OptParams)
    model = params.model
    N = params.N
    planet = params.planet
    ρ = planet.ρ
    g0 = planet.g0
    targets = params.targets
    Δt = params.Δt

    rocket = params.rocket
    M = rocket.M
    Cd = rocket.Cd
    Tmax = rocket.Tmax
    max_rot = rocket.max_rot*Δt

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
    JuMP.@constraint(model, targets.Vx - targets.ϵ_Vx ≤ Vx[N] ≤ targets.Vx + targets.ϵ_Vx)
    JuMP.@constraint(model, y[1] == 0.0)
    JuMP.@constraint(model, targets.y - targets.ϵ_y ≤ y[N] ≤ targets.y + targets.ϵ_y)
    JuMP.@constraint(model, Vy[1] == 0.0)
    JuMP.@constraint(model, targets.Vy - targets.ϵ_Vy ≤ Vy[N] ≤ targets.Vy + targets.ϵ_Vy)
    JuMP.@constraint(model, θ[1] == π/2)

    JuMP.@NLconstraint(model,
        [i=2:N], x[i] == x[i-1] + Vx[i-1]Δt + 0.5*(T[i-1] - Fd(0.01*ρ(y[i-1]), Vx[i-1], Vy[i-1]))*cos(θ[i-1])*Δt^2
    )

    JuMP.@NLconstraint(model,
        [i=2:N], Vx[i] == Vx[i-1] + ((T[i-1] - Fd(0.01*ρ(y[i-1]), Vx[i-1], Vy[i-1]))*cos(θ[i-1])/M)*Δt
    )

    JuMP.@NLconstraint(model,
        [i=2:N], y[i] == y[i-1] + Vy[i-1]Δt + (1/(2*M))*((T[i-1] - Fd(0.01*ρ(y[i-1]), Vx[i-1], Vy[i-1]))*sin(θ[i-1]) - M*g0)*Δt^2
    )

    JuMP.@NLconstraint(model,
        [i=2:N], Vy[i] == Vy[i-1] + ((T[i-1] - Fd(0.01*ρ(y[i-1]), Vx[i-1], Vy[i-1]))*sin(θ[i-1])/M - g0)*Δt
    )

    JuMP.@constraint(model,
        [i=2:N], θ[i] == θ[i-1] + θ_dot[i-1]*Δt
    )

    JuMP.@objective(model, Min,
        sum(T.^2)
    )
end

#=
function instantiate!(params::OptParams)
    model = params.model
    N = params.N
    planet = params.planet
    ρ = planet.ρ
    g0 = planet.g0
    targets = params.targets
    Δt = params.Δt

    rocket = params.rocket
    M = rocket.M
    Cd = rocket.Cd
    Tmax = rocket.Tmax
    max_rot = rocket.max_rot

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
    JuMP.@constraint(model, targets.Vx - targets.ϵ_Vx ≤ Vx[N] ≤ targets.Vx + targets.ϵ_Vx)
    JuMP.@constraint(model, y[1] == 0.0)
    JuMP.@constraint(model, targets.y - targets.ϵ_y ≤ y[N] ≤ targets.y + targets.ϵ_y)
    JuMP.@constraint(model, Vy[1] == 0.0)
    JuMP.@constraint(model, targets.Vy - targets.ϵ_Vy ≤ Vy[N] ≤ targets.Vy + targets.ϵ_Vy)
    JuMP.@constraint(model, θ[1] == π/2)

    JuMP.@NLconstraint(model,
        [i=2:N],
            x[i] == x[i-1]
                + Vx[i-1]Δt
                + 0.5*(
                    T[i-1]*cos(θ[i-1])
                    - Fd(ρ(y[i-1]),Cd,Vx[i-1])*sign(Vx[i-1])
                    )
                *Δt^2/M
        )
        # begin
        # _fdx = Fd(ρ(y[i-1]),Cd,Vx[i-1])*sign(Vx[i-1])
        # ∑Fx = T[i-1]*cos(θ[i-1]) - _fdx
        # x[i] == x[i-1] + Vx[i-1]Δt + 0.5*(∑Fx/M)*Δt^2
        # end


    JuMP.@NLconstraint(model,
        [i=2:N],
            Vx[i] == Vx[i-1]
                +(
                    T[i-1]*cos(θ[i-1])
                    #- Fd(ρ(y[i-1]),Cd,Vx[i-1])*sign(Vx[i-1])
                )
                *Δt/M
        )
        # _fdx = Fd(ρ(y[i-1]),Cd,Vx[i-1])*sign(Vx[i-1])
        # ∑Fx = T[i-1]*cos(θ[i-1]) - _fdx
        # Vx[i] == Vx[i-1] + (∑Fx/M)*Δt


    JuMP.@NLconstraint(model,
        [i=2:N],
            y[i] == y[i-1]
                + Vy[i-1]Δt
                + (
                    T[i-1]*sin(θ[i-1])
                    #- Fd(ρ(y[i-1]),Cd,Vy[i-1])*sign(Vy[i-1])
                    - M*g0
                )/M
                *Δt^2
        )
        # begin
        # _fdy = Fd(ρ(y[i-1]),Cd,Vy[i-1])*sign(Vy[i-1])
        # ∑Fy = T[i-1]*sin(θ[i-1]) - _fdy - M*g0
        # y[i] == y[i-1] + Vy[i-1]Δt + (∑Fy/M)*Δt^2
        # end


    JuMP.@NLconstraint(model,
        [i=2:N],
            Vy[i] == Vy[i-1]
                + (
                    T[i-1]*sin(θ[i-1])
                    #- Fd(ρ(y[i-1]),Cd,Vy[i-1])*sign(Vy[i-1])
                    - M*g0
                )
                *Δt/M
        )
        # begin
        # _fdy = Fd(ρ(y[i-1]),Cd,Vy[i-1])*sign(Vy[i-1])
        # ∑Fy = T[i-1]*sin(θ[i-1]) - _fdy - M*g0
        # Vy[i] == Vy[i-1] + (∑Fy/M)*Δt
        # end

    JuMP.@constraint(model,
        [i=2:N], θ[i] == θ[i-1] + θ_dot[i-1]*Δt
    )

    JuMP.@objective(model, Min,
        sum(T.^2)
    )
end
=#

JuMP.optimize!(params::OptParams) = optimize!(params.model)
