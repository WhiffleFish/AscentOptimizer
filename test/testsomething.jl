using JuMP, Ipopt
using AscentOptimizer

planet = Kerbin
rocket = Rocket(100,1.0,5_000,deg2rad(20))
targets = Targets(
    10_000,
    10,
    500,
    10,
    0,
    5,
)

model = JuMP.Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "max_iter", 10_000)
params = OptParams(
    model,
    targets,
    planet,
    rocket,
    1.0, # Δt
    200 # N
)

instantiate!(params)

optimize!(params)

using Plots
fig = plot(params)

plot(value.(model[:θ_dot]))
