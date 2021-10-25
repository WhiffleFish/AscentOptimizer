module AscentOptimizer

using JuMP
using DataFrames, CSV

include("Typedef.jl")
export Planet, Rocket, Targets, OptParams

include("Planets.jl")
export Earth, Kerbin

include("Optimize.jl")
export instantiate!

include("Plot.jl")

end # module
