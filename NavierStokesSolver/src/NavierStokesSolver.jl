module NavierStokesSolver

include("core.jl")
include("euler.jl")
include("rk2.jl")
include("rk4.jl")
include("pressure.jl")
include("boundary.jl")
include("diagnostic.jl")
include("vizualization.jl")

export navier_stokes_solver, animate_flow

end
