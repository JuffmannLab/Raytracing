module Raytracing

# export the different functions and structs that can be used
export Free, Lens, PondInteraction, Setup, propagation!
export Electron, getintensity, mcp
export loadintensity, LightField

# define some nature constants
global const m_e = 9.1093837015e-31  # the electron mass
global const q = 1.602176634e-19     # electron charge
global const c = 299792458           # the speed of light
global const ε_0 = 8.8541878128e-12  # vacuum permitivity

# include all the code that is used
include("./helperFunctions.jl")
include("./electron.jl")
include("./lightField.jl")
include("./components.jl")

end # module
