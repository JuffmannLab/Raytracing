module Raytracing

# export the different functions and structs that can be used
export Free, Lense, PondInteraction, Setup, propagation!
export Electron, getintensity
export loadintensity, LightField

# define some nature constants
global const m_e = 9.1093837015e-31  # the electron mass
global const q = 1.602176634e-19     # electron charge
global const c = 299792458           # the speed of light

# include all the code that is used
include("./electron.jl")
include("./lightField.jl")
include("./components.jl")

end # module
