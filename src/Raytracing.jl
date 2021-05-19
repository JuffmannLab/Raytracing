module Raytracing

export Free, Lense, Setup, propagation!
export Wave

abstract type AbstractWave end

include("electronBeam.jl")
include("components.jl")

end # module
