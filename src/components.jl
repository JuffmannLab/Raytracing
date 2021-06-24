
abstract type Component end

include("./components/free.jl")
include("./components/lense.jl")
include("./components/pondInteraction.jl")

struct Setup
    setup::Vector{<:Component}
end

"""
    Setup(comps::Component...)::Setup

Return the Setup struct.

Create and return the setup struct. The varargs `comps` describe
the different components that are used in the setup.
"""
Setup(comps::Component...) = Setup(collect(comps))

"""
    propagation!(rays::AbstractElectron, setup::Setup)

Propagate the electrons in the setup.

Take the electrons in `rays` and apply all the different alterations
that are carried out by the setup denoted by `setup`.
"""
function propagation!(rays::AbstractElectron, setup::Setup)
    
    # iterate over all the components, apply the calculate function
    for component in setup.setup
        calculate!(rays, component)
    end
end
