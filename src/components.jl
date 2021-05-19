abstract type Component end

include("./components/free.jl")


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
