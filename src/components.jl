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


function propagation!(wave::Wave, setup::Setup)
    
    # create the setup matrix that contains all the 
    setup_mat = setup.setup[end].mat

    for i = size(setup.setup, 1)-1:-1:1
        setup_mat *= setup.setup[i].mat
    end

    for i = 1:size(wave.x, 1)
        wave.ψ[i] *= setup_mat
    end
end
