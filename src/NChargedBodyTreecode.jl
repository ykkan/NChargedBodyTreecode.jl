module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools

    export Particles
    export update_particle_field!, update_particle_field_brutal!

    include("particles.jl")
    include("cluster.jl")
    include("update_particles_field.jl")
end # module
