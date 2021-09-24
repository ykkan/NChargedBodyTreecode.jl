module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools
    using TimerOutputs

    export Particles
    export ClusterTree, num_cluster
    export update_particle_field!, update_particle_field_brutal!

    include("particles.jl")
    include("cluster.jl")
    include("update_particles_field.jl")
end # module
