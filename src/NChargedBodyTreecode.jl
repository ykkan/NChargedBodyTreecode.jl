module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools
    using TimerOutputs

    export Particles
    export ClusterTree, num_cluster
    export update_particle_field!, update_particle_field_brutal!
    export update_particle_field_stretch!, update_particle_field_AVGRF!
    export particle2frame!, lab2avgframe!, avg2labframe!, field2labframe, field2avgframe
    export field2frame, fields2frame!
    export cluster2p

    include("particles.jl")
    include("cluster.jl")
    include("cluster_tree.jl")
    include("cluster2p.jl")
    include("update_particles_field.jl")
    include("force_kernel.jl")
    include("relativity.jl")
end # module
