module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools
    using TimerOutputs

    export Particles
    export ClusterTree, num_cluster

    export BruteForce, TreecodeUnstretch, TreecodeStretch, TreecodeAvgRestFrame
    export updateParticlesField!

    export transformParticlesDistribution!, transformParticlesField!

    include("algorithm.jl")
    include("particles.jl")
    include("cluster.jl")
    include("cluster_tree.jl")
    include("cluster2p.jl")
    include("update_particles_field.jl")
    include("force_kernel.jl")
    include("relativity.jl")
end # module
