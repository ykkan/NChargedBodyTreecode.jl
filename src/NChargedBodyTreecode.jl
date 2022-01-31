module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra
    using BenchmarkTools
    using TimerOutputs

    export Particles
    export ClusterTree, num_cluster

    export updateParticlesField!
    export BruteForce, TreecodeUnstretch, TreecodeStretch, TreecodeAvgRestFrame

    export boostParticlesPosition!
    export transformParticlesMomentums!
    export transtransformParticlesField!

    export relerror

    include("algorithm.jl")
    include("particles.jl")
    include("cluster.jl")
    include("cluster_tree.jl")
    include("cluster2p.jl")
    include("update_particles_field.jl")
    include("force_kernel.jl")
    include("relativity.jl")
    include("diagnostic.jl")
end # module
