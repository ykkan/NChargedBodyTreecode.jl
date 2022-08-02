module NChargedBodyTreecode
    using StaticArrays
    using LinearAlgebra

    export Particles, Beam
    export MacroParticles
    export ClusterTree

    export updateParticlesField!
    export BruteForce, TreecodeUniform, TreecodeStretch, TreecodeAvgRestFrame

    export boostParticlesPosition!
    export transformParticlesMomentums!
    export transtransformParticlesField!

    export relerror

    include("utils.jl")
    include("algorithm.jl")
    include("particles.jl")
    include("macroparticles.jl")
    include("cluster.jl")
    include("cluster_tree.jl")
    include("cluster2p.jl")
    include("upwardpass.jl")
    include("update_particles_field.jl")
    include("force_kernel.jl")
    include("relativity.jl")
    include("diagnostic.jl")
end # module
