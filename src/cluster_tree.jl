struct ClusterTree{T}
    npar::Int64
    parindices::Vector{Int64}
    root::Cluster{T}
end

function ClusterTree(particles::Particles; n, N0, stretch=SVector(1.0,1.0,1.0)) where {T}
    npar = particles.npar
    parindices = [1:npar;]
    root = subdivide(particles, parindices, 1, npar, 1; n=n, N0=N0, stretch=stretch)
    return ClusterTree(npar, parindices, root)
end

function print_bbox(ct::ClusterTree{T}) where {T}
    print_bbox(ct.root)
end
