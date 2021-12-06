struct ClusterTree{T}
    npar::Int64
    parindices::Vector{Int64}
    root::Cluster{T}
end

function ClusterTree(particles::Particles; n, threshold, stretch=SVector(1.0,1.0,1.0)) where {T}
    npar = particles.npar
    parindices = [1:npar;]
    root = make_cluster(particles, parindices, 1, npar, 1; n=n, threshold=threshold, stretch=stretch)
    return ClusterTree(npar, parindices, root)
end

function print_bbox(ct::ClusterTree{T}) where {T}
    print_bbox(ct.root)
end
