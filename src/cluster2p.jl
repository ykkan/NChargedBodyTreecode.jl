function cluster2p(x::SVector{3,T}, cluster::Cluster{T}, particles::Particles, parindices::Vector{Int64}; p_kernel, n, eta, stretch=SVector(1.0,1.0,1.0)) where {T}
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    if cluster.children === nothing
        efield, bfield = cluster2p_brutal(x, cluster, particles, parindices)
    elseif admissible(x, cluster.bbox; eta=eta, stretch=stretch)
        xcoords = cluster.xcoords
        ycoords = cluster.ycoords
        zcoords = cluster.zcoords
        gamma_hat = cluster.gammas
        mom_hat = cluster.momenta
        for k in 1:(n+1)
            sz = zcoords[k]
            for j in 1:(n+1)
                sy = ycoords[j]
                for i in 1:(n+1)
                    sx = xcoords[i]
                    p_hat = mom_hat[i,j,k]
                    K = kernel_relativity(x, SVector(sx,sy,sz), p_kernel)
                    efield += gamma_hat[i,j,k]*K
                    bfield += cross(p_hat, K)
                end
            end
        end
    else
        efield1, bfield1 = cluster2p(x, cluster.children[1], particles, parindices; p_kernel=p_kernel, n=n, eta=eta, stretch=stretch)
        efield2, bfield2 = cluster2p(x, cluster.children[2], particles, parindices; p_kernel=p_kernel, n=n, eta=eta, stretch=stretch)
        efield = efield1 + efield2
        bfield = bfield1 + bfield2
    end
    return efield, bfield
end

function admissible(x::SVector{3,T}, bbox::BBox{T}; stretch::SVector{3,T}=SVector(1.0,1.0,1.0), eta::T) where {T}
    bmin = bbox.bmin
    bmax = bbox.bmax
    xc = (bmin + bmax)/2.0
    r = norm(stretch .* (bmax - xc))
    R = norm(stretch .* (x - xc))
    return r/R < eta
end

function cluster2p_brutal(x::SVector{3,T}, cluster::Cluster{T}, particles::Particles, parindices::Vector{Int64}) where {T}
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    pindex_lo = cluster.pindex_lo
    pindex_hi = cluster.pindex_hi
    pos = particles.positions
    mom = particles.momenta
    for k in pindex_lo:pindex_hi
        xj = pos[parindices[k]]
        pj = mom[parindices[k]]
        K = kernel_relativity(x, xj, pj)
        efield += sqrt(1.0 + dot(pj, pj))*K
        bfield += cross(pj, K)
    end
    return efield, bfield
end
