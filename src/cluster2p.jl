function cluster2p(x::SVector{3,T}, cluster::Cluster{T}, particles::Particles, parindices::Vector{Int64}; p_kernel, eta, stretch=SVector(1.0,1.0,1.0)) where {T}
    if cluster.children === nothing
        P2P(x, cluster, particles, parindices)
    elseif admissible(x, cluster; eta=eta, stretch=stretch)
        M2P(x, cluster; p_kernel=p_kernel)
    else
        efield1, bfield1 = cluster2p(x, cluster.children[1], particles, parindices; p_kernel=p_kernel, eta=eta, stretch=stretch)
        efield2, bfield2 = cluster2p(x, cluster.children[2], particles, parindices; p_kernel=p_kernel, eta=eta, stretch=stretch)
        efield = efield1 + efield2
        bfield = bfield1 + bfield2
        return efield, bfield
    end
end

function admissible(x::SVector{3,T}, cluster::Cluster{T}; stretch::SVector{3,T}=SVector(1.0,1.0,1.0), eta::T) where {T}
    bmin = cluster.bbox.bmin
    bmax = cluster.bbox.bmax
    xc = (bmin + bmax)/2.0
    r = norm(stretch .* (bmax - xc))
    R = norm(stretch .* (x - xc))
    return r/R < eta
end

function M2P(x::SVector{3,T}, cluster::Cluster{T}; p_kernel) where {T}
    macroparticles = cluster.macroparticles
    xcoords = macroparticles.xcoords
    ycoords = macroparticles.ycoords
    zcoords = macroparticles.zcoords
    cgamma = macroparticles.gammas
    cmomenta = macroparticles.momenta
    n = length(xcoords) - 1
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    for k in 1:(n+1)
        sz = zcoords[k]
        for j in 1:(n+1)
            sy = ycoords[j]
            for i in 1:(n+1)
                sx = xcoords[i]
                K = kernel_relativity(x, SVector(sx,sy,sz), p_kernel)
                efield += cgamma[i,j,k]*K
                bfield += cross(cmomenta[i,j,k], K)
            end
        end
    end
    return efield, bfield
end

function P2P(x::SVector{3,T}, cluster::Cluster{T}, particles::Particles, parindices::Vector{Int64}) where {T}
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    parindex_lo = cluster.parindex_lo
    parindex_hi = cluster.parindex_hi
    pos = particles.positions
    mom = particles.momenta
    for k in parindex_lo:parindex_hi
        xj = pos[parindices[k]]
        pj = mom[parindices[k]]
        K = kernel_relativity(x, xj, pj)
        efield += sqrt(1.0 + dot(pj, pj))*K
        bfield += cross(pj, K)
    end
    return efield, bfield
end
