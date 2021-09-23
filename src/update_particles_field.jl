function kernel(xi::SVector, xj::SVector{3,T}, mom_avg::SVector{3,T}) where {T}
    xij = xi - xj
    return xij / sqrt(dot(xij, xij) + dot(mom_avg, xij)^2 + eps())^3
end

function admissible(x::SVector{3,T}, bbox::BBox{T}; ita::T) where {T}
    bmin = bbox.bmin
    bmax = bbox.bmax
    xc = (bmin + bmax)/2.0
    r = norm(bmax - xc)
    R = norm(x - xc)
    return r/R < ita
end

function cluster_field_p2p(x::SVector{3,T}, particles::Particles, parindices::Vector{Int64}, index_start, index_end, mom_avg::SVector{3,T}) where {T}
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    pos = particles.positions
    mom = particles.momentums
    for k in index_start:index_end
        xj = pos[parindices[k]]
        pj = mom[parindices[k]]
        K = kernel(x, xj, mom_avg)
        efield += sqrt(1.0 + dot(pj, pj))*K
        bfield += cross(pj, K)
    end
    return efield, bfield
end

function cluster_field(x::SVector{3,T}, cluster::Cluster{T}, particles::Particles, parindices::Vector{Int64}, mom_avg::SVector{3,T}; n, ita) where {T}
    efield = SVector(0.0,0.0,0.0)
    bfield = SVector(0.0,0.0,0.0)
    if cluster.children === nothing
        index_start = cluster.pindex_start
        index_end = index_start + cluster.npar - 1
        efield, bfield = cluster_field_p2p(x, particles, parindices, index_start, index_end, mom_avg)
    elseif admissible(x, cluster.bbox; ita=ita)
        xcoords = cluster.xcoords
        ycoords = cluster.ycoords
        zcoords = cluster.zcoords
        gamma_hat = cluster.gamma_hat
        mom_hat = cluster.mom_hat
        for k in 1:(n+1)
            for j in 1:(n+1)
                for i in 1:(n+1)
                    sx = xcoords[i]
                    sy = ycoords[j]
                    sz = zcoords[k]
                    K = kernel(x, SVector(sx,sy,sz),mom_avg)
                    efield += gamma_hat[i,j,k]*K
                    bfield += cross(mom_hat[i,j,k],K)
                end
            end
        end
    else
        efield1, bfield1 = cluster_field(x, cluster.children[1], particles, parindices, mom_avg; n=n, ita=ita)
        efield2, bfield2 = cluster_field(x, cluster.children[2], particles, parindices, mom_avg; n=n, ita=ita)
        efield = efield1 + efield2
        bfield = bfield1 + bfield2
    end
    return efield, bfield
end


function update_particle_field!(particles::Particles; n, ita, threshold, lambda) where {T}
    q = particles.charge
    mom_avg = sum(particles.momentums) / particles.npar
    ct = ClusterTree(particles; threshold=threshold)
    for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster_field(x, ct.root, particles, ct.parindices, mom_avg; n=n, ita=ita)
        particles.self_efields[i] = 2.8179403699772166e-15 * q / lambda * efield
        particles.self_bfields[i] = 2.8179403699772166e-15 * q / lambda * bfield
    end
end

function update_particle_field_brutal!(beam::Particles; lambda)
    q = beam.charge
    @inbounds for i in 1:beam.npar
        beam.self_efields[i] = SVector(0.0, 0.0, 0.0)
        beam.self_bfields[i] = SVector(0.0, 0.0, 0.0)
        xi = beam.positions[i]
        for j in 1:beam.npar
            xj = beam.positions[j]
            pj = beam.momentums[j]
            xij = xi - xj
            kernel = xij / sqrt(dot(xij, xij) + dot(pj, xij)^2 + eps())^3
            beam.self_efields[i] += 2.8179403699772166e-15 * q / lambda * sqrt(1.0 + dot(pj, pj)) * kernel
            beam.self_bfields[i] += 2.8179403699772166e-15 * q / lambda * cross(pj, kernel)
        end
    end
end
