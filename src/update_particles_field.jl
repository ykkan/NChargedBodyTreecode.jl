function update_particle_field!(particles::Particles; n, ita, threshold, lambda) where {T}
    q = particles.charge
    ct = ClusterTree(particles; n=n, threshold=threshold)
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices; n=n, ita=ita)
        particles.self_efields[i] = (2.8179403699772166e-15 * q / lambda) * efield
        particles.self_bfields[i] = (2.8179403699772166e-15 * q / lambda) * bfield
    end
end

function update_particle_field_booster!(particles::Particles; n, ita, threshold, lambda) where {T}
    q = particles.charge
    mom_avg = sum(particles.momentums) / particles.npar
    g_avg = sqrt(1.0 + dot(mom_avg, mom_avg))
    lab2avgframe!(particles, mom_avg)

    ct = ClusterTree(particles; n=n, threshold=threshold)
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices; n=n, ita=ita)
        (efield, bfield) = field2labframe(efield, bfield, g_avg, mom_avg) 
        particles.self_efields[i] = (2.8179403699772166e-15 * q / lambda) * efield
        particles.self_bfields[i] = (2.8179403699772166e-15 * q / lambda) * bfield
    end

    avg2labframe!(particles, mom_avg)
end

function update_particle_field_brutal!(particles::Particles; lambda)
    q = particles.charge
    npar = particles.npar
    @inbounds for i in 1:npar
        particles.self_efields[i] = SVector(0.0, 0.0, 0.0)
        particles.self_bfields[i] = SVector(0.0, 0.0, 0.0)
        xi = particles.positions[i]
        for j in 1:particles.npar
            xj = particles.positions[j]
            pj = particles.momentums[j]
            K = kernel_relativity(xi, xj, pj)
            particles.self_efields[i] += 2.8179403699772166e-15 * q / lambda * sqrt(1.0 + dot(pj, pj)) * K
            particles.self_bfields[i] += 2.8179403699772166e-15 * q / lambda * cross(pj, K)
        end
    end
end
