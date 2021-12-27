function update_particle_field!(particles::Particles; n, eta, threshold, lambda) where {T}
    q = particles.charge
    p_avg = sum(particles.momentums) / particles.npar

    ct = ClusterTree(particles; n=n, threshold=threshold)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices, p_avg; n=n, eta=eta)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end
end

function update_particle_field_stretch!(particles::Particles; n, eta, threshold, lambda) where {T}
    q = particles.charge
    p_avg = sum(particles.momentums) / particles.npar
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    stretch = SVector(1.0,1.0,g_avg)
    ct = ClusterTree(particles; n=n, threshold=threshold, stretch=stretch)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices, p_avg; n=n, eta=eta, stretch=stretch)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end
end

function update_particle_field_AVGRF!(particles::Particles; n, eta, threshold, lambda) where {T}
    q = particles.charge
    p_avg = sum(particles.momentums) / particles.npar
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))

    # transform particles to rest-frame
    transformParticleDistribution!(particles, p_avg)

    # evaluate particles field in the rest-frame
    ct = ClusterTree(particles; n=n, threshold=threshold)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices, SVector(0.0,0.0,0.0); n=n, eta=eta)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end

    # transform particle field back to the lab-frame
    transformParticlesField!(particles, -p_avg)
end

function update_particle_field_brutal!(particles::Particles; lambda)
    q = particles.charge
    npar = particles.npar
    @inbounds for i in 1:npar
        particles.self_efields[i] = SVector(0.0, 0.0, 0.0)
        particles.self_bfields[i] = SVector(0.0, 0.0, 0.0)
        xi = particles.positions[i]
        amp = 2.8179403699772166e-15 * q / lambda
        for j in 1:npar
            xj = particles.positions[j]
            pj = particles.momentums[j]
            K = kernel_relativity(xi, xj, pj)
            particles.self_efields[i] += amp * sqrt(1.0 + dot(pj, pj)) * K
            particles.self_bfields[i] += amp * cross(pj, K)
        end
    end
end
