"""
    updateParticlesField!(particles::Particles, alg::BruteForce; lambda)

Evaluate and update the space charge field experienced by each particle in `beam` using
`alg` method.
"""
function updateParticlesField!(particles::Particles{T}, alg::BruteForce; lambda) where {T}
    q = particles.charge
    npar = particles.npar
    @inbounds for i in 1:npar
        particles.self_efields[i] = SVector(0.0, 0.0, 0.0)
        particles.self_bfields[i] = SVector(0.0, 0.0, 0.0)
        xi = particles.positions[i]
        amp = 2.8179403699772166e-15 * q / lambda
        for j in 1:npar
            xj = particles.positions[j]
            pj = particles.momenta[j]
            K = kernel_relativity(xi, xj, pj)
            particles.self_efields[i] += amp * sqrt(1.0 + dot(pj, pj)) * K
            particles.self_bfields[i] += amp * cross(pj, K)
        end
    end
end

function updateParticlesField!(particles::Particles{T}, alg::TreecodeStretch{T}; lambda) where {T}
    (;n, N0, eta) = alg

    q = particles.charge
    p_avg = sum(particles.momenta) / particles.npar
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    stretch = SVector(1.0,1.0,g_avg)

    ct = ClusterTree(particles; n=n, N0=N0, stretch=stretch)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices; p_kernel=p_avg, eta=eta, stretch=stretch)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end
end

function updateParticlesField!(particles::Particles{T}, alg::TreecodeUniform{T}; lambda) where {T}
    (;n, N0, eta) = alg

    q = particles.charge
    p_avg = sum(particles.momenta) / particles.npar

    ct = ClusterTree(particles; n=n, N0=N0)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices; p_kernel=p_avg, eta=eta)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end
end

function updateParticlesField!(particles::Particles{T}, alg::TreecodeAvgRestFrame{T}; lambda) where {T}
    (;n, N0, eta) = alg

    q = particles.charge
    p_avg = sum(particles.momenta) / particles.npar
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))

    # transform particles to rest-frame
    boostParticlesPosition!(particles, p_avg, g_avg)
    transformParticlesMomentum!(particles, p_avg)

    # evaluate particles field in the rest-frame
    ct = ClusterTree(particles; n=n, N0=N0)
    amp = 2.8179403699772166e-15 * q / lambda
    @inbounds for i in 1:particles.npar
        x = particles.positions[i]
        (efield, bfield) = cluster2p(x, ct.root, particles, ct.parindices; p_kernel=SVector(0.0,0.0,0.0), eta=eta)
        particles.self_efields[i] = amp * efield
        particles.self_bfields[i] = amp * bfield
    end

    # transform particle distribution and field back to the lab-frame
    # deboost particles position
    boostParticlesPosition!(particles, -p_avg, 1.0/g_avg)
    transformParticlesMomentum!(particles, -p_avg)
    transformParticlesField!(particles, -p_avg)
end
