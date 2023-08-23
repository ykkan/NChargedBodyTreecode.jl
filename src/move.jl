function move!(particles::Particles{T}, dt::T, alg::A; lambda::T) where {A, T}
    q = particles.charge
    m = particles.mass
    npar = particles.npar

    @inbounds for i in 1:npar
        p = particles.momenta[i]
        particles.positions[i] += 0.5 * dt * p / sqrt(1.0 + dot(p, p))
    end

    updateParticlesField!(particles, alg; lambda=lambda)

    @inbounds for i in 1:npar
        efield = particles.efields[i]
        bfield = particles.bfields[i]

        p = particles.momenta[i]
        p_minus = p + 0.5 * dt * (q / m) * efield
        gamma_minus = sqrt(1.0 + dot(p_minus, p_minus))
        t_vec = 0.5 * dt * (q / m) * bfield / gamma_minus
        s_vec = 2.0 * t_vec / ( 1.0 + dot(t_vec, t_vec) )
        p_plus = p_minus + cross(p_minus + (cross(p_minus, t_vec)), s_vec)
        particles.momenta[i] = p_plus + 0.5 * dt * (q / m) * efield
    end

    println("$(particles.momenta[1])")

    @inbounds for i in 1:npar
        p = particles.momenta[i]
        particles.positions[i] += 0.5 * dt * p / sqrt(1.0 + dot(p, p))
    end
end
