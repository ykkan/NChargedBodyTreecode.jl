function lab2avgframe!(beam::Particles, p_avg::SVector{3,T}) where {T}
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    p_avg_u = p_avg / norm(p_avg)
    pos = beam.positions
    mom = beam.momentums
    for k in 1:beam.npar
        x = pos[k]
        p = mom[k]
        g = sqrt(1.0 + dot(p, p))
        pos[k] = x + (g_avg - 1.0)*dot(x, p_avg_u)*p_avg_u 
        mom[k] = p + (g_avg - 1.0)*dot(p, p_avg_u)*p_avg_u - g*p_avg
    end
end

function avg2labframe!(beam::Particles, p_avg::SVector{3, T}) where {T}
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    p_avg_u = p_avg / norm(p_avg)
    pos = beam.positions
    mom = beam.momentums
    for k in 1:beam.npar
        x = pos[k]
        p = mom[k]
        g = sqrt(1.0 + dot(p, p))
        pos[k] = x + (1.0/g_avg - 1.0)*dot(x, p_avg_u)*p_avg_u 
        mom[k] = p + (g_avg - 1.0)*dot(p, p_avg_u)*p_avg_u + g*p_avg
    end
end

function field2labframe(e_p::SVector{3, T} , b_p::SVector{3, T}, g_avg::T, p_avg::SVector{3, T}) where {T}
    efield = g_avg*e_p - cross(p_avg, b_p) - 1.0/(g_avg + 1.0)*dot(p_avg, e_p)*p_avg
    bfield = g_avg*b_p + cross(p_avg, e_p) - 1.0/(g_avg + 1.0)*dot(p_avg, b_p)*p_avg
    return efield, bfield
end
