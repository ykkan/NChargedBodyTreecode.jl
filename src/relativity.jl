function particle2frame!(beam::Particles, p_avg::SVector{3,T}) where {T}
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    pos = beam.positions
    mom = beam.momentums
    for k in 1:beam.npar
        x = pos[k]
        p = mom[k]
        g = sqrt(1.0 + dot(p, p))
        pos[k] = x + 1.0/(g_avg + 1.0) * dot(p_avg, x)*p_avg
        mom[k] = p + 1.0/(g_avg + 1.0) * dot(p_avg, p)*p_avg - g*p_avg
    end
end

function lab2avgframe!(beam::Particles, p_avg::SVector{3,T}) where {T}
    g_avg = sqrt(1.0 + dot(p_avg, p_avg))
    pos = beam.positions
    mom = beam.momentums
    for k in 1:beam.npar
        x = pos[k]
        p = mom[k]
        g = sqrt(1.0 + dot(p, p))
        pos[k] = x + 1.0/(g_avg + 1.0) * dot(p_avg, x)*p_avg
        mom[k] = p + 1.0/(g_avg + 1.0) * dot(p_avg, p)*p_avg - g*p_avg
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
        mom[k] = p + 1.0/(g_avg + 1.0)*dot(p_avg, p)*p_avg + g*p_avg
    end
end

function field2frame(efield::SVector{3, T}, bfield::SVector{3, T}, g_f::T, p_f::SVector{3, T}) where {T}
  efield_f = g_f*efield + cross(p_f, bfield) - 1.0/(g_f + 1.0)*dot(p_f, efield)*p_f
  bfield_f = g_f*bfield - cross(p_f, efield) - 1.0/(g_f + 1.0)*dot(p_f, bfield)*p_f
  return efield_f, bfield_f
end


function fields2frame!(beam::Particles, p_f::SVector{3, T}) where {T}
  g_f = sqrt(1.0 + dot(p_f, p_f))
  efields = beam.self_efields
  bfields = beam.self_bfields
  for i in 1:beam.npar
    efield, bfield = field2frame(efields[i], bfields[i], g_f, p_f)
    efields[i] = efield
    bfields[i] = bfield
  end
end

function field2labframe(e_p::SVector{3, T} , b_p::SVector{3, T}, g_avg::T, p_avg::SVector{3, T}) where {T}
    efield = g_avg*e_p - cross(p_avg, b_p) - 1.0/(g_avg + 1.0)*dot(p_avg, e_p)*p_avg
    bfield = g_avg*b_p + cross(p_avg, e_p) - 1.0/(g_avg + 1.0)*dot(p_avg, b_p)*p_avg
    return efield, bfield
end

function field2avgframe(e_p::SVector{3, T} , b_p::SVector{3, T}, g_avg::T, p_avg::SVector{3, T}) where {T}
    efield = g_avg*e_p + cross(p_avg, b_p) - 1.0/(g_avg + 1.0)*dot(p_avg, e_p)*p_avg
    bfield = g_avg*b_p - cross(p_avg, e_p) - 1.0/(g_avg + 1.0)*dot(p_avg, b_p)*p_avg
    return efield, bfield
end
