"""
    transformParticleDistribution!(beam::Particles, fp::SVector{3,T}) where {T}

Transform the position and the momentum of each particle contained in `beam`
to the inertial frame moving with the momentum `fp`.
"""
function transformParticlesDistribution!(beam::Particles, fp::SVector{3,T}) where {T}
    fg = sqrt(1.0 + dot(fp, fp))
    pos = beam.positions
    mom = beam.momentums
    for k in 1:beam.npar
        x = pos[k]
        p = mom[k]
        g = sqrt(1.0 + dot(p, p))
        pos[k] = transformPosition(x, fg, fp)
        mom[k] = transformMomentum(p, fg, fp)
    end
end

function transformPosition(x::SVector{3,T}, fg::T, fp::SVector{3,T}) where {T}
    return x + 1.0/(fg + 1.0) * dot(fp, x)*fp
end

function transformMomentum(p::SVector{3,T}, fg::T, fp::SVector{3,T}) where {T}
    g = sqrt(1.0 + dot(p,p))
    return  p + 1.0/(fg + 1.0) * dot(fp, p)*fp - g*fp
end


"""
    transformParticlesField!(beam::Particles, fp::SVector{3,T}) where {T}

Transform the EM field experienced by each particle contained in `beam`
to the inertial frame moving with the momentum `fp`.
"""
function transformParticlesField!(beam::Particles, fp::SVector{3, T}) where {T}
    fg = sqrt(1.0 + dot(fp, fp))
    efields = beam.self_efields
    bfields = beam.self_bfields
    for i in 1:beam.npar
      efield_f, bfield_f = transformEMField(efields[i], bfields[i], fg, fp)
      efields[i] = efield_f
      bfields[i] = bfield_f
    end
  end

function transformEMField(efield::SVector{3, T}, bfield::SVector{3, T}, fg::T, fp::SVector{3, T}) where {T}
  efield_f = fg*efield + cross(fp, bfield) - 1.0/(fg + 1.0)*dot(fp, efield)*fp
  bfield_f = fg*bfield - cross(fp, efield) - 1.0/(fg + 1.0)*dot(fp, bfield)*fp
  return efield_f, bfield_f
end
