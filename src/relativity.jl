"""
    boostParticlesPosition!(beam::Particles, fp::SVector{3,T}, factor::T) where {T}

Boost the position of each particle contained in `beam` to the
inertial frame moving with the momentum `fp` by `factor`.
"""
function boostParticlesPosition!(beam::Particles, fp::SVector{3,T}, factor::T) where {T}
    pos = beam.positions
    @inbounds for k in 1:beam.npar
        x = pos[k]
        pos[k] = boostPosition(x, fp, factor)
    end
end

function boostPosition(x::SVector{3, T}, fp::SVector{3,T}, factor::T) where {T}
    fp_hat = fp / norm(fp)
    return x + (factor - 1.0) * dot(x, fp_hat) * fp_hat
end


"""
    transformParticlesMomentum!(beam::Particles, fp::SVector{3,T}) where {T}

Transform the momentum of each particle contained in `beam`
to the inertial frame moving with the momentum `fp`.
"""
function transformParticlesMomentum!(beam::Particles, fp::SVector{3,T}) where {T}
    fg = sqrt(1.0 + dot(fp, fp))
    mom = beam.momenta
    @inbounds for k in 1:beam.npar
       p = mom[k]
       mom[k] = transformMomentum(p, fg, fp)
    end
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
    @inbounds for i in 1:beam.npar
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
