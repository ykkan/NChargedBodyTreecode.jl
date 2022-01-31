
"""
    relerror(beam1::Particles, beam2::Particles)

Compute the relative error of particle field from `beam1` comparing to `beam2`.
"""
function relerror(beam1::Particles, beam2::Particles)
    beam1.npar == beam2.npar || error("npar mismatches")
    npar = beam1.npar
    efield_error = relerror(beam1.self_efields, beam2.self_efields, npar)
    bfield_error = relerror(beam1.self_bfields, beam2.self_bfields, npar)
    return max(efield_error, bfield_error)
end

"""
    relerror(A::AbstractVector{SVector{3,T}}, B::AbstractVector{SVector{3,T}}, n) where {T}

Compute the relative error of an array of Vector3 `A` comparing to another array of Vector3
`B`.
"""
function relerror(A::AbstractVector{SVector{3,T}}, B::AbstractVector{SVector{3,T}}, n) where {T}
    denominator = 0.0
    numerator = 0.0
    for i in 1:n
        a = A[i]
        b = B[i]
        numerator += norm(a - b)^2
        denominator += norm(b)^2
    end
    return sqrt(numerator / denominator)
end
