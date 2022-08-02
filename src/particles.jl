using StaticArrays
struct Particles{T,V <: AbstractVector{SVector{3,T}}}
    npar::Int
    charge::T
    mass::T
    positions::V
    momenta::V
    efields::V
    bfields::V
end

# general interface for creating particle bunch
function Particles(;pos, mom, charge=-1.0, mass=1.0)
    return Particles(pos, mom; charge=charge, mass=mass)
end

function Particles(pos::Vector{SVector{3,T}}, mom::Vector{SVector{3,T}}; charge=-1.0, mass=1.0) where {T}
    @assert length(pos) == length(mom)
    npar = length(pos)
    return Particles(npar, charge, mass, pos, mom, fill(SVector(0.0, 0.0, 0.0), npar), fill(SVector(0.0, 0.0, 0.0), npar))
end

function Particles(pos::Matrix{T}, mom::Matrix{T}; charge=-1.0, mass=1.0) where {T}
    @assert size(pos, 1) == 3
    @assert size(mom, 1) == 3
    @assert size(pos, 2) == size(mom, 2)
    npar = size(pos, 2)
    pos_new = reinterpret(reshape, SVector{3,T}, pos)
    mom_new = reinterpret(reshape, SVector{3,T}, mom)
    return Particles(npar, charge, mass, pos_new, mom_new, reinterpret(reshape, SVector{3,T}, zeros(3, npar)), reinterpret(reshape, SVector{3,T}, zeros(3, npar)))
end

# type alias for Particles
const Beam{T} = Particles{T}

function output(beam::Particles, fname)
    pos = beam.positions
    mom = beam.momenta
    open(fname, "w") do io
        writedlm(io, pos)
    end
end
