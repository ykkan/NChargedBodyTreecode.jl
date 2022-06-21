struct MacroParticles{T}
    xcoords::Vector{T}
    ycoords::Vector{T}
    zcoords::Vector{T}
    gammas::Array{T,3}
    momenta::Array{SVector{3,T},3}
    efields::Array{SVector{3,T},3}
    bfields::Array{SVector{3,T},3}
end

function MacroParticles(bbox::BBox, n)
    xmin, ymin, zmin = bbox.bmin
    xmax, ymax, zmax = bbox.bmax
    xcoords = cheb2(n, xmin, xmax)
    ycoords = cheb2(n, ymin, ymax)
    zcoords = cheb2(n, zmin, zmax)

    gammas = zeros(n + 1, n + 1, n + 1)
    momenta = fill(SVector(0.0, 0.0, 0.0), (n + 1, n + 1, n + 1))
    efields = fill(SVector(0.0, 0.0, 0.0), (n + 1, n + 1, n + 1))
    bfields = fill(SVector(0.0, 0.0, 0.0), (n + 1, n + 1, n + 1))
    return MacroParticles(xcoords, ycoords, zcoords, gammas, momenta, efields, bfields)
end
