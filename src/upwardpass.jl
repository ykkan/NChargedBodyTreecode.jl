# values of Lagragian polynomials on chevyschev point evaluated at x
function weight(x::T, nodes::Vector{T}) where {T}
    n = length(nodes) - 1
    w = [(-1.0)^(j) for j in 0:n]
    w[1] = 0.5
    w[n + 1] = 0.5 * w[n + 1]
    diffs = x .- nodes
    flag = findfirst(x->abs(x) < eps(T), diffs)
    values = zeros(n+1)
    if flag == nothing
        values .= w ./ diffs
        sum_v = sum(values)
        return values ./ sum_v
    else
        values[flag] = 1.0
        return values
    end
end

"""
    deposit!(mp::MacroParticles{T}, pos::SVector{3,T}, mom::SVector{3,T})
deposite the momentum and gamma of one particle (with position `pos` and momentum `mom`)
to a group of macro-particles `mp`
"""
function deposit!(mp::MacroParticles{T}, pos::SVector{3,T}, mom::SVector{3,T}) where {T}
    x, y, z = pos
    weights_x = weight(x, mp.xcoords)
    weights_y = weight(y, mp.ycoords)
    weights_z = weight(z, mp.zcoords)

    n = length(mp.xcoords) - 1
    gamma = sqrt(1.0 + dot(mom, mom))
    for k in 1:(n + 1)
        for j in 1:(n + 1)
            for i in 1:(n + 1)
                weight = weights_x[i] * weights_y[j] * weights_z[k]
                mp.gammas[i,j,k] += gamma * weight
                mp.momenta[i,j,k] += mom * weight
            end
        end
    end
end

function P2M!(cluster::Cluster{T}, particles::Particles{T}, parindices::Vector{Int}) where {T}
    pos = particles.positions
    mom = particles.momenta
    lo = cluster.parindex_lo
    hi = cluster.parindex_hi
    for i = lo:hi
        pindex = parindices[i]
        deposit!(cluster.macroparticles, pos[pindex], mom[pindex])
    end
end

function M2M!(parent::Cluster{T}, child::Cluster{T}) where {T}
    n = length(parent.macroparticles.xcoords) - 1
    child_xcoords = child.macroparticles.xcoords
    child_ycoords = child.macroparticles.ycoords
    child_zcoords = child.macroparticles.zcoords
    child_momenta = child.macroparticles.momenta
    for k in 1:(n+1)
        for j in 1:(n+1)
            for i in 1:(n+1)
                x = child_xcoords[i]
                y = child_ycoords[j]
                z = child_zcoords[k]
                p = child_momenta[i,j,k]
                deposit!(parent.macroparticles, SVector(x,y,z), p)
            end
        end
    end
end

function upwardpass!(cluster::Cluster{T}, particles::Particles{T}, parindices::Vector{Int}) where {T}
    if cluster.children == nothing
        P2M!(cluster, particles, parindices)
    else
        upwardpass!(cluster.children[1], particles, parindices)
        upwardpass!(cluster.children[2], particles, parindices)
        M2M!(cluster, cluster.children[1])
        M2M!(cluster, cluster.children[2])
    end
end

function upwardpass!(ct::ClusterTree)
    upwardpass!(ct.root, ct.particles, ct.parindices)
end
