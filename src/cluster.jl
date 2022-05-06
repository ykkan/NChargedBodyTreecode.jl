struct BBox{T}
    bmin::SVector{3,T}
    bmax::SVector{3,T}
end

function BBox(particles::Particles{T}, parindices::Vector{Int}, lo::Int, hi::Int) where {T}
    pos = particles.positions
    tmin = typemin(T)
    tmax = typemax(T)
    xmin, xmax = tmax, tmin
    ymin, ymax = tmax, tmin
    zmin, zmax = tmax, tmin
    for i in lo:hi
        (x, y, z) = pos[parindices[i]]
        xmin = x < xmin ? x : xmin
        xmax = x > xmax ? x : xmax
        ymin = y < ymin ? y : ymin
        ymax = y > ymax ? y : ymax
        zmin = z < zmin ? z : zmin
        zmax = z > zmax ? z : zmax
    end
    bmin = SVector(xmin, ymin, zmin)
    bmax = SVector(xmax, ymax, zmax)
    return BBox(bmin, bmax)
end

struct Cluster{T}
    level::Int
    pindex_lo::Int
    pindex_hi::Int
    bbox::BBox{T}
    macroparticles::MacroParticles{T}
    children::Union{Nothing,Tuple{Cluster{T},Cluster{T}}}
end

function subdivide(particles::Particles{T}, parindices::Vector{Int}, lo, hi, level; n, N0, stretch=SVector(1.0,1.0,1.0)) where {T}
    npar = (hi - lo + 1)
    pos = particles.positions
    bbox = BBox(particles, parindices, lo, hi)
    macroparticles = MacroParticles(bbox, n)
    xcoords, ycoords, zcoords = cluster_coord(bbox; n=n)
    cgammas, cmomenta = cluster_weight(particles, parindices, lo, hi, bbox; n=n)
    macroparticles.xcoords .= xcoords
    macroparticles.ycoords .= ycoords
    macroparticles.zcoords .= zcoords
    macroparticles.gammas .= cgammas
    macroparticles.momenta .= cmomenta
    if npar <= N0
        return Cluster(level, npar, bbox, lo, macroparticles, nothing)
    else
        splitdir = argmax( stretch .* (bbox.bmax - bbox.bmin) )
        k = split_median!(parindices, pos, splitdir, lo, hi)
        children = (subdivide(particles, parindices, lo, k, level + 1; n=n, N0=N0, stretch=stretch), subdivide(particles, parindices, k + 1, hi, level + 1; n=n, N0=N0, stretch=stretch))
        return Cluster(level, lo, hi, bbox, macroparticles, children)
    end
end


function cluster_coord(bbox::BBox{T}; n) where {T}
    # transform chebpts to [0,1]
    chebpts_shifted = [ (cos(pi * i / n) + 1.0) / 2.0 for i in 0:n]
    xmin, ymin, zmin = bbox.bmin
    xmax, ymax, zmax = bbox.bmax
    xcoords = chebpts_shifted .* (xmax - xmin) .+ xmin
    ycoords = chebpts_shifted .* (ymax - ymin) .+ ymin
    zcoords = chebpts_shifted .* (zmax - zmin) .+ zmin
    return xcoords, ycoords, zcoords
end

function cluster_weight(particles::Particles, parindices::Vector{Int}, lo, hi, bbox::BBox{T}; n) where {T}
    pos = particles.positions
    mom = particles.momentums
    bmin = bbox.bmin
    bmax = bbox.bmax
    chebpts = [cos(pi * i / n) for i in 0:n]
    w = [(-1.0)^(j) for j in 0:n]
    w[1] = 0.5
    w[n + 1] = 0.5 * (-1.0)^n
    diffx = zeros(n + 1)
    diffy = zeros(n + 1)
    diffz = zeros(n + 1)
    ax = zeros(n + 1)
    ay = zeros(n + 1)
    az = zeros(n + 1)
    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0
    gamma_hat = zeros(n + 1, n + 1, n + 1)
    mom_hat  = fill(SVector(0.0, 0.0, 0.0), (n + 1, n + 1, n + 1))
    for i = lo:hi
        pindex = parindices[i]
        x, y, z = (pos[pindex] - bmin) ./ (bmax - bmin) * 2 .- 1.0
        diffx .= (x .- chebpts)
        diffy .= (y .- chebpts)
        diffz .= (z .- chebpts)
        indx = findfirst(x -> abs(x) < eps(0.0), diffx)
        indy = findfirst(x -> abs(x) < eps(0.0), diffy)
        indz = findfirst(x -> abs(x) < eps(0.0), diffz)
        if indx === nothing
            ax .= w ./ diffx
            sum_x = sum(ax)
        else
            ax .= 0.0
            ax[indx] = 1.0
            sum_x = 1.0
        end
        if indy === nothing
            ay .= w ./ diffy
            sum_y = sum(ay)
        else
            ay .= 0.0
            ay[indy] = 1.0
            sum_y = 1.0
        end
        if indz === nothing
            az .= w ./ diffz
            sum_z = sum(az)
        else
            az .= 0.0
            az[indz] = 1.0
            sum_z = 1.0
        end
        p = mom[pindex]
        gamma = sqrt(1.0 + dot(p, p))
        for k in 1:(n + 1)
            for j in 1:(n + 1)
                for i in 1:(n + 1)
                    weight = ax[i] * ay[j] * az[k] / (sum_x * sum_y * sum_z)
                    gamma_hat[i,j,k] += gamma * weight
                    mom_hat[i,j,k] += p * weight
                end
            end
        end
    end
    return gamma_hat, mom_hat
end


function print_bbox(cluster::Cluster{T}) where {T}
    bbox = cluter.bbox
    level = cluster.level
    xmin, ymin, zmin = bbox.bmin
    xmax, ymax, zmax = bbox.bmax
    npar = cluster.npar
    println("$(xmin) $(ymin) $(zmin) $(xmax) $(ymax) $(zmax) $(level) $(npar)")
    children = cluster.children
    if children !== nothing
        print_bbox(children[1])
        print_bbox(children[2])
    end
end

macro swap!(a, b)
    return quote
        $(esc(a)), $(esc(b)) = $(esc(b)), $(esc(a))
    end
end

function partition!(parindices::Vector{Int}, pos::AbstractVector{SVector{3,T}}, dim, lo, hi, pindex) where {T}
pval = pos[parindices[pindex]][dim]
    @swap!(parindices[pindex], parindices[hi])
    j = lo
    for i in lo:(hi - 1)
        if pos[parindices[i]][dim] < pval
            @swap!(parindices[j], parindices[i])
            j = j + 1
        end
    end
    @swap!(parindices[hi], parindices[j])
    return j
end

function split_median!(parindices::Vector{Int}, pos::AbstractVector{SVector{3,T}}, dim, lo, hi) where {T}
    k = floor(Int, (lo + hi) / 2.0)
    while true
        if lo == hi
            return lo
        end
        pindex = rand(lo:hi) # floor(Int, (lo + hi) / 2.0)
        pindex = partition!(parindices, pos, dim, lo, hi, pindex)
        if pindex == k
            return k
        elseif pindex > k
            hi = pindex - 1
        else
            lo = pindex + 1
        end
    end
end

function num_cluster(c::Cluster{T}) where {T}
    children = c.children
    if (children === nothing)
        return 1
    else
        return 1 + num_cluster(children[1]) + num_cluster(children[2])
    end
end
