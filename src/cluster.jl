struct Cluster{T}
    level::Int
    parindex_lo::Int
    parindex_hi::Int
    bbox::BBox{T}
    macroparticles::MacroParticles{T}
    children::Union{Nothing,Tuple{Cluster{T},Cluster{T}}}
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

function subdivide(particles::Particles{T}, parindices::Vector{Int}, lo, hi, level; n, N0, stretch=SVector(1.0,1.0,1.0)) where {T}
    npar = (hi - lo + 1)
    pos = particles.positions
    bbox = BBox(particles, parindices, lo, hi)
    macroparticles = MacroParticles(bbox, n)
    if npar <= N0
        return Cluster(level, lo, hi, bbox, macroparticles, nothing)
    else
        splitdir = argmax( stretch .* (bbox.bmax - bbox.bmin) )
        k = split_median!(parindices, pos, splitdir, lo, hi)
        children = (subdivide(particles, parindices, lo, k, level + 1; n=n, N0=N0, stretch=stretch), subdivide(particles, parindices, k + 1, hi, level + 1; n=n, N0=N0, stretch=stretch))
        return Cluster(level, lo, hi, bbox, macroparticles, children)
    end
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

function partition!(parindices::Vector{Int}, pos::AbstractVector{SVector{3,T}}, dim, lo, hi, parindex) where {T}
pval = pos[parindices[parindex]][dim]
    @swap!(parindices[parindex], parindices[hi])
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
        parindex = rand(lo:hi) # floor(Int, (lo + hi) / 2.0)
        parindex = partition!(parindices, pos, dim, lo, hi, parindex)
        if parindex == k
            return k
        elseif parindex > k
            hi = parindex - 1
        else
            lo = parindex + 1
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
