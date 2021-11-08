function kernel_relativity(xi::SVector, xj::SVector{3,T}, pj::SVector{3,T}) where {T}
    R = xi - xj
    return R / sqrt(dot(R, R) + dot(pj, R)^2 + eps())^3
end

function kernel_coulomb(xi::SVector, xj::SVector{3,T}) where {T}
    R = xi - xj
    return R / sqrt(dot(R, R) + eps())^3
end


