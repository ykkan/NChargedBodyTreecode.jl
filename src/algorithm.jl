import Base.@kwdef

struct BruteForce end

@kwdef struct TreecodeStretch{T}
    n::Int
    N0::Int
    eta::T
end

@kwdef struct TreecodeAvgRestFrame{T}
    n::Int
    N0::Int
    eta::T
end

@kwdef struct TreecodeUnstretch{T}
    n::Int
    N0::Int
    eta::T
end
