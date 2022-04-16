# NChargedBodyTreecode.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6464509.svg)](https://doi.org/10.5281/zenodo.6464509)

by Yi-Kai Kan (<yikai.kan@desy.de>)

This package provides treecode algorithms for calculating the relativistic space-charge field. The treeocodes are implemented based on the tensor-product interpolation with Lagragian polynomials.

# Installation
The package can be installed using Julia's REPL
```julia
julia> using pkg
julia> pkg.add("https://github.com/ykkan/NChargedBodyTreecode.git.git")
```
or with Pkg mode (hitting `]` in the command prompt)
```julia
pkg> add("https://github.com/ykkan/NChargedBodyTreecode.git") 
```

# Usage
## Creating a Charged Particle Beam
A charged particle beam can be created by providing an array of particle position and an array of particle momentum.
``` julia
using NChargedBodyTreecode

# number of particles
const N = 1000    

# create position and momentum distribution of N particles
positions = rand(3, N)
momenta = zeros(3, N)

# create a particle beam in which each particle's charge and mass are -1 and 1 
beam = Particles(; pos=positions, mom=momenta, charge=-1.0, mass=1.0) 
```

## Space-Charge Field Calculation
The space-charge field of a particle beam can be evaluated using `updateParticlesField!` with different algorithms:
* `BruteForce`: brute-forece method
* `TreecodeStretch`: treecode with stretched admissibility condition
* `TreecodeAvgRestFrame`: treecode using average rest-frame technique
* `TreecodeUniform`: treecode with conventional admissibility condition 
``` julia
using NChargedBodyTreecode
# update particle field by brute-force method (lambda is a characteristic length for the normalization of length quantity)

updateParticlesField!(beam0, BruteForce(); lambda=1.0)

# update particle field by different treecode methods with the following parameters: 
#   n: degree of interpolation
#   N0: maximum number of particles in the leaf cluster
#   eta: admissibility parameter 
#   lambda: a characteristic length for the normalization of length quantity
####
const n = 2 
const N0 = 20  
const eta = 0.5
updateParticlesField!(beam1, TreecodeStretch(eta=eta, N0=N0, n=n); lambda=1.0)
updateParticlesField!(beam2, TreecodeAvgRestFrame(eta=eta, N0=N0, n=n); lambda=1.0)
updateParticlesField!(beam3, TreecodeUniform(eta=eta, N0=N0, n=n); lambda=1.0)

# relative errors of different treecodes (comparing to brute-force)
error1 = relerror(beam1, beam0) 
error2 = relerror(beam2, beam0)
error3 = relerror(beam3, beam0)
```
The space-charge field (both E- and B-field) of `beam` is stored in `beam.self_efields` and `beam.self_bfields` (both have a type `Vector{SVector{3,T}}`, see [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) for more details).

After the update of space-charge field of `beam`, the space-charge field on i-th particle can be accessed by, for example
```julia
efield = beam.self_efields[i]
bfield = beam.self_bfields[i]
```

 More details can found in the demo file from `examples/`.
