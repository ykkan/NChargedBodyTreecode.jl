# This example computes the space-charge field of a particle beam with single momentum.
# The particle beam moves along z-direction with γ=30.

using NChargedBodyTreecode

# number of particles
const N = 10000

# create position and momentum distribution of N particles
positions = rand(3, N)
momenta = zeros(3, N)
momenta[3, :] .= 29.9833287

# create a particle beam in which each particle's charge and mass are -1 and 1
beam = Particles(; pos=positions, mom=momenta, charge=-1.0, mass=1.0)

# create 4 copies of a particle beam
beam0 = deepcopy(beam)
beam1 = deepcopy(beam)
beam2 = deepcopy(beam)
beam3 = deepcopy(beam)

# update particle field by brute-force method (lambda is a characteristic length for the normalization of length quantity)
updateParticlesField!(beam0, BruteForce(); lambda=1.0)

# update particle field by different treecode methods with the following parameters:
#   n: degree of interpolation
#   N0: maximum number of particles in the leaf cluster
#   eta: admissibility parameter
#   lambda: a characteristic length for the normalization of length quantity
####
const n = 4
const N0 = 128
const eta = 0.5
updateParticlesField!(beam1, TreecodeStretch(eta=eta, N0=N0, n=n); lambda=1.0)
updateParticlesField!(beam2, TreecodeAvgRestFrame(eta=eta, N0=N0, n=n); lambda=1.0)
updateParticlesField!(beam3, TreecodeUniform(eta=eta, N0=N0, n=n); lambda=1.0)

# relative errors of different treecodes (comparing to brute-force)
error1 = relerror(beam1, beam0)
error2 = relerror(beam2, beam0)
error3 = relerror(beam3, beam0)
println("----------------")
println("relative errors")
println("----------------")
println("Treecode-Stretch: $(error1)")
println("Treecode-AVGRF: $(error2)")
println("Treecode-Uniform: $(error3)")
