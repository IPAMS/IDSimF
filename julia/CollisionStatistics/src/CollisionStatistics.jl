module CollisionStatistics
#export particle

#using Distributed
#using LinearAlgebra
#using Random
#using Distributions
#using Plots
#using StatsBase, StatsPlots, KernelDensity
#using Profile
#using ProfileSVG
#using SharedArrays
#using DelimitedFiles


const π2 = π*2
const sigma_maxwellboltzmann = 0.627

include("core.jl")
include("simulation.jl")

end # module
