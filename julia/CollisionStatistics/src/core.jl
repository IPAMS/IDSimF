using Distributed
using LinearAlgebra
using Random
using Distributions
using StatsBase, StatsPlots, KernelDensity
using SharedArrays
using DelimitedFiles



"""
Simple particle class
"""
mutable struct particle
    R::Float32 #radius
    m::Float32 #mass
    p::Array{Float32,1} #position
    v::Array{Float32,1} #velocity
end

"Particle Constructor: Constructs a particle at the coordinate origin"
particle(R,m,v) = particle(R,m,[0.0,0.0,0.0],v)
"Copy constructor"
particle(p::particle) = particle(p.R,p.m,p.p,p.v)


struct collisionstatistics
    n_particles # the number of particles in the simulation of the icdfs 
    n_collisions # the number of collisions in the simulation of the icdfs 
    mass_ratios # a list of mass ratios for the individual icdfs
    icdf_data #the icdf collision statistics data: list of individual icdf functions
end

"""
Calculate the post collision velocities (velocities after a fully elastic collision) from 
the positions, masses and velocities of two particles P1 and P2
"""
function postcollisionvelocities(P1,P2)        
    # calculate the normalized connecting vector: 
    e_para = (P1.p - P2.p)/(P1.R+P2.R)
    #normconnectingvector(P1,P2)
    
    # calculate parallel impact velocities
    v_i_para_1 = dot(P1.v,e_para)
    v_i_para_2 = dot(P2.v,e_para)
    
    # calculate perpendicular impact velocities
    v_i_perp_1 = P1.v - v_i_para_1 * e_para
    v_i_perp_2 = P2.v - v_i_para_2 * e_para
    
    # calculate parallel final velocities
    m_total = P1.m + P2.m
    v_f_para_1 = ((P1.m - P2.m)*v_i_para_1 + 2*P2.m*v_i_para_2) / m_total
    v_f_para_2 = ((P2.m - P1.m)*v_i_para_2 + 2*P1.m*v_i_para_1) / m_total    
    
    # calculate final velocities in lab frame
    v_f_1 = v_f_para_1 * e_para + v_i_perp_1 
    v_f_2 = v_f_para_2 * e_para + v_i_perp_2     
    
    return [v_f_1,v_f_2]
end

"""
Find a random collision position in the lab frame for particle P2 for a collision of particle P1 with particle P2. 
    The "collision position" is a random position of P2 so that the particles touch, with given velocities of both particles. 
    The distribution of collision positions essentialy form a half sphere oriented in the collision direction with 
    the sum of the radii of both particles as radius. 
"""
function randomcollisionposition(P1,P2)
    
    R_sum = P1.R + P2.R
    
    # get a random collision in the standard frame 
    #(with a collision direction as y axis in a 2d coordinate system)
    p_impact_stdframe = randomimpactposition(P1.R,P2.R)
    
    vrel = P1.v - P2.v
    vrel_norm = vrel / norm(vrel)
    
    #rotate impact position to lab frame coordinate system
    P_rot = coordinatetransform(p_impact_stdframe, vrel_norm)
    
    #translate found rotated vector to lab frame: 
    p_impact = P_rot + P1.p

    return(p_impact)
end

"""
Find a random collision position in a standard frame for particle P2 for a collision of particle P1 with particle P2. 
    The "collision position" is a random position of P2 so that the particles touch. The standard frame has the y direction 
    as collision axis. 
    The distribution of collision positions essentialy form a half sphere oriented towards the y axis with 
    the sum of the radii of both particles as radius. 
"""
function randomimpactposition(R1,R2)
    Rsum = R1+R2
    
    #normalized impact position:
    impact_R = sqrt(rand(Uniform(0, 1)))
    φ = rand(Uniform(0, π2))
    
    #impactangle = acos(impactradius)
    cp_x = impact_R * sin(φ) * Rsum
    cp_y = impact_R * cos(φ) * Rsum
    
    cp_z = sqrt(1-abs(impact_R)^2)*Rsum
    
    collisionpos = [cp_x, cp_y, cp_z]
    
    return (collisionpos)
end

"""
Rotates a vector "vec" so that newbasevec and the result vector have the same configuration 
as vec and the y basevector of the coordinate system. 
In other words: The position of "vec" in the normal frame is found when "newbasevec" becomes 
the y base vector of the rotated coordinate system.
E.g., a collision position in a standard frame with y as collision direction is transformed to 
the lab frame. 
"""
function coordinatetransform(vec, newbasevec)
    r = norm(newbasevec)
    φ = atan(newbasevec[2], newbasevec[1])
    θ = acos(newbasevec[3] / r)

    cosφ = cos(φ)
    sinφ = sin(φ)
    cosθ = cos(θ)
    sinθ = sin(θ)

    Rx = [1 0 0; 
          0 cosθ -sinθ;
          0 sinθ cosθ]
    
    Rz = [cosφ -sinφ 0; 
          sinφ cosφ  0;
          0 0 1]
    
    resultvec = Rz * Rx * vec
end

"""
Performa a random collision between a tracked particle P1 and an random impactor particle. 
The position and the velocity of the tracked particle is modified to the post collision conditions. 
The impactor velocity is assumed to be Maxwell-Boltzmann distributed. The velocity in every spatial direction 
is a normal distribution with a standard deviation `sigma`

Parameters: 
  + P1: The tracked particle
  + impact_R: Radius of impactor particle
  + impactor_m: Mass of impactor particle
  + sigma: Width / standard deviation of onedimensional velocity distribution of the impactor gas 
  (sigma parameter of the normal distribution)
  + impactor_v_factor: Scaling factor of the threedimensional velocity of the impactor gas  
  + ts_length: Length of the simulation timestep
"""
function dorandomcollision!(P1, impactor_R, impactor_m, sigma, impactor_v_factor, ts_length= 1.0)
    impactor_v = rand(Normal(0.0,sigma),3) * impactor_v_factor
    impactor = particle(impactor_R, impactor_m,impactor_v)
    impactor.p = randomcollisionposition(P1,impactor)
    vf1,vf2 = postcollisionvelocities(P1,impactor)
    P1.p = P1.p + (vf1 * ts_length)
    P1.v = vf1
    
    if isnan(P1.v[1])
        println(impactor_v, " f:",impactor_v_factor, " ", impactor.p, " ", P1.v)
    end
end


"""
Calculate the x value for a given probability from 
an empirical cumulative density function
  + cf: an empirical cumulative density function 
  + prob_val: a probability to return the x value for 
"""
function inverse_ecdf(cf, prob_val)

    cf_vals = cf.sorted_values
    
    if prob_val <= 0.0
        return cf_vals[1]
    end
    
    if prob_val > 1.0
        return cf_vals[end]
    end
    
    len_cf = length(cf_vals)
    
    fi = Int(floor(prob_val * len_cf))
    if fi == 0
        fi = 1
    end
    return cf_vals[fi]
end

"""Writes an inverse empirical cumulative distribution function `ecdf` to a text file
with name `filename`. 
The IECDF is usually written in 0.1% probability steps (1000 steps) but the step number 
is variable (`n_psteps` parameter)
"""
function write_iecdf(ecdf, filename, n_psteps = 1000)
    xsamples = range(0.0, 1, length= n_psteps )
    cf_samples = [inverse_ecdf(ecdf,x) for x in xsamples]
    open(filename, "w") do io
           writedlm(io, cf_samples)
    end
end

"""Writes an simion style mbmr file from a collisionstatistics struct `cs` to a structured text file with name `filename`. 
The ECDFs in cs are written as inverse cumulative density functions in ~0.1% probability steps (1002 steps) by default, 
but the step number is variable (`n_psteps` parameter)
"""
function write_collision_statistics_file(cs::collisionstatistics, filename, n_psteps = 1002)
    xsamples = range(0.0, 1, length=n_psteps )
    icdfs = cs.icdf_data

    open(filename, "w") do io
        write(io, "; MBMR collision statistic array file\n; Written by CollisionStatistics.jl\n")
        write(io, "; n_statistics="*string(length(cs.mass_ratios))*"\n")
        write(io, "; n_particles="*string(cs.n_particles)*"\n")
        write(io, "; n_collisions="*string(cs.n_collisions)*"\n")
        write(io, "; n_dist_points="*string(n_psteps)*"\n")
        for i= 1:length(icdfs)
            df = icdfs[i]
            mr = cs.mass_ratios[i]
            cf_samples = [inverse_ecdf(df,x) for x in xsamples]
            write(io, "\n; ICDF_massratio="*string(mr)*"\n")
            writedlm(io, cf_samples)
        end
    end
end