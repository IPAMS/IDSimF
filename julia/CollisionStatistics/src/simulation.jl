using Plots

"""
Plots a single collision between particles P1 and P2
"""
function plotcollision(P1,P2)
    vf1,vf2 = postcollisionvelocities(P1,P2)
    println(vf1)
    println(vf2)

    plot( [(P1.p-P1.v)[1],P1.p[1]], [(P1.p-P1.v)[2],P1.p[2]],)
    plot!([P1.p[1], (P1.p+vf1)[1]],[P1.p[2], (P1.p+vf1)[2]])
    plot!( [(P2.p-P2.v)[1],P2.p[1]], [(P2.p-P2.v)[2],P2.p[2]],)
    plot!([P2.p[1], (P2.p+vf2)[1]],[P2.p[2], (P2.p+vf2)[2]])
end

"""
Tracks and plots a random collision path of a particle P1. 
The particles in the background gas P1 is colliding with (impactors) have Maxwell Boltzmann distributed 
velocities. 
The timestep length of the simulation is uniformly 1. 

Parameters: 
+ impactor_R: Impactor / background gas radius
+ impactor_m: Impactor / background gas mass
+ sigma: Standard deviation of the Maxwell Boltzmann distribution of the background gas particles 
+ nsamples: The number of collision events 
"""
function tracecollisionpath(P1,impactor_R,impactor_m,sigma,nsamples)
    p1pos = zeros(nsamples,2)
    
    for i in 1:nsamples
        p1pos[i,:] = P1.p        
        dorandomcollision!(P1, impactor_R, impactor_m, sigma)
    end
    Plots.plot(p1pos[:,1],p1pos[:,2],markershape=:circle, markersize = 2)
end


"""
Single Diffusion migration simulation (for one particle) with a background ("impactor") gas with Maxwell Boltzmann
ditributed velocities. The time steps can be  poisson distributed (with a mean legth of 1) or uniformly 1. 

Parameters: 
  + P1: The colliding particle
  + impactor_R: Impactor particles radius
  + impactor_m: Impactor particles mass
  + sigma: Standard deviation of the one dimensional velocity distribution 
  (sigma of the normal distributed velocities in the individual spatial dimensions)
  + ncollisions: number collisions to simulate
  + mode: either :poisson for poisson distributed timesteps or :normalized for timesteps with length 1
"""
function sim_diffusionmigration(P1,impactor_R,impactor_m,sigma,ncollisions::Int64,mode)  
    ncollisions_performed = 1
    velocities = zeros(ncollisions)
    freepaths = zeros(ncollisions)
    impactor_v_factor = sqrt(P1.m / impactor_m)
    m_ratio = sqrt(P1.m / impactor_m)
    while ncollisions_performed < ncollisions
        
        #get number of collisions in this time step: 
        if mode == :poisson
            ncollisions_ts = rand(Poisson())
            if (ncollisions_ts + ncollisions_performed) > ncollisions 
                ncollisions_ts = ncollisions - ncollisions_performed 
            end
        elseif mode == :normalized
            ncollisions_ts = 1
        end
            
        
        if ncollisions_ts == 0
            P1.p = P1.p + P1.v
        elseif ncollisions_ts == 1
            #impactor_v_factor = m_ratio * norm(P1.v)
            dorandomcollision!(P1, impactor_R, impactor_m, sigma, impactor_v_factor, 1.0)
            velocities[ncollisions_performed] = norm(P1.v)
            freepaths[ncollisions_performed] = velocities[ncollisions_performed]
            ncollisions_performed += 1
        else
            ts_length = 1.0/ ncollisions
            for i=1:ncollisions_ts
                #impactor_v_factor = m_ratio * norm(P1.v)
                dorandomcollision!(P1, impactor_R, impactor_m, sigma, impactor_v_factor, ts_length)
                velocities[ncollisions_performed] = norm(P1.v)
                freepaths[ncollisions_performed] = velocities[ncollisions_performed] * ts_length
                ncollisions_performed += 1
            end
        end
    end
    
    return [P1, velocities, freepaths]
end

"""
Simulates the distribution of diffusion migration of an ensemble of particles and returns the final states of 
the individual simulated particles after the diffusion migration. 

The the velocities of the background gas ("impactor" gas) particles are assumed to be Maxwell Boltzmann distributed. 
The time step lengths of the simulation can be Poisson distributed (with a mean of 1) or uniformly with 
length 1. 

Parameters: 
  + P_start: A particle which is the initialization of the simulated particles
  + impactor_R: Impactor particles radius
  + impactor_m: Impactor particles mass
  + sigma: Standard deviation of the one dimensional velocity distribution 
  (sigma of the normal distributed velocities in the individual spatial dimensions)
  + ncollisions: number of collisions to simulate 
  + nparticles: number of particles to simulate 
  + plot_kde: if true, a kernel density estimation of the resulting distribution is plotted
  + mode: either :poisson for poisson distributed timesteps or :normalized for timesteps with length 1
"""
function sim_diffusionmigrationdist(P_start,impactor_R,impactor_m,sigma,nparticles,ncollisions,plot_kde;mode= :poisson)
    finalpositions = SharedArray{Float64}(nparticles,6)
    @sync @distributed for i = 1:nparticles
        P = particle(P_start)
        P.v = rand(Normal(0.0,sigma),3) #distribute the ion with a MB velocity
        P, velocities, freepaths = sim_diffusionmigration(P,impactor_R,impactor_m,sigma,ncollisions,mode)
        finalpositions[i,1:3] = P.p
        finalpositions[i,4] = norm(P.p)
        finalpositions[i,5] = mean(velocities)
        finalpositions[i,6] = mean(freepaths)
    end
    
    #get kernel density estimation:
    if plot_kde
        dens = kde((finalpositions[:,1],finalpositions[:,2]))
        plt = plot(dens)
        display(plt)
        dens = kde((finalpositions[:,2],finalpositions[:,3]))
        plt = plot(dens)
        display(plt)
    end
    return finalpositions
end

"""
Calculates the r50% radius, which is the radius with 50% probability density with background ("impactor") gas 
with radius 1 and mass 1. 

+ P1: The colliding / tracked particle 
+ ncollisions: number collisions to simulate
+ sigma: Standard deviation of the one dimensional velocity distribution 
(sigma of the normal distributed velocities in the individual spatial dimensions)
+ mode: either :poisson for poisson distributed timesteps or :normalized for timesteps with length 1
"""
function sim_run_r50percent(P1, ncollisions,sigma,nparticles;mode=:poisson)
    f_pos = sim_diffusionmigrationdist(P1, 1, 1, sigma,nparticles,ncollisions,false,mode=mode)
    ecdf = StatsBase.ecdf(f_pos[:,4])
    return inverse_ecdf(ecdf,0.5)
end

"""
Calculate and export collision statistics as inverse cumulative density function
(e.g. as input parameter in SDS). 
The collision induced diffusive drift of number of particles is simulated and an 
empirical inverse cumulative density function is calculated from the results. 

Parameter: 
  + filebasename: basename pattern of the result files
  + nparticles: number of particles to simulate
  + ncollisions: number of collisions to simulate
  + ionmasses: Array of impactor masses. The simulation is run for every mass in this list and an 
  individual result file is exported
"""
function calculate_collision_statistics(filebasename, nparticles = 3000, ncollisions = 100000, 
                                        ionmasses = [1,10,100,1000,10000]; filemode = :mbmrfile) 

    result_ecdfs = []
    
    for im in ionmasses
        println(im)
        P1 =particle(0.5, im,  [0,0,0], [1,0,0])
        f_pos = sim_diffusionmigrationdist(P1, 0.5, 1, sigma_maxwellboltzmann, nparticles, ncollisions, false; mode=:poisson)
        ecdf = StatsBase.ecdf(f_pos[:,4])
        filename = filebasename*string(im)*".txt"

        if filemode == :singlefiles
            write_iecdf(ecdf,filename)
        end
        if filemode == :mbmrfile
            append!(result_ecdfs, [ecdf])
        end
    end

    if filemode == :mbmrfile
        cs = collisionstatistics(nparticles, ncollisions, ionmasses, result_ecdfs)
        write_collision_statistics_file(cs, filebasename*".dat")
    end

end