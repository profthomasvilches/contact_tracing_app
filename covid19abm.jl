module covid19abm

using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@enum HEALTH SUS LAT PRE ASYMP INF REC DED UNDEF 

# TO DO: remove fields not used, comment on others, rename variable to be more description
Base.@kwdef mutable struct Human
    idx::Int64 = 0 # Each individual is identified by an index ranging from 1 to 10000.
    health::HEALTH = SUS # Records the health status of a given individual. 
    health_status::HEALTH = SUS # Records the health status of a given individual. 
    swap::HEALTH = UNDEF # Records the health status that an individual adopts after leaving their current state.
    swap_status::HEALTH = UNDEF # Records the health status that an individual adopts after leaving their current state.
    sickfrom::HEALTH = UNDEF # Records the health status of the individual that infected a given person.
    sickby::Int64 = -1 # Records the index of the individual that infected a given person.
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16 = 0 # in years. don't really need this but left it incase needed later
    ag::Int16 = 0 # Identifies the age bracket to which an individual belongs.
    tis::Int16 = 0   # time in state
    exp::Int16 = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16 = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)   
    wentto::Int8 = 0 # Used to compute the reduced sensitivity of a test for asymptomatic individuals. This is done in matrices_code.jl.
    incubationp::Int16 = 0 # Records the incubation period of an indivdual.
    got_inf::Bool = false # Tracks whether an individual is infected.
    recovered::Bool = false  # what is this? why not use the health_status field? 
    # The health_status is updated to REC when an individual recovers. The recovered field simply indicates which individuals have gotten over the infection. It plays no role in the simulations but could be useful for local testing.
    daysisolation::Int64 = 999 # Records the number of days that an individual spends in isolation.
    daysinf::Int64 = 0 # Tracks the amount of time that an individual has been infected.
    totaldaysiso::Int32 = 0 # Tracks the number of days that an individual spends in isolation.
    has_app::Bool = false # Records whether an individual downloads the application
    contacts::Vector{Vector{Int64}} = [[0; 0]] # Records the indexes of the app-using contacts of an individual over the past three days.
    ncontacts_day::Int8 = 0 # Records the number of app-using contacts that an individual accumulates on a given day.
    testedpos::Bool = false # Records whether an individual has tested positive for the virus.
    notified::Bool = false # Records whether an individual has received a notification indicating that they were in close proximity to an infected person.
    timetotest::Int64 = 9999 # The number of days remaining before an individual tests for the infection.
    n_tests_perf::Int64 = 0 # Records the number of tests that an indiviudal performs.
    time_since_testing::Int64 = 0 # The amount of time that has passed since an individual last tested.
    n_neg_tests::Int64 = 0 # Records the number of negative tests that an individual has received for a given notification.
    quar::Int64 = 0 # Records a 1 if an individual isolates. This is used to write the totalquar.dat file that records the number of isolations that take place during each simulation.
    symp_inf::Bool = false # what is this? why not use the health_status field?
    # This field was used in local testing to record whether an individual became symptomatic. The health_status field was not used because of the potential for it to change to either REC or DED if a person recovers or dies, respectively.
    reported::Bool = false # Records whether an individual reported their positive test results on the application
    tstd::Bool = false # Used to record whether an individual tested at some point during the simulation. This was used in local testing and has no effect on the simulation.
    isolation_days::Int64 = 8 # Assigns the minimum number of days that an individual can spend in isolation. Symptomatic individuals isolate at the onset of symptoms for 5 days, while those without symptoms isolate for 7. Should an individual develop symptoms in isolation, they remain for another 5 days after this event.
    # The reason why isolation_days is set to 8 as opposed to 7 is that the number of days in isolation increases by 1 in the same time step that they are isolated. This approach ensures that individuals are isolated for 7 time steps and released on the eighth. For newly isolated symtomatic individuals, at the end of the time step, they have recorded 0 days in isolation. They remain in isolation up to and including daysisolation = 4, meaning that they isolated for five time steps, and are released on the 6th when daysiolation = 5. 
    # The latter scenario occurs because daysisolation increases by 1 in time_update before move_to_inf and the internal summoning of _set_isolation_ are called.
    inf_asp::Bool = false # Records whether an individual was infected at any point during the simulation. Used in local testing and should have no effect on the simulations.
end



# TO DO: comment on fields, rename variable to be more description
# if a field is only used inside a function, make it constant inside that function
# if a field remains constant, pull it out and make it a global const variable
# the only fields that should be here are the ones used for scenario analysis (i.e., ones that we need to change)
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.1 # This is the transmission probability of the virus. It is altered based on the health status of an individual in the _get_betavalue function.
    popsize::Int64 = 10000 # We define the size of the population being considered in the simulation.
    prov::Symbol = :ontario # We use this scenario to define the age distribution of the population using the get_province_ag and initialize functions.
    modeltime::Int64 = 365 # This parameter defines the duration of the simulation in days.
    initialinf::Int64 = 1 # This parameter defines the number of infected individuals when the simulation starts.
    start_testing::Int64 = 1 # This is the day that testing starts in the simulation.
    test_for::Int64 = 365 # The number of days that testing is performed during the simulation.
    file_index::Int16 = 0 # Can be used to distinguish between file names when results are sent to the cluster.
    app_coverage = 1.0 # The proportion of the eligible population that downloads the application.
    time_until_testing::Int64 = 1 # The number of days between the receipt of a notification and performance of a test. 
    n_tests::Int64 = 2 # The maximum number of tests that an individual can perform for a single notification.
    testing::Bool = false # This field dictates whether testing is in effect.
    not_swit::Bool = false # This field indicates whether app-using individuals that tested positive send notifications to their contacts during the simulation. 
    num_sims::Int64 = 500 # This is the number of simulations being performed.
    iso_con::Int64 = 0 # This refers to the number of contacts that isolated individuals encounter. If the value is 0, there will be no interaction and a random number of contacts otherwise. 
    test_sens::Int64 = 1 # This is the sensitivity of the test. If the value is 1, infected individuals always test positive and otherwise, the get_matrix function, defined in matrices.jl and called in matrices_code.jl is used to determine the rate.
    comp_bool::Bool = true # This controls the probability of individuals complying with the notification regulations. If true, individuals always comply and if false, individuals have a 50% chance of obeying. 
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z) # Displays the HUMAN struct.

# The probabilities of individuals testing positive is determined using the matrices.jl and matrices_code.jl scripts.
include("matrices_code.jl")
include("matrices.jl")



## constants 
const humans = Array{Human}(undef, 0) # This array is used to define the unqiue characteristics of each individual.
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99] # We define our 5 age brackets.
const time_between_tests = 3 # The number of days that pass between the first and second test that an individual takes for a single notification.
export ModelParameters, HEALTH, Human, humans, time_between_tests

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix, hh1, nra = main(ip,simnum) 

    #Get the R0    
    R01 = length(findall(k -> k.sickby in hh1,humans))/length(hh1)

    all1 = _collectdf(hmatrix) # This matrix records information related to the incidence and prevalence of infections throughout the population.

    age_groups = [0:14, 15:24, 25:34, 35:44, 45:54, 55:64, 65:999] # We separate the population into 7 age groups, each of which will be used to generate files relevant to infection.

    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    
    # We proceed to divide the population into the age groups.
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    ag7 = _collectdf(spl[7])

    # We divide the population into two groups: those that use the application and those that do not.
    v_has = [x.has_app ? 1 : 2 for x in humans]
    spl = _splitstate(hmatrix, v_has)
    use_1 = _collectdf(spl[1])
    use_2 = _collectdf(spl[2])

    # We prepare the following dataframes to generate the desired .dat files in the cluster.
    insertcols!(all1, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum); insertcols!(use_1, 1, :sim => simnum); insertcols!(use_2, 1, :sim => simnum);  

    # We record the indexes of deceased individuals at the end of the simulation. The number of deceased individuals for each age are stored in vector_ded and this will be used to create the year_of_death.dat files in simulations_cluster.jl.
    pos = findall(y-> y == 7, hmatrix[:,end])
    vector_ded::Vector{Int64} = zeros(Int64,100)
    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    ### total days of isolation per age group
    geniso_gr = map(y->findall(x-> x.age in y,humans),age_groups)
    giso = map(y-> sum([ii.totaldaysiso for ii in humans[y]]),geniso_gr)

    quar_tot = sum([Int(x.quar) for x in humans]) # We record the total number of isolation for a given simulation. This variable will be used to generate the totalquar.dat files.

    return (a=all1, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7,
    vector_dead=vector_ded,nra=nra,R0 = R01, giso = giso, u1 = use_1, u2 = use_2, 
    quar_tot = quar_tot)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726) # We set the random seed.
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given") # Ensuring that a valid population has been given.
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime) # Records the health status of each individual every day of the simulation.
    initialize() # initialize population
    
    # split population in agegroups 
    grps = get_ag_dist()

    nra::Vector{Int64} = zeros(Int64,p.modeltime) # We use this vector to run the time loop via time_update, which propels the simulation through the various time points.

    #insert one infected in the latent status in age group 4
    insert_infected(LAT, p.initialinf, 4)

    h_init1 = findall(x->x.health_status  in (LAT,INF,PRE,ASYMP), humans) # Records the indexes of individuals infected at the beginning of the simulation.
    
    # We generate the interactions in which each individual participates the following day.
    for x in humans
        get_nextday_counts(x)
    end
    
    # distributing app_coverage
    dist_app(humans, p, p.num_sims)  

    # start the time loop relating to the period of the simulation before testing protocols are implemented
    for st = 1:min((p.start_testing-1),p.modeltime)
        for x in humans
            if x.iso && !(x.health_status == DED) # Increasing the count for the number of days that an individual is isolated, provided that they are alive.
                x.totaldaysiso += 1
            end
        end

        # Used to create hmatrix, which contains the health status of each individual, recorded every day of the simulation.
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        
        dyntrans(st, grps,sim) # This function determines the number of contacts of each individual based on their age.
        
        sw = time_update(st,sim) ###update the system
    end
    
    setfield!(p,:testing,true) # We initiate testing for the simulation.    

    # start the time loop for the period in which testing is being implemented
    for st = p.start_testing:min((p.start_testing+p.test_for-1),p.modeltime)
        for x in humans

            if x.iso && !(x.health_status == DED) # Increasing the count for the number of days that an individual is isolated, provided that they are alive.
                x.totaldaysiso += 1
            end
        end

        # Used to create hmatrix, which contains the health status of each individual, recorded every day of the simulation.
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop

        dyntrans(st, grps,sim) # This function determines the number contacts of each individual based on their age.

        sw = time_update(st,sim) ###update the sys_time

        nra[st] += sw

        # end of day
    end

    setfield!(p,:testing,false) # We stop testing after the designated period.

    # We run the time loop for the period after the testing period has ended.
    for st = (p.start_testing+p.test_for):p.modeltime
        for x in humans

            if x.iso && !(x.health_status == DED) # Increasing the count for the number of days that an individual is isolated, provided that they are alive.
                x.totaldaysiso += 1
            end
        end

        # Used to create hmatrix, which contains the health status of each individual, recorded every day of the simulation.
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop

        dyntrans(st, grps,sim) # This function determines the number contacts of each individual based on their age.

        sw = time_update(st,sim) ###update the system
        
        # end of day
    end
    return hmatrix, h_init1, nra
end
export main


# The dist_app function distributes the application across the eligible population.
function dist_app(humans, p, sim)
    ageintapp::Vector{Int64} = [18; 65] # The age range for the application.

    rng = MersenneTwister(2421*sim) # We create a simulation-dependent random number generator, which will be used to determine the individuals that receive the application.

    pos = findall(x->x.age in ageintapp[1]:ageintapp[2], humans)
    pos = shuffle(rng, pos)
    pos = pos[1:Int(round(p.app_coverage*length(pos)))]

    for i in pos
        humans[i].has_app = true
    end
end


function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    resize!(humans, p.popsize) # Creates the population.
end
export reset_params, reset_params_default


## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state


# This function converts the incidence and prevalence of individuals at each time step into a dataframe.
function _collectdf(hmatrix) 
    ## takes the output of the humans x time matrix and processes it into a dataframe
    
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev) 
    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix) # This function stores the incidence and prevalence in arrays.
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol) # Determines the incidence based on health statuses and time points.
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findall(x-> r[x] == inth && r[x] != r[x-1],2:length(r))
        idx = idx .+ 1
        if idx !== nothing
            for i in idx 
                timevec[i] += 1
            end
        end
    end
    return timevec
end

function _get_column_prevalence(hmatrix, hcol) # Determines the prevalence based on health statuses and time points.
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence

## initialization functions 
function get_province_ag(prov) # Assists in the determination of the age distribution of the simulation based on the given location.
    ret = @match prov begin
       #= :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])
        :newyorkcity   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])=#
        :ontario => Distributions.Categorical(@SVector [0.04807822, 0.10498712, 0.12470340, 0.14498051, 0.13137129, 0.12679091, 0.13804896, 0.10292032, 0.05484776, 0.02327152])
        :canada => Distributions.Categorical(@SVector [0.04922255,0.10812899,0.11792442,0.13956709,0.13534216,0.12589012,0.13876094,0.10687438,0.05550450,0.02278485])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag


function initialize() # We initialize the population.
    agedist = get_province_ag(p.prov)
    agebraksnew = [0:4,5:14,15:24,25:34,35:44,45:54,55:64,65:74,75:84,85:99]
    
    track_days::Int8 = 3 # The number of preceding days for which the contacts of an infected individual are considered.

    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        agn = rand(agedist)
        x.age = rand(agebraksnew[agn]) 
        x.ag = findfirst(y-> x.age in y, agebraks)
        x.exp = 999  ## susceptible people don't expire.
        x.time_since_testing = time_between_tests

        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        x.contacts = repeat([[0]], track_days)
    end
end
export initialize

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false) 
        @inbounds for i in h 
            x = humans[i]
            x.dur = sample_epi_durations(x) # Determing the amount of time that an individual spends in each state.

            if health == PRE
                    x.swap = health
                    x.swap_status = PRE
                    x.daysinf = x.dur[1]+1
                    x.wentto = 1
                    move_to_pre(x) 

            elseif health == LAT
                    x.swap = health
                    x.swap_status = LAT
                    x.daysinf = rand(1:x.dur[1])
                    move_to_latent(x)

            elseif health == INF
                    x.swap = health
                    x.swap_status = INF
                    x.wentto = 1
                    move_to_inf(x)

            elseif health == ASYMP
                    x.swap = health
                    x.swap_status = ASYMP
                    x.wentto = 2
                    move_to_asymp(x)
                
            elseif health == REC 
                    x.swap = health
                    x.swap_status = REC
                    x.wentto = rand(1:2)
                    move_to_recovered(x)

            else 
                    error("can not insert human of health $(health)")
            end

            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function time_update(st,sim) # This function updates the system at each time point.
    # counters to calculate incidence

    nra::Int64 = 0

    test_ra = 2 # This signifies that Abbott_PanBio tests are being performed. (1 - RT-PCR, 2 - Abbott_PanBio, 3 - BTNX_Rapid_Response, 4 - Artron)    
    
    if p.testing 
        for x in humans 
            x.timetotest -= 1
            x.time_since_testing += 1 # Taiye: We could measure this in days.

            if x.notified && !x.testedpos && x.n_tests_perf <= p.n_tests && x.timetotest <= 0 && x.time_since_testing >= time_between_tests # These are the conditions required for testing to take place: The individual must be notified either by the app or because they are symptomatic and use the application, they have not tested positive yet, the maximum number of tests has not been reached, timetotest is immediately and the time since testing is at least as long as the 3 day period between the taking the first and second test. 
                testing_infection(x, test_ra)                
                x.time_since_testing = 0
                x.n_tests_perf += 1
                if x.n_tests_perf == p.n_tests # If the maximum number of tests has been reached, individuals will not test again unless they receive another notification.
                    x.notified = false
                    x.n_tests_perf = 0
                end
            end
        end
    end
    
    rng = MersenneTwister(246*st*sim) # Creating a random number generator to for the send_notification function.
    for x in humans
        if x.testedpos && !x.reported # Individuals will send notifications if they test positive and have not communicated their contraction before.
                send_notification(x,p.not_swit,st,sim,rng)
        end
    end

    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
         
        x.daysinf += 1
        x.daysisolation += 1

        

        if x.tis >= x.exp # If the time in state reaches the duration that individual is supposed to spend in a state, their health status changes.
            @match Symbol(x.swap_status) begin
                :LAT  => begin 
                    move_to_latent(x); 
                end
                :PRE  => begin move_to_pre(x); end
                :ASYMP => begin move_to_asymp(x);  end              
                :INF  => begin move_to_inf(x); end # Taiye (2025.07.04)
                :REC  => begin move_to_recovered(x); end
                :DED  => begin move_to_dead(x); end
                _    => begin dump(x); error("swap expired, but no swap set."); end
            end
        end

        #if the individual recovers, we need to set they free. This loop must be here
        if x.iso && x.daysisolation >= x.isolation_days && !(x.health_status == DED) 
            _set_isolation(x,false,:null) 
            
            x.n_tests_perf = 0 # Taiye
            x.n_neg_tests = 0 # Taiye            
        end

        # get the meet counts for the next day 
        get_nextday_counts(x)
        
    end

   return nra
end
export time_update


@inline function _set_isolation(x::Human, iso, via) # This functions places and releases individuals from isolation.
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    if x.isovia == :null || via == :symp 
        x.iso = iso 
        x.isovia = via
        x.daysisolation = 0

        x.isolation_days = via == :test ? 8 : 5 # Taiye (2025.10.14)

    elseif !iso
        x.iso = iso 
        x.isovia = via
    end

    # Taiye (2025.08.05):
    if x.iso == true
        x.quar = 1
    end

end
function sample_epi_durations(y::Human)
    # when a person is sick, samples the duration of their time in a given state

    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # The latent distribution is truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)# The presymptomatic distribution is truncated between 0.8 and 3

    asy_dist = Gamma(5, 1) # We compute the asymptomatic distribution.
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2) # We compute the symptomatic distribution.

    # We narrow the distributions down so that they can be assigned at the individual level.
    latents = Int.(round.(rand(lat_dist)))
    y.incubationp = latents
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = max(Int.(ceil.(rand(asy_dist))),1)
    infs = max(Int.(ceil.(rand(inf_dist))),1)
    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration

    x.inf_asp = true

    x.health = x.swap
    x.health_status = x.swap_status
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] # Age-specific values used to determine whether latent individuals enter the asymptomatic or presymptomatic state.
    
    age_thres = [4, 19, 49, 64, 79, 999] # The oldest age in each group being considered.
  
    g = findfirst(y-> y >= x.age, age_thres) # Used to sync the age thresholds with their respective brackets.
 
    if rand() < (symp_pcts[g])
        x.swap = PRE
        x.swap_status = PRE
        x.wentto = 1
    else
        x.swap = ASYMP
        x.swap_status = ASYMP
        x.wentto = 2
    end
    
    x.got_inf = true # Infected individuals always enter the latent stage first.
    ## in calibration mode, latent people never become infectious.
    
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = x.swap  
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
   
    x.swap = REC
    x.swap_status = REC
    
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    # This function moves individuals to the presymptomatic state.

    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period

    x.swap = INF
    x.swap_status = INF    
end
export move_to_pre

function testing_infection(x::Human, teste)
    x.tstd = true
    pp = _get_prob_test(x,teste,p.test_sens) # This is the probability is produced by the matrices_code.jl and matrices.jl scripts.

    if rand() < pp # This is the likelihood of an infected individual testing positive.
       
        x.testedpos = true
        if !x.iso # Provided that an individual is not already in isolation, they enter after testing positive. We included this condition to avoid restarting the isolation period.
            _set_isolation(x, true, :test)
        end
    else # Taiye: counting the number of negative tests performed.
          x.n_neg_tests += 1
    end
end

function send_notification(x::Human,switch,st,sim,rng) 
    if switch # The function only performs its objectives if notifications are turned on.
      
        v = vcat(x.contacts...) # Prepares the indexes of the app-using contacts recorded by the individual over the last three days to receive notifications.
       
        for i in v
            if 1 <= i <= length(humans) && !humans[i].notified # Taiye: To avoid new notifications resetting times.
               
                humans[i].notified = rand(rng) < round(0.5,digits = 1) && !p.comp_bool ? false : true # These conditions are only relevant if we are considering 50% compliance since individuals in essence flip a coin to determine whether they will respond to a notification appropriately.
          
                humans[i].timetotest = p.time_until_testing # We set the delay between the receipt of a notification and performance of a test.
            end
        end
        x.reported = true # We record the fact that an individual reported the infection.
    end

end

function move_to_inf(x::Human) # Transfers individuals to the symptomatic stage.

    x.symp_inf = true # Recording that the individual was symptomatic at some point during the simulation. Mainly used in local testing.

    ## for swap, check if person will die, or recover
 
    # We define the age groups being considered, which will aid in the determination of whether an individual recovers or passes away. 
    groups = [0:34,35:54,55:69,70:84,85:100] 
    gg = findfirst(y-> x.age in y,groups)

    # These elements will be used to determine whether individuals recover or pass away.
    θ = (0.05, 0.1, 0.15, 0.4, 0.8) # This array was originally used to store the rates at which individuals transitioned from the presymptomatic to severe class. We repurpose it here to determine the death and symptomatic recovery rates.
    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
    mh_2 = round.(θ.*mh,digits=5) 

    # Individuals are registered as symptomatic with an undefined future progression.
    x.health = x.swap
    x.health_status = x.swap_status
    x.swap = UNDEF
    
    x.tis = 0 # The time in state is set to 0 days.
    
    if p.testing && !x.testedpos && x.has_app # We assume that while testing is in place, symptomatic app-users will test for the virus, even if they did not receive a notification from a contact.

        x.notified = true

        x.timetotest = 0 # Taiye (2025.10.21): Infected individuals test as soon as possible after the onset of symptoms.
    end

    _set_isolation(x, true, :symp) # We assume that symptomatic individuals always isolate, regardless of whether they use the app or test positive.
       
    # We determine whether an individual recovers or passes away.
    if rand() < mh_2[gg] 
            x.exp = x.dur[4] 
            x.swap = DED
            x.swap_status = DED
    else 
            x.exp = x.dur[4]  
            x.swap = REC
            x.swap_status = REC
    end
end


function move_to_dead(h::Human) # Symptomatic individuals have the potential of passing away.
    h.health = h.swap
    h.health_status = h.swap_status
    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
end

function move_to_recovered(h::Human) # Asymptomatic individuals always recover, while symptomatic individuals do so with a certain probability.
    h.health = h.swap
    h.health_status = h.swap_status

    h.recovered = true # Recording that an individual got over the virus. This has no effect on the simulation but could be useful in local testing.

    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    h.got_inf = false 
end


@inline function _get_betavalue(xhealth) # We determine the relative transmissibility, β.

    frelasymp = 0.26 ## relative transmission of asymptomatic
    
    bf = p.β # Presymptomatic individuals have a relative transmissibility equal to p.β, while the value differs for asymptomatic and symptomatic cases.

    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * frelasymp 
    elseif xhealth == INF 
        bf = bf * 0.89
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human) # We determine the number of contacts that an individual encounters the next day.
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday   
    
    contact_change_rate = 1.0 # The first contact change rate to compute the number of recorded contacts the following day.
    
    contact_change_2 = 1.0 # The second contact change rate to compute the number of recorded contacts the following day.
    
    cnt = 0 # We initialize the number of contacts on a given day to be 0.

    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    
    if !x.iso # We first determine the number of contacts of indiviudals that are not in isolation.
        aux =  contact_change_rate*contact_change_2

        cnt = rand(negative_binomials(ag,aux)) ##using the contact average for shelter-in
        
        x.nextday_meetcnt = cnt # We determine the number of contacts that an individual has the following day.
    
    elseif !(x.health_status == (DED)) # We determine the contacts of individuals in isolation, provided that they are not deceased.

        # If p.iso_con = 0, we assume that isolated individuals do not have any human interaction and otherwise, we use the negative_binomials_shelter function to determine the number of contacts.
        if p.iso_con == 0
            cnt = p.iso_con
        else
            cnt = rand(negative_binomials_shelter(ag,contact_change_2))
        end
      
        x.nextday_meetcnt = cnt # We determine the number of contacts that an individual has the following day.
    end
    
   # We assume that deceased individuals do not have contacts.
    if x.health_status == DED
        x.nextday_meetcnt = 0
    end

    # We record the contacts of an individual recorded over the last three days, which are separated by the indexes of the contacts field by day.
    # The contacts from each preceding day are passed to the following when time updates and contacts[1] is initialized as an array of zeros.
    x.contacts[3] = deepcopy(x.contacts[2])
    x.contacts[2] = deepcopy(x.contacts[1])
    x.contacts[1] = repeat([0], max(x.nextday_meetcnt, 1))
    x.ncontacts_day = 0

    return cnt
end


function dyntrans(sys_time, grps,sim) # We use this function to prepare the population for the formulation of contacts.

    pos = shuffle(1:length(humans)) # We randomize the order in whihc individuals will be assessed.
    # go through every infectious person

    for x in humans[pos] # Iterating through the population using the new order.
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts
           
            #general population contact
            gpw = Int.(round.(cm[x.ag]*cnts)) # We use the contact matrix, cm and the relevant age distribution to determine the contacts of the indvidual.

            #cnts = number of contacts on that day

            perform_contacts(x,gpw,grps,xhealth) # We determine the contacts that an indiviudal encounters using the perform_contacts function.                    
    end
end
export dyntrans

function perform_contacts(x,gpw,grp_sample,xhealth)
    for (i, g) in enumerate(gpw) 
        meet = rand(grp_sample[i], g)   # sample the people from each group
       
        # go through each person
        for j in meet 
            y = humans[j]
            ycnt = y.nextday_meetcnt  
            xcnt = x.nextday_meetcnt 
            
            # We remove a contact to avoid an individual sending a notification to the person that initially informed them.
            y.nextday_meetcnt = y.nextday_meetcnt - min(1,ycnt) 
            x.nextday_meetcnt = x.nextday_meetcnt - min(1,xcnt) 

            # We ensure that x and y recorded contacts.
            xcnt == 0 && continue 
            ycnt == 0 && continue 
            
            y.idx == x.idx && continue # verifies that indexes are not the same

            # Assuming that both x and y have the app, they will record each other as contacts, increasing the number of contacts that they obtained on this day.
            if y.has_app && x.has_app
                x.ncontacts_day = x.ncontacts_day+1
                if length(x.contacts[1]) >= x.ncontacts_day
                    x.contacts[1][x.ncontacts_day] = y.idx
                end
                y.ncontacts_day = y.ncontacts_day+1
                if length(y.contacts[1]) >= y.ncontacts_day
                  y.contacts[1][y.ncontacts_day] = x.idx
                end
            end
            
            # We determine the relative transmissibility of the infected individual in an interaction.
            if x.health_status in (PRE, ASYMP, INF) && y.health == SUS && y.swap == UNDEF 
                beta = _get_betavalue(x.health)
                if rand() < beta
                
                    y.exp = y.tis   ## force the move to latent in the next time step.

                    y.sickfrom = x.health ## stores the infector's status to the infectee's sickfrom

                    y.sickby = y.sickby < 0 ? x.idx : y.sickby # Recording the index of the infector.

                    y.swap = LAT
                    y.swap_status = LAT
                    y.daysinf = 0
                    y.dur = sample_epi_durations(y)
                end  
            elseif y.health_status in (PRE, ASYMP, INF) && x.health == SUS && x.swap == UNDEF
                beta = _get_betavalue(y.health)
                if rand() < beta
                
                    x.exp = x.tis   ## force the move to latent in the next time step.

                    x.sickfrom = y.health ## stores the infector's status to the infectee's sickfrom
                  
                    x.sickby = x.sickby < 0 ? y.idx : x.sickby # Recording the index of the infector.

                    x.swap = LAT
                    x.swap_status = LAT
                    x.daysinf = 0
                    x.dur = sample_epi_durations(x)
                end            
            end        
        end
    end  

    x.nextday_meetcnt = 0 # We reset this field to 0 for the following time step.

end

function contact_matrix() # We generate the contact matrix, which is divided according the age of an individual.
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
    CM[1] = [0.25,0.132,0.44,0.144,0.034]
    CM[2] = [0.0264,0.43,0.404,0.108,0.0316]
    CM[3] = [0.03,0.13,0.602,0.179,0.059]
    CM[4] = [0.026,0.086,0.456,0.3,0.132]
    CM[5] = [0.012,0.052,0.303,0.266,0.367]  
   
    return CM
end

function negative_binomials(ag,mult)  # We use the negative_binomials function to determine the number of contacts that individuals obtain on a given day using age distributions.
    ## the means/sd here are calculated using _calc_avgag
    # [0:4, 5:19, 20:49, 50:64, 65:99]
    means = [6.97;9.54;10.96;8.05;4.41]
    sd = [5.22;6.66;8.35;6.86;3.83]
    means = means*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]
end

const cm = contact_matrix() # We define the contact matrix as a constant.

export negative_binomials


function negative_binomials_shelter(ag,mult) # We use the negative_binomials function to determine the number of contacts that individuals obtain in isolation on a given day using age distributions.
    ## the means/sd here are calculated using _calc_avgag
    #72% reduction
    means = [1.95; 2.67;3.07; 2.255;1.234]
    sd = [1.461518;1.863781;2.337172;1.920688;1.12]
    means = means*mult
    #sd = sd*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]   
end
end