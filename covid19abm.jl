module covid19abm

# Thomas:
# - parameter how many tests someone will perform after notification (if negative)
# - if someone tested negative, they will test again and again until the number is reached or is positive
# - be careful: new notification cannot set the times to zero if someone is in a series of testing

# Edit: 2025.06.24
# Any edits that I make will include "#Taiye:".

# Taiye (2025.05.27):
import Pkg
Pkg.add("Parameters")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("StaticArrays")
Pkg.add("Match")
Pkg.add("Random")

using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
# @enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED UNDEF
@enum HEALTH SUS LAT PRE ASYMP INF REC DED UNDEF # Taiye: MILD MISO IISO HOS ICU
# Taiye: I commented out hospital and ICU patients because my understanding is that the tests
# should be self-administered. Seyed also pointed out that their inclusion would overcomplicate 
# the model. Consequently, MILD, MISO and IISO are also removed because the severity of the 
# symptoms is not relevant.
# Taiye: If we assume that the disease can be contracted more than once, recovered individuals
# could be included in the symptomatic class and otherwise, removed from the population. We could
# create a REM (removed) class for these two groups.
# Taiye (update): Recovered individuals will still be considered contacts. We also want to record the individuals that are deceased.

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    health_status::HEALTH = SUS
    swap::HEALTH = UNDEF
    swap_status::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    

    #comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    # Taiye: We are not considering comorbidities at this stage.

    #vac_status::Int8 = 0 ##
    wentto::Int8 = 0
    incubationp::Int16 = 0

    got_inf::Bool = false
    
    # herd_im::Bool = false
    # Taiye: We are not considering herd immunity at this stage.

    ag_new::Int16 = -1
    
    # days_vac::Int16 = -1
    # Taiye: We are not considering vaccinations at this stage.

    #strain::Int16 = -1
    index_day::Int16 = 1
   
    recovered::Bool = false 
    # vaccine::Symbol = :none
    # vaccine_n::Int16 = 0
    # protected::Int16 = 0
    # boosted::Bool = false
    # n_boosted::Int8 = 0
    # recvac::Int8 = 0 # 1 - rec , 2 - vac ... this field shows which immunity will be used for protection

    # vac_eff_inf::Array{Array{Array{Float32,1},1},1} = [[[0.0]]]
    # vac_eff_symp::Array{Array{Array{Float32,1},1},1} = [[[0.0]]]
    # vac_eff_sev::Array{Array{Array{Float32,1},1},1} = [[[0.0]]]

    # workplace_idx::Int64 = -1

    #### for testing

    daysisolation::Int64 = 999 

    # days_after_detection::Int64 = 999
    # positive::Bool = false

    # Taiye (2025.06.12):
    days_after_detection::Int64 = 0

    daysinf::Int64 = 999 

    tookpcr::Bool = false
    
    nra::Int16 = 0
    pcrprob::Float64 = 0.0
    test::Bool = false

    isofalse::Bool = false
    
    # proportion_contacts_workplace::Float64 = 0.0

    totaldaysiso::Int32 = 0  
    has_app::Bool = false
    contacts::Vector{Vector{Int16}} = [[0; 0]]
    ncontacts_day::Int8 = 0
    testedpos::Bool = false
    notified::Bool = false
    timetotest::Int64 = -1
    n_tests_perf::Int64 = 0
    time_since_testing::Int64 = 0

    n_neg_tests::Int64 = 0 # Taiye

end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :ontario
    calibration::Bool = false
    calibration2::Bool = false 
    start_several_inf::Bool = true
    modeltime::Int64 = 435
    initialinf::Int64 = 1
    fmild::Float64 = 0.5  ## percent of people practice self-isolation
    # Taiye: Could be useful later for keeping track of the population in isolation.

    start_testing::Int64 = 0
    test_for::Int64 = 0
    fsevere::Float64 = 1.0 #
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    
    # herd::Int8 = 0 #typemax(Int32) ~ millions
    # Taiye: We are not considering herd immunity at this stage.

    file_index::Int16 = 0
    
    app_coverage = 0.0
    track_days::Int8 = 3
    #the cap for coverage should be 90% for 65+; 95% for HCW; 80% for 50-64; 60% for 16-49; and then 50% for 12-15 (starting from June 1).
    #comor_comp::Float64 = 0.7 #prop comorbidade tomam
    
    #one waning rate for each efficacy? For each strain? I can change this structure based on that
    # reduce_days::Int64 = 0 # Taiye: We are not considering waning.
    # waning::Int64 = 1 # Taiye: We are not considering waning.

    ### after calibration, how much do we want to increase the contact rate... in this case, to reach 70%
    ### 0.5*0.95 = 0.475, so we want to multiply this by 1.473684211

    # hosp_red::Float64 = 3.1 # Taiye: We can add this if we decide to include hospitalizations.
    isolation_days::Int64 = 5
    ageintapp::Vector{Int64} = [10; 60]
    ##for testing

    test_ra::Int64 = 1 # Taiye (2025.06.24): 1 - PCR, 2 - Abbott_PanBio 3 - 	BD VERITO	4 - SOFIA
    # Taiye: I believe that PCR tests are the only ones being considered.

    time_until_testing::Int64 = 1
    #prop_working::Float64 = 0.65 #https://www.ontario.ca/document/ontario-employment-reports/april-june-2021#:~:text=Ontario's%20overall%20labour%20force%20participation,years%20and%20over%20at%2038.8%25.
    n_tests::Int64 = 2
    time_between_tests::Int64 = 0

    #n_neg_tests::Int64 = 0 # Taiye

    # Taiye (2025.06.09): Tentative contact change rates
    contact_change_rate::Float64 = 0.80
    contact_change_2::Float64 = 0.50

    # Taiye (2025.06.12): Attempting to correct 'ERROR: type ModelParameters has no field testing'
    testing::Bool = false

    # Taiye (2025.06.12): Defining initial_day_week.
    initial_day_week::Int64 = 1

    # Taiye (2025.06.24): asymp_red was not defined in matrices_code.jl.
    asymp_red::Float64 = 2 # Taiye (2025.06.24): tentative value
end


Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

include("matrices_code.jl")
include("matrices.jl")
## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
#const agebraks_vac = @SVector [0:0,1:4,5:14,15:24,25:44,45:64,65:74,75:100]
const BETAS = Array{Float64, 1}(undef, 0) ## to hold betas (whether fixed or seasonal), array will get resized

# const waning_factors = waning_factor() # Taiye: We are not considering waning.
# const waning_factors_rec = waning_factor() # Taiye: We are not considering waning.

export ModelParameters, HEALTH, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    # Taiye (2025.06.06): hmatrix, hh1, nra, npcr, nleft = main(ip,simnum)           
    hmatrix, hh1, nra = main(ip,simnum)  

    #Get the R0
    
    R01 = length(findall(k -> k.sickby in hh1,humans))/length(hh1)
    
    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    #ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    
    all1 = _collectdf(hmatrix)

    # Taiye (2025.06.12): We are not considering workplaces
    #spl = _splitstate(hmatrix, ags)
    #work = _collectdf(spl[1])
    
    age_groups = [0:14, 15:24, 25:34, 35:44, 45:54, 55:64, 65:999]
    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    ag7 = _collectdf(spl[7])
    insertcols!(all1, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum); 
    
    # Taiye (2025.06.12): We are not considering workplaces.
    #insertcols!(work, 1, :sim => simnum);
    

    pos = findall(y-> y in (11,22,33),hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,100)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    ### total days of isolation per age group, both working and general
    geniso_gr = map(y->findall(x-> x.age in y,humans),age_groups)
    #workiso_gr = map(y->findall(x-> x.age in y && x.workplace_idx > 0,humans),age_groups)

    giso = map(y-> sum([ii.totaldaysiso for ii in humans[y]]),geniso_gr)
    #wiso = map(y-> sum([ii.totaldaysiso for ii in humans[y]]),workiso_gr)

    # Taiye (2025.06.12):
    # return (a=all1, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7, work = work,vector_dead=vector_ded,nra=nra,npcr=npcr, R0 = R01, niso_t_p=niso_t_p, nleft=nleft,giso = giso, wiso = wiso)
    return (a=all1, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7,
    vector_dead=vector_ded,nra=nra,R0 = R01, giso = giso)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    
    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 

    #create herd immunity
    # herd_immu_dist_4(sim,1) # Taiye: We are not considering herd immunity.

    # split population in agegroups 
    grps = get_ag_dist()
    
    initial_dw::Int64 = 0

    nra::Vector{Int64} = zeros(Int64,p.modeltime)
    # (Taiye 2025.06.06): npcr::Vector{Int64} = zeros(Int64,p.modeltime)
    # (Taiye 2025.06.06): nleft::Vector{Int64} = zeros(Int64,p.modeltime)

    #insert one infected in the latent status in age group 4
    insert_infected(LAT, p.initialinf, 4)

    # h_init1 = findall(x->x.health_status  in (LAT,MILD,INF,PRE,ASYMP), humans) # Taiye: MILD is unnecessary.
    h_init1 = findall(x->x.health_status  in (LAT,INF,PRE,ASYMP), humans)
    
    ## save the preisolation isolation parameters
    #we need the workplaces to get the next days counts
    for x in humans
        get_nextday_counts(x)
    end
    
    # distributing app_coverage
    dist_app(humans, p)
     
    # start the time loop
    for st = 1:min((p.start_testing-1),p.modeltime)
        
        for x in humans
    #        if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
            if x.iso && !(x.health_status == DED) #&& !(x.health_status in (HOS,ICU,DED))
                x.totaldaysiso += 1

            end
        end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        sw = time_update() ###update the system
        
        # end of day
    end
    
    setfield!(p,:testing,true) 

    # Taiye: Attempt at implementation.
    #? Thomas: we cannot use this
    # for x in humans
    #     if x.time_since_testing < p.time_between_tests
    #         setfield!(p,:testing,false)
    #     else
            #setfield!(p,:testing,true)
    #     end
    # end


    # start the time loop
    for st = p.start_testing:min((p.start_testing+p.test_for-1),p.modeltime)
        
        for x in humans
            # if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
            if x.iso && !(x.health_status == DED)
                x.totaldaysiso += 1
            end
        end

        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        sw = time_update() ###update the system

        # Taiye (2025.06.12): sw might be a scalar
        # nra[st]+= sw[6]
        nra[st] += sw[length(sw)]

        # end of day
    end

    setfield!(p,:testing,false)

    for st = (p.start_testing+p.test_for):p.modeltime
        initial_dw = st+(p.initial_day_week-1)-7*Int(floor((st-1+(p.initial_day_week-1))/7))
        for x in humans
            # if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
            if x.iso && !(x.health_status == DED) #&& !(x.health_status in (HOS,ICU,DED))
                x.totaldaysiso += 1
                
            end
        end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        sw = time_update() ###update the system
        
        # end of day
    end
    
    
    # Taiye (2025.06.06): return hmatrix, h_init1, nra, npcr, nleft## return the model state as well as the age groups. 
    return hmatrix, h_init1, nra
end
export main

# Taiye (2025.06.09): We are not considering workplaces.
##### creating workplaces

# Taiye (2025.06.09)
#function work_size() #https://www150.statcan.gc.ca/t1/tbl1/en/cv.action?pid=3310039501
 #   breaks = [1:4,5:9,10:19,20:49,50:99,100:199,200:499,500:1000]
  
 #s = [748387, 233347, 152655, 99732, 32889, 14492, 7119, 2803]
    #s/=sum(s)
   # s = [0.5795052593106524,0.18068968828208243,0.11820672374061501,0.07722637956240554,0.025467236167207672,0.01122172113883589,0.00551251951334341,0.0021704722848576454]
   
   # Taiye (2025.06.09)
   #aux = Distributions.Categorical(@SVector [0.5795052593106524,0.18068968828208243,0.11820672374061501,0.07722637956240554,0.025467236167207672,0.01122172113883589,0.00551251951334341,0.0021704722848576454])

    #return aux,breaks
#end


function dist_app(humans, p)
    # Taiye (2025.06.12): pos = findall(x.age in p.ageintapp[1]:p.ageintapp[2], humans)
    pos = findall(x->x.age in p.ageintapp[1]:p.ageintapp[2], humans)
    pos = sample(pos, Int(round(p.app_coverage*p.popsize)))

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

    # reset the contact tracing data collection structure
    #= for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    =#    # resize and update the BETAS constant array
    #init_betas()

    # resize the human array to change population size
    resize!(humans, p.popsize)
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

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev) 
    
    # Taiye (2025.06.09): Replacing instances with iterate.
    #_names_inc = Symbol.(string.((Symbol.(iterate(HEALTH)[1:end - 1])), "_INC"))
    #_names_prev = Symbol.(string.((Symbol.(iterate(HEALTH)[1:end - 1])), "_PREV"))
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

function _get_incidence_and_prev(hmatrix)
    # Taiye (2025.06.09): Replacing instances with iterate.
    #cols = iterate(HEALTH)[1:end - 1]
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end
# Checkpoint

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findall(x-> r[x] == inth && r[x] != r[x-1],2:length(r))
        idx = idx .+ 1
        #idx = findfirst(x -> x == inth, r)
        if idx !== nothing
            for i in idx 
                timevec[i] += 1
            end
        end
    end
    return timevec
end

# Taiye: Since we are not considering herd immunity at this stage, this function might not be necessary.
#function herd_immu_dist_4(sim::Int64,strain::Int64)
 #   rng = MersenneTwister(200*sim)
  #  vec_n = zeros(Int32,6)
   # N::Int64 = 0
    #if p.herd == 5
     #   vec_n = [9; 148; 262;  68; 4; 9]
      #  N = 5

   # elseif p.herd == 10
    #    vec_n = [32; 279; 489; 143; 24; 33]

     #   N = 9

    #elseif p.herd == 20
     #   vec_n = [71; 531; 962; 302; 57; 77]

      #  N = 14
    #elseif p.herd == 30
     #   vec_n = [105; 757; 1448; 481; 87; 122]

      #  N = 16
   # elseif p.herd == 50
    #    vec_n = map(y->y*5,[32; 279; 489; 143; 24; 33])

     #   N = 16
    #elseif p.herd == 0
     #   vec_n = [0;0;0;0;0;0]
       
    #else
     #   vec_n = map(y->Int(round(y*p.herd/10)),[32; 279; 489; 143; 24; 33])
      #  N = 16
   # end

    #vprob::Vector{Float64} = vector_probs()
    #dprob = Distributions.Categorical(vprob)
   # for g = 1:6
       # pos = findall(y->y.ag_new == g && y.health == SUS,humans)
        #n_dist = min(length(pos),Int(floor(vec_n[g]*p.popsize/10000)))
        #pos2 = sample(rng,pos,n_dist,replace=false)
        #for i = pos2
         #   humans[i].strain = strain
          #  humans[i].swap = strain == 1 ? REC : REC2
           # humans[i].swap_status = REC
            #move_to_recovered(humans[i])
           # r = rand()
            #day = rand(dprob)
          #  humans[i].sickfrom = INF
           # humans[i].herd_im = true
  #      end
   # end
    #return N
#end

function _get_column_prevalence(hmatrix, hcol)
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
function get_province_ag(prov) 
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

# Taiye: We are not considering comorbidities at this stage.
#function comorbidity(ag::Int16)

 #   a = [4;19;49;64;79;999]
  #  g = findfirst(x->x>=ag,a)
   # prob = [0.05; 0.1; 0.28; 0.55; 0.74; 0.81]

    #com = rand() < prob[g] ? 1 : 0

#    return com    
#end
#export comorbidity


function initialize() 
    agedist = get_province_ag(p.prov)
    agebraksnew = [0:4,5:14,15:24,25:34,35:44,45:54,55:64,65:74,75:84,85:99]
    
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        agn = rand(agedist)
        x.age = rand(agebraksnew[agn]) 
        x.ag = findfirst(y-> x.age in y, agebraks)
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g
        x.exp = 999  ## susceptible people don't expire.
        
        #x.dur = sample_epi_durations() # sample epi periods   
      
      #  x.comorbidity = comorbidity(x.age) # Taiye: We are not considering comorbidities at this stage.
        x.time_since_testing = p.time_between_tests
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        x.contacts = repeat([[0]], p.track_days)
        
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
            # Taiye (2025.06.08): x.strain = strain
            x.dur = sample_epi_durations(x)

            # Taiye (2025.06.08): Removed if x.strain > 0
            #if x.strain > 0
            if health == PRE
                    x.swap = health
                    x.swap_status = PRE
                    x.daysinf = x.dur[1]+1
                    x.wentto = 1
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
            elseif health == LAT
                    x.swap = health
                    x.swap_status = LAT
                    x.daysinf = rand(1:x.dur[1])
                    move_to_latent(x)
                
                # Taiye: We are not currently considering MILD**.
                #elseif health == MILD
                 #   x.swap = health
                  #  x.swap_status = MILD
                   # x.wentto = 1
                    #x.daysinf = x.dur[2]+1
                    #move_to_mild(x)
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
            #else
             #   error("no strain insert inf")
            #end
            
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence

    nra::Int64 = 0
    
    if p.testing
        for x in humans
            if x.notified && x.n_tests_perf <= p.n_tests && x.timetotest == 0 && x.time_since_testing >= p.time_between_tests && x.n_neg_tests <= x.n_tests_perf # Taiye
                testing_infection(x, p.test_ra)
                
                x.time_since_testing = 0
                x.n_tests_perf += 1
                if x.n_tests_perf == p.n_tests
                    x.notified = false
                end
            end
        end
    end

    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
         
        x.daysinf += 1
        x.days_after_detection += 1 #we don't care about this untill the individual is detected
        x.daysisolation += 1
        x.timetotest -= 1

        x.time_since_testing += 1 # Taiye: We could measure this in days.
        

        if x.tis >= x.exp             
            @match Symbol(x.swap_status) begin
                :LAT  => begin 
                    move_to_latent(x); 
                end
                :PRE  => begin move_to_pre(x); end
                :ASYMP => begin move_to_asymp(x);  end
              # Taiye:  :MILD => begin nra+=move_to_mild(x); end
              
                :INF  => begin move_to_inf(x); end    
              
              # Taiye:   :HOS  => begin move_to_hospicu(x); end 
              # Taiye:   :ICU  => begin move_to_hospicu(x); end
                :REC  => begin move_to_recovered(x); end
                :DED  => begin move_to_dead(x); end
                _    => begin dump(x); error("swap expired, but no swap set."); end
            end
        end
        #if the individual recovers, we need to set they free. This loop must be here

        # if x.iso && x.daysisolation >= p.isolation_days && !(x.health_status in (HOS,ICU,DED))
        if x.iso && x.daysisolation >= p.isolation_days && !(x.health_status == DED) 
            _set_isolation(x,false,:null)
            
            x.n_tests_perf = 0 # Taiye
            x.n_neg_tests = 0 # Taiye

            # if x.testedpos # if the individual was tested and the days of isolation is finished, we can return the tested to false
            #     x.testedpos = false
            # end
            
        end
        # run covid-19 functions for other integrated dynamics. 
        #ct_dynamics(x)
        # get the meet counts for the next day 
        get_nextday_counts(x)
        
    end

   return nra
end
export time_update

#@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    if x.isovia == :null || via == :sev
        x.iso = iso 
        x.isovia = via
        x.daysisolation = 0
        x.days_after_detection = 0
    elseif !iso
        x.iso = iso 
        x.isovia = via
    end

end
function sample_epi_durations(y::Human)
    # when a person is sick, samples the 

    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3



    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    y.incubationp = latents
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    
    # Taiye (2025.06.09): Removing aux.
    asymps = max(Int.(ceil.(rand(asy_dist))),1)
    infs = max(Int.(ceil.(rand(inf_dist))),1)
    #asymps = max(Int.(ceil.(rand(asy_dist)))-aux,1)
    #infs = max(Int.(ceil.(rand(inf_dist)))-aux,1)

    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = x.swap
    x.health_status = x.swap_status
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
    age_thres = [4, 19, 49, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)

 
    if rand() < (symp_pcts[g])

       
        x.swap = PRE
        x.swap_status = PRE
        x.wentto = 1
        
    else
        x.swap = ASYMP
        x.swap_status = ASYMP
        x.wentto = 2
        
    end
    
    x.got_inf = true
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
    # Taiye: Removing if-statement
    #if x.strain == 1 
     #   θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    #elseif x.strain == 2 || x.strain == 3
     #   θ = (0.89, 0.78, 0.67, 0.48, 0.04)
    #else
     #   error("no strain in move to pre")
    #end  # percentage of sick individuals going to mild infection stage

    # Taiye (2025.06.12):
    θ = (0.95, 0.9, 0.85, 0.6, 0.2)
    
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period
    ##########

    
    if rand() < (1-θ[x.ag])
        x.swap = INF
        x.swap_status = INF
    
    # Taiye: This loop might be unnecessary.
    #else 
     #   x.swap = MILD
      #  x.swap_status = MILD
    end
    
end
export move_to_pre

function testing_infection(x::Human, teste)
    pp = _get_prob_test(x,teste)
    if rand() < pp
        x.testedpos = true
        _set_isolation(x, true, :test)

        # Taiye (2025.06.24): send_notifications(x)
        send_notification(x)

    else # Taiye: counting the number of negative tests performed.
          x.n_neg_tests += 1
    end

end

function send_notification(x::Human) # Taiye (2025.05.22): added an 's' to 'human'; Update: 'humans' -> 'Human'
    v = vcat(x.contacts...)
    for i in v
        if 1 <= i <= length(humans) && !humans[i].notified # Taiye: To avoid new notifications resetting times.
            humans[i].notified = true
            humans[i].timetotest = p.time_until_testing
        end
        #humans[i].time_since_testing = 0#p.time_between_tests # Taiye
    end

end

# Taiye: This function might be unnecessary.
#function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
 
 # Taiye:   
 #   x.health = x.swap 
  #  x.health_status = x.swap_status
   # x.tis = 0 
    #x.exp = x.dur[4]
 #   x.swap = REC
  #  x.swap_status = REC
    
    #x.swap = x.strain == 1 ? REC : REC2
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
   
   # Taiye:
   # nra = 0
    
    #if p.testing && !x.tested && x.has_app
     #   testing_infection(x, p.test_ra)
   # end

    #return nra
#end
#export move_to_mild


function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)

    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
 
    # Taiye: The following lines might be unnecessary.
    #comh = 0.98
    #h = x.comorbidity == 1 ? comh : 0.04 #0.376
    #c = x.comorbidity == 1 ? 0.396 : 0.25
    #time_to_hospital = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
   	
    x.health = x.swap
    x.health_status = x.swap_status
    x.swap = UNDEF
    
    x.tis = 0 
    
    #? Thomas:
    if p.testing && !x.testedpos && x.has_app
        #testing_infection(x, p.test_ra)
        x.notified = true

       # Taiye (2025.06.23): humans[i].timetotest = 1
        x.timetotest = 1
    end

    # This if-statement might be unnecessary.
   # if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
    #    x.exp = time_to_hospital
     #   if rand() < c
      #      x.swap = ICU
       #     x.swap_status = ICU
            
       # else
        #    x.swap = HOS
         #   x.swap_status = HOS
            
      #  end
       
   # else ## no hospital for this lucky (but severe) individual 
    if rand() < mh[gg]
            x.exp = x.dur[4] 
            x.swap = DED
            x.swap_status = DED
    else 
            x.exp = x.dur[4]  
            x.swap = REC
            x.swap_status = REC
    end
         
   # end
    ## before returning, check if swap is set 
    #x.swap == UNDEF && error("agent I -> ?")
end

# Taiye: This function might be unnecessary.
#function move_to_hospicu(x::Human)   
    #death prob taken from https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#Comorbidities
    # on May 31th, 2020
    #= age_thres = [24;34;44;54;64;74;84;999]
    g = findfirst(y-> y >= x.age,age_thres) =#
    #https://www.medrxiv.org/content/10.1101/2021.08.24.21262415v1
 #   aux = [0:4, 5:19, 20:44, 45:54, 55:64, 65:74, 75:84, 85:99]
   

  #      mh = [0.001, 0.001, 0.0015, 0.0065, 0.01, 0.02, 0.0735, 0.38]
   #     mc = [0.002,0.002,0.0022, 0.008, 0.022, 0.04, 0.08, 0.4]

    #gg = findfirst(y-> x.age in y,aux)

 #   psiH = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17))))
  #  psiC = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
   # muH = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15))))
    #muC = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

 #   swaphealth = x.swap_status 
  #  x.health = x.swap ## swap either to HOS or ICU
   # x.health_status = x.swap_status
    #x.swap = UNDEF
  #  x.tis = 0
   # _set_isolation(x, true, :hosp) # do not set the isovia property here.  

    #if swaphealth == HOS
         
     #   if rand() < mh[gg] ## person will die in the hospital 
      #      x.exp = muH 
       #     x.swap = DED
        #    x.swap_status = DED
           
      #  else 
       #     x.exp = psiH 
        #    x.swap = REC
         #   x.swap_status = REC
            
      #  end    
   # elseif swaphealth == ICU
              
    #    if rand() < mc[gg] ## person will die in the ICU 
     #       x.exp = muC
      #      x.swap = DED
       #     x.swap_status = DED
           
#        else 
 #           x.exp = psiC
  #          x.swap = REC
   #         x.swap_status = REC
            
    #    end
   # else
    #    error("error in hosp")
   # end
    
    ## before returning, check if swap is set 
 #   x.swap == UNDEF && error("agent H -> ?")    
#end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = h.swap
    h.health_status = h.swap_status
    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely

   # h.iso = true # a dead person is isolated
   # _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = h.swap
    h.health_status = h.swap_status
    
    h.recovered = true

    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely

    #h.iso = false ## a recovered person has ability to meet others
   
    #h.daysinf = 999

    h.got_inf = false 
    
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function _get_betavalue(xhealth) 
    #bf = p.β ## baseline PRE
    #length(BETAS) == 0 && return 0
    bf = p.β#BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.11

# Taiye: This elseif statement might not be necessary.
   # elseif xhealth == MILD || xhealth == MISO 
    #    bf = bf * 0.44

   # Taiye: elseif xhealth == INF || xhealth == IISO 
    elseif xhealth == INF #|| xhealth == IISO 
        bf = bf * 0.89
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    
    if !x.iso 
        aux =  p.contact_change_rate*p.contact_change_2
        cnt = rand(negative_binomials(ag,aux)) ##using the contact average for shelter-in
        
        x.nextday_meetcnt = cnt
    # elseif !(x.health_status  in (HOS,ICU,DED)) # Taiye
        cnt = rand(negative_binomials_shelter(ag,p.contact_change_2))  # expensive operation, try to optimize
      
    # Taiye (2025.06.09): nextday_meetcnt_w is only used here and could refer to workplaces, which would make it unnecessary.
      #  x.nextday_meetcnt_w = 0
      
        x.nextday_meetcnt = cnt
    #else (Taiye 2025.06.10)
    end
    
   # Taiye (2025.06.12): if x.health_status in (DED)
    if x.health_status == DED
        x.nextday_meetcnt = 0
    end

    # Taiye (2025.06.10)
    # if x.health_status in (HOS,ICU,DED) # Taiye
    #if x.health_status in (DED)
     #   x.nextday_meetcnt = 0
    #end
   
    for i in 2:p.track_days
        x.contacts[i] = deepcopy(x.contacts[i-1])
    end
    x.contacts[1] = repeat([0], max(x.nextday_meetcnt, 1))
    x.ncontacts_day = 0

    return cnt
end

function dyntrans(sys_time, grps,sim)
    #totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    #totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)
    pos = shuffle(1:length(humans))
    # go through every infectious person
    for x in humans[pos]        
    # Taiye:     if x.health_status in (PRE, ASYMP, MILD, MISO, INF, IISO)
       # if x.health_status in (PRE, ASYMP, INF)
            
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts
            #general population contact
            gpw = Int.(round.(cm[x.ag]*cnts))
            
            #cnts = number of contacts on that day

            
            perform_contacts(x,gpw,grps,xhealth)
                      
        
    end
    #return totalmet, totalinf
end
export dyntrans

function perform_contacts(x,gpw,grp_sample,xhealth)

    for (i, g) in enumerate(gpw) 
        meet = rand(grp_sample[i], g)   # sample the people from each group
        # go through each person
        for j in meet 
            y = humans[j]

        
            ycnt = y.nextday_meetcnt  
            
            y.nextday_meetcnt = y.nextday_meetcnt - min(1,ycnt) # remove a contact

            ycnt == 0 && continue
            
            if y.has_app && x.has_app
                x.ncontacts_day = x.ncontacts_day+1
                # Taiye (2025.06.27): Attempting to correct BoundsError: attempt to access 1-element Vector{Int16} at index [2]
                if length(x.contacts[1]) >= x.ncontacts_day
                    #x.contacts[1] = y.idx
                #else
                    x.contacts[1][x.ncontacts_day] = y.idx
                end

             # Taiye (2025.06.27): Attempting to correct BoundsError: attempt to access 1-element Vector{Int16} at index [2]
                y.ncontacts_day = y.ncontacts_day+1
                if length(y.contacts[1]) >= y.ncontacts_day
                    #y.contacts[1] = x.idx
                #else    
                    y.contacts[1][y.ncontacts_day] = x.idx
                end
            end
            #adj_beta = 0 # adjusted beta value by strain and vaccine efficacy
            if x.health_status in (PRE, ASYMP, INF) && y.health == SUS && y.swap == UNDEF 
                
                beta = _get_betavalue(x.health)
                
                if rand() < beta
                
                    y.exp = y.tis   ## force the move to latent in the next time step.
                    y.sickfrom = x.health ## stores the infector's status to the infectee's sickfrom
                    y.sickby = y.sickby < 0 ? x.idx : y.sickby
                # Taiye (2025.06.08): y.strain = x.strain       
                    y.swap = LAT
                    y.swap_status = LAT
                    y.daysinf = 0
                    y.dur = sample_epi_durations(y)
                    #y.swap = y.strain == 1 ? LAT : LAT2
                end  
            elseif y.health_status in (PRE, ASYMP, INF) && x.health == SUS && x.swap == UNDEF
                beta = _get_betavalue(y.health)
                if rand() < beta
                
                    x.exp = x.tis   ## force the move to latent in the next time step.
                    x.sickfrom = y.health ## stores the infector's status to the infectee's sickfrom
                    x.sickby = x.sickby < 0 ? y.idx : x.sickby
                # Taiye (2025.06.08): y.strain = x.strain       
                    x.swap = LAT
                    x.swap_status = LAT
                    x.daysinf = 0
                    x.dur = sample_epi_durations(x)
                    #y.swap = y.strain == 1 ? LAT : LAT2
                end  
            
            end
        
        end
    end  

end
function contact_matrix()
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

# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials(ag,mult) 
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
#const nbs = negative_binomials()
const cm = contact_matrix()
#export negative_binomials, contact_matrix, nbs, cm

export negative_binomials


function negative_binomials_shelter(ag,mult) 
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

#const vaccination_days = days_vac_f()
#const vac_rate_1 = vaccination_rate_1()
#const vac_rate_2 = vaccination_rate_2()
## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end