

ip = ModelParameters(app_coverage = 0.6, testing = true,
not_swit = true, test_sens = 1)

sim = 1

# checking coverage
n = length(findall(x -> x.age >= 18 && x.age <=65, humans))
length(findall(x -> x.has_app, humans))/n


#! Main code

Random.seed!(sim*726)

# Taiye (2025.11.18):
rng = MersenneTwister(246*sim)

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
dist_app(humans, p, p.num_sims,rng)





#! Testing if a symptomatic isolates
#? Let's create an individual who will be set to symptomatic
# they must have the app
x = humans[findall(x -> x.has_app, humans)[1]]

#? Now, let's create the sample epi dur
 x.dur = sample_epi_durations(x)

 #? change the swap to INF
 x.swap = x.swap_status = INF
 #? And change the health to INF
 move_to_inf(x)
x

#? Now, let's check how long the individual stays isolated

st = 1

for x in humans
    # if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
    if x.iso && !(x.health_status == DED)
        x.totaldaysiso += 1
    end
end

_get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
dyntrans(st, grps,sim)
sw = time_update(st,sim,rng) ###update the system

#* That's working fine

#! Testing if a symptomatic isolates and send notification

x = humans[findall(x -> x.has_app, humans)[1]]

#? Now, let's create the sample epi dur
 x.dur = sample_epi_durations(x)

 #? change the swap to INF
 x.swap = x.swap_status = INF
 #? And change the health to INF
 move_to_inf(x)
x

x.contacts[1] = findall(x -> x.has_app, humans)[2:4]

#? Let's force one of them to become infected
y = humans[x.contacts[1][1]]

y.exp = y.tis   ## force the move to latent in the next time step.
y.sickfrom = x.health ## stores the infector's status to the infectee's sickfrom
y.sickby = y.sickby < 0 ? x.idx : y.sickby

y.swap = LAT
y.swap_status = LAT
y.daysinf = 0
y.dur = sample_epi_durations(x)

move_to_latent(y)

send_notification(x,p.not_swit,st,sim,rng)

#! Testing if a notified person tests

for x in humans
    # if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
    if x.iso && !(x.health_status == DED)
        x.totaldaysiso += 1
    end
end

_get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
dyntrans(st, grps,sim)
sw = time_update(st,sim,rng) ###update the system

#* They are testing
[humans[i].pp for i in [6, 9, 13]]

humans[6].testedpos
humans[6].iso
humans[6].isolation_days
humans[6].exp


#! Testing if a the notified person who tested positive also isolates
 _get_prob_test(y,teste,p.test_sens) #p.test_sens == 1 set prob to 1

#! Testing if the notified person who is isolated reset isolation if symptomatic
#? Forcing symptoms


humans[6].swap = PRE
humans[6].swap_status = PRE



st = 1

for x in humans
    # if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
    if x.iso && !(x.health_status == DED)
        x.totaldaysiso += 1
    end
end

_get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
dyntrans(st, grps,sim)
sw = time_update(st,sim,rng) ###update the system

y.health

humans[6].isolation_days
humans[6].daysisolation