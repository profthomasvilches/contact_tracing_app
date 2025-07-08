# Taiye (2025.06.23)

# Taiye (2025.05.30): run_param_scen_cal has 16 inputs here but 21 in simulations_cluster.jl.

#for the paper

#run_param_scen_cal(b,province,h_i,ic1,strains,index,scen,tra,eb,wpt,dayst,rc,dc,mt,vac,nsims)

#0.35 for omicron generates a R0 of 0.84
#b::Float64,province::String="ontario",ic1::Int64=1,index::Int64 = 0,scen::Int64 = 0,test_time::Int64=0,test_dur::Int64=0,mt::Int64=300,nsims::Int64=500)

a_c = 0.05 # Taiye (2025.06.23): app_coverage
cv.ModelParameters(popsize=10)

# Taiye (2025.07.01): Uncomment after rectifying coverage.
#for i = 0:20
    # run_param_scen_cal(0.1,"ontario",1,0,10,50,200,2,1000,a_c*i)
 #   run_param_scen_cal(1.0,"ontario",1,0,1,200,200,500,1000,a_c*i)
#end

run_param_scen_cal(0.1,"ontario",1,0,1,200,300,500,1000,0.0)
run_param_scen_cal(0.1,"ontario",1,0,1,200,300,500,1000,0.5)
run_param_scen_cal(0.1,"ontario",1,0,1,200,300,500,1000,0.8)
run_param_scen_cal(0.1,"ontario",1,0,1,200,300,500,1000,1.0)

# run_param_scen_cal(0.1,"ontario",1,0,10,50,200,2,1000)