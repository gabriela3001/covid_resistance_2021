import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import time as tm

period = -1
period_length = 1

params = int(sys.argv[1])
run_number = int(sys.argv[2])

with open('paramgrid_appearance_barchart_linearrealdata.txt', 'rb') as f: 
    param_file = pickle.load(f)
    
param_dict = param_file[params]

with open('vaccination_data.txt', 'rb') as f:
    vaccination_data = pickle.load(f)
    
vaccination_vector = vaccination_data[param_dict['countries']]
    
a = 0.25
N = param_dict['N']
s_threshold = param_dict['s_threshold']*period_length
d = param_dict['d']*a
N = param_dict['N']
i_threshold = param_dict['i_threshold']
s = 1/3
beta1 = 3*a/N
beta2 = 3*a/N
mu = param_dict['mu']
q = param_dict['q']
hesitancy = param_dict['hesitancy']

population = [N-i_threshold, i_threshold, 0, 0, 0, 0, 0, 0]
pop_evol = [np.array(population)]

def calculate_reaction_rates(population, s, beta1, beta2, mu, c, a, d, q):
    
    # 0:x; 1:y1; 2:y2x; 3:y2w; 4:z1; 5:z2; 6:w1; 7:w2
    
    if (population[0]+population[4]+population[5])>0:

        rates = np.array([beta1*s*population[0]*population[1],  # infection of susceptible with WT
                          beta2*s*population[0]*population[2],  # infection of susceptible with VR (susceptible)
                          q*beta2*s*population[4]*population[2],  # infection of WT recovered with VR (susceptible)
                          q*beta2*s*population[6]*population[2],  # infection of vaccinated with VR (susceptible)
                          beta2*s*population[0]*population[3],  # infection of susceptible with VR vaccinated
                          q*beta2*s*population[4]*population[3],  # infection of WT recovered with VR vaccinated
                          q*beta2*s*population[6]*population[3],  # infection of vaccinated with VR vaccinated
                          beta1*s*population[0]*population[1]*mu,# mutation occurring in an individual
                          a*population[1],                      # recovered from WT
                          a*population[2],                      # recovered from VR, was susceptible
                          a*population[3],                      # recovered from VR, was vaccinated
                          d*population[1],                      # death
                          d*population[2],                      # death
                          d*population[3],                      # death
                          (c*population[0])/(population[0]+population[4]+population[5]), # vaccination of susceptible
                          (c*population[4])/(population[0]+population[4]+population[5]), # vaccination of WT recovered
                          (c*population[5])/(population[0]+population[4]+population[5])])# vaccination of WT and VR recovered

    else:
        rates = np.array([beta1*s*population[0]*population[1],  # infection of susceptible with WT
                  beta2*s*population[0]*population[2],  # infection of susceptible with VR (susceptible)
                  q*beta2*s*population[4]*population[2],  # infection of WT recovered with VR (susceptible)
                  q*beta2*s*population[6]*population[2],  # infection of vaccinated with VR (susceptible)
                  beta2*s*population[0]*population[3],  # infection of susceptible with VR vaccinated
                  q*beta2*s*population[4]*population[3],  # infection of WT recovered with VR vaccinated
                  q*beta2*s*population[6]*population[3],  # infection of vaccinated with VR vaccinated
                  beta1*s*population[0]*population[1]*mu,# mutation occurring in an individual
                  a*population[1],                      # recovered from WT
                  a*population[2],                      # recovered from VR, was susceptible
                  a*population[3],                      # recovered from VR, was vaccinated
                  d*population[1],                      # death
                  d*population[2],                      # death
                  d*population[3],                      # death
                  0., # vaccination of susceptible
                  0., # vaccination of WT recovered
                  0.])# vaccination of WT and VR recovered



        
    return(rates)
    
reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, 0, a, d, q)

                  # infections by WT
update_vectors = [[-1,1,0,0,0,0,0,0],
                  
                  # infections by y2x
                  [-1,0,1,0,0,0,0,0],  
                  [0,0,1,0,-1,0,0,0], 
                  [0,0,0,1,0,0,-1,0],
                  
                  # infections by y2w
                  [-1,0,1,0,0,0,0,0], 
                  [0,0,1,0,-1,0,0,0], 
                  [0,0,0,1,0,0,-1,0],
                  
                  # mutation
                  [-1,0,1,0,0,0,0,0],
                  
                  # recoveries
                  [0,-1,0,0,1,0,0,0],
                  [0,0,-1,0,0,1,0,0], 
                  [0,0,0,-1,0,0,0,1],
                  
                  # deaths
                  [0,-1,0,0,0,0,0,0], 
                  [0,0,-1,0,0,0,0,0], 
                  [0,0,0,-1,0,0,0,0],
                  
                  # vaccination
                  [-1,0,0,0,0,0,1,0], 
                  [0,0,0,0,-1,0,1,0], 
                  [0,0,0,0,0,-1,0,1]]
                  
appearance_VR = []
timing_VR = []
deaths = []

for nrun in range(10):

    start_time = tm.time()

    population = [N-int(i_threshold), int(i_threshold), 0, 0, 0, 0, 0, 0]
    pop_evol = [np.array(population)]

    time = 0.
    times = [time]
    counter = 0
    

    c = vaccination_vector[0]
    reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, c, a, d, q)

    vaccination_period = 0.

    VR_infected = False
    year_passed = False

    total_infections = [0.]
    infections_per_day = [0.]
    s = 1/3
    sevol = [s]
   

    infection_counter = 0

    while not all([population[cat] <= 0 for cat in [1,2,3]]):

        counter += 1

        # Generate random numbers
        r1 = np.random.uniform(0,1)
        r2 = np.random.uniform(0,1)

        # Calculate alpha
        cumsum = reaction_rates.cumsum()
        alpha = cumsum[-1]

        # Reaction Time        
        tau = (1.0/alpha)*np.log(float(1.0/r1))
        time += tau
        index = np.searchsorted(cumsum,r2*alpha)
        population = np.array(population) + np.array(update_vectors[index])

        if counter % 1000 == 0:
            pop_evol.append(population)
            times.append(time)
            sevol.append(s)

        if population[2]+population[3]>population[1] and not VR_infected:
            timing_VR.append(time)
            print(time, population[2]+population[3], population[1])
            VR_infected = True

        if index in [0,1,2,3,4,5,6]:
            infection_counter += 1

        if time // period_length != period:
            new_infections = (infection_counter-total_infections[-1])/period_length
            infections_per_day.append(new_infections)

            #s_threshold = np.random.uniform(0.,0.1)
            
            # tightening lockdown 
            if s  > s_threshold and new_infections > i_threshold:
                s -= s_threshold

            # loosening lockdown
            if s < 1. and new_infections < i_threshold:
                s += s_threshold
                
            if s < 0:
                s = 0.05
            elif s > 1:
                s = 1
            
            total_infections.append(infection_counter)

            increment = time // period_length - period
            period = time // period_length

        # Update Reactions Rates
        if time > len(vaccination_vector):
            break
            c = np.min(np.nonzero(vaccination_vector))
            if 1-(population[6]+population[7])/N<hesitancy:
                reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, 0, a, d, q)
            else:
                reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, c, a, d, q)
        else:
            c = vaccination_vector[int(np.floor(time))]
            if 1-(population[6]+population[7])/N<hesitancy:
                reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, 0, a, d, q)
            else:
                reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, c, a, d, q)     
                
        if time > 365 and not year_passed:
          year_passed = True
          deaths.append(N-np.sum(population))       
            
    appearance_VR.append(VR_infected)
    
    if year_passed == False:
      year_passed = True
      deaths.append(N-np.sum(population))        

    end_time = tm.time() - start_time

    print(end_time) 
    
output = COVID_vaccine_resistance_review/bar_chart/vaccine_real_data_plugin/'
output_argv = 'p'+str(params)+'_r'+str(run_number)

with open(output + 'results_gillespie_paramgrid_' + output_argv + '.txt', 'wb') as f:
    pickle.dump({'params':params, 'times':timing_VR, 'presence_VR':appearance_VR, 'death_rates':deaths}, f)  

    