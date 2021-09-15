import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import time as tm

names = ['x', 'y1', 'y2', 'y3x', 'y3w', 'z1', 'z2', 'z3', 'w1', 'w2']

params = int(sys.argv[1])
run_number = sys.argv[2]

with open('paramgrid_appearance_barchart_twogenevaccine_mrates.txt', 'rb') as f: 
    param_file = pickle.load(f)

param_dict = param_file[params]

def calculate_reaction_rates(population, s, beta1, beta2, beta3, mu, c, a, d, eps):
    
    # 0:x; 1:y1; 2:y2; 3:y3x; 4:y3w; 5:z1; 6:z2; 7:z3; 8:w1; 9:w2;
    
    if (population[0]+population[5]+population[6]+population[7])>0:

        rates = np.array([beta1*s*population[0]*population[1],  # infection of susceptible with WT
                          beta2*s*population[0]*population[2],  # infection of susceptible with MT 
                          beta3*s*population[0]*population[3],  # infection of susceptible with EM (susceptible)
                          beta3*s*population[0]*population[4],  # infection of susceptible with EM (vaccinated)
                          beta3*s*population[5]*population[3],  # infection of recovered from WT with EM (susceptible)
                          beta3*s*population[5]*population[4],  # infection of recovered from WT with EM (vaccinated)
                          beta3*s*population[6]*population[3],  # infection of recovered from MT with EM (susceptible)
                          beta3*s*population[6]*population[4],  # infection of recovered from MT with EM (vaccinated)
                          beta3*s*population[8]*population[3],  # infection of vaccinated from MT with EM (susceptible)
                          beta3*s*population[8]*population[4],  # infection of vaccinated from MT with EM (vaccinated) 
                          beta1*s*population[0]*population[1]*mu,# mutation to MT during infection with WT
                          beta2*s*population[0]*population[2]*mu,# mutation to EM during infection with MT
                          a*population[1],                      # recovered from WT
                          a*population[2],                      # recovered from MT
                          a*population[3],                      # recovered from EM, was susceptible
                          a*population[4],                      # recovered from EM, was vaccinated
                          d*population[1],                      # death
                          d*population[2],                      # death
                          d*population[3],                      # death
                          d*population[4],                      # death
                          (c*population[0])/(population[0]+population[5]+population[6]+population[7]), # vaccination of susceptible
                          (c*population[5])/(population[0]+population[5]+population[6]+population[7]), # vaccination of WT recovered
                          (c*population[6])/(population[0]+population[5]+population[6]+population[7]), # vaccination of MT recovered
                          (c*population[7])/(population[0]+population[5]+population[6]+population[7])])# vaccination of EM recovered

    else:
        rates = np.array([beta1*s*population[0]*population[1],  # infection of susceptible with WT
                          beta2*s*population[0]*population[2],  # infection of susceptible with MT 
                          beta3*s*population[0]*population[3],  # infection of susceptible with EM (susceptible)
                          beta3*s*population[0]*population[4],  # infection of susceptible with EM (vaccinated)
                          beta3*s*population[5]*population[3],  # infection of recovered from WT with EM (susceptible)
                          beta3*s*population[5]*population[4],  # infection of recovered from WT with EM (vaccinated)
                          beta3*s*population[6]*population[3],  # infection of recovered from MT with EM (susceptible)
                          beta3*s*population[6]*population[4],  # infection of recovered from MT with EM (vaccinated)
                          beta3*s*population[8]*population[3],  # infection of recovered from MT with EM (susceptible)
                          beta3*s*population[8]*population[4],  # infection of recovered from MT with EM (vaccinated) 
                          beta1*s*population[0]*population[1]*mu,# mutation to MT during infection with WT
                          beta2*s*population[0]*population[2]*mu,# mutation to EM during infection with MT
                          a*population[1],                      # recovered from WT
                          a*population[2],                      # recovered from MT
                          a*population[3],                      # recovered from EM, was susceptible
                          a*population[4],                      # recovered from EM, was vaccinated
                          d*population[1],                      # death
                          d*population[2],                      # death
                          d*population[3],                      # death
                          d*population[4],                      # death
                          0., # vaccination of susceptible
                          0., # vaccination of WT recovered
                          0., # vaccination of MT recovered
                          0.])# vaccination of EM recovered
    return(rates)
    
# 0:x; 1:y1; 2:y2; 3:y3x; 4:y3w; 5:z1; 6:z2; 7:z3; 8:w1; 9:w2;
                  # infections of susceptible
update_vectors = [[-1,1,0,0,0,0,0,0,0,0],
                  [-1,0,1,0,0,0,0,0,0,0],
                  [-1,0,0,1,0,0,0,0,0,0],
                  [-1,0,0,1,0,0,0,0,0,0],
                  
                  # infections of recovered from WT
                  [0,0,0,1,0,-1,0,0,0,0],
                  [0,0,0,1,0,-1,0,0,0,0],
                  
                  # infections of recovered from MT
                  [0,0,0,1,0,0,-1,0,0,0],
                  [0,0,0,1,0,0,-1,0,0,0],
                  
                  # infections of vaccinated
                  [0,0,0,0,1,0,0,0,-1,0],
                  [0,0,0,0,1,0,0,0,-1,0],
                                    
                  # mutations
                  [-1,0,1,0,0,0,0,0,0,0],
                  [-1,0,0,1,0,0,0,0,0,0],
                  
                  # recoveries
                  [0,-1,0,0,0,1,0,0,0,0],
                  [0,0,-1,0,0,0,1,0,0,0],
                  [0,0,0,-1,0,0,0,1,0,0],
                  [0,0,0,0,-1,0,0,0,0,1],
                  
                  # deaths
                  [0,-1,0,0,0,0,0,0,0,0],
                  [0,0,-1,0,0,0,0,0,0,0],
                  [0,0,0,-1,0,0,0,0,0,0],
                  [0,0,0,0,-1,0,0,0,0,0],

                  # vaccinations
                  [-1,0,0,0,0,0,0,0,1,0],
                  [0,0,0,0,0,-1,0,0,1,0],
                  [0,0,0,0,0,0,-1,0,1,0],
                  [0,0,0,0,0,0,0,-1,0,1]]
                  
a = param_dict['a']
d = param_dict['d']*a
N = param_dict['N']
s_threshold = 0.1
i_threshold = param_dict['i_threshold']
s = 1/3
desired_R =  param_dict['R']
beta1 =(a*desired_R)/(N)
beta2 =(a*desired_R)/(N)
beta3 = (a*desired_R)*(param_dict['delta']/N)
mu = param_dict['mu']

c = param_dict['c']
eps = param_dict['q']

period = -1
period_length = 1.

population = [N-0-i_threshold, i_threshold, 0, 0, 0, 0, 0, 0,0,0]
pop_evol = [np.array(population)]

reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, beta3, mu, c, a, d, eps)

appearance_VR = []
timing_VR = []

for nrun in range(10):

    start_time = tm.time()

    population = [N-i_threshold, i_threshold, 0, 0, 0, 0, 0, 0,0,0]
    pop_evol = [np.array(population)]

    time = 0.
    times = [time]
    counter = 0

    reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, beta3, mu, c, a, d, eps)

    vaccination_period = 0.

    VR_infected = False

    total_infections = [0.]
    infections_per_day = [0.]
    s = 1/3
    sevol = [s]
   

    infection_counter = 0

    while not all([population[cat] <= 0 for cat in [1,2,3,4]]):
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

        if population[4]+population[3]>population[1] and not VR_infected:
            timing_VR.append(time)
            print(time, population[4]+population[3], population[1])
            VR_infected = True

        if index in [0,1,2,3,4,5,6,7,8,9,10,11]:
            infection_counter += 1

        if time // period_length != period:
            new_infections = infection_counter-total_infections[-1]
            infections_per_day.append(new_infections)

            s_threshold = np.random.uniform(0,0.1)
            
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

            sevol.append(s)
            total_infections.append(infection_counter)

            increment = time // period_length - period
            period = time // period_length

        # Update Reactions Rates
        reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, beta3, mu, c, a, d, eps)

    appearance_VR.append(VR_infected)

    end_time = tm.time() - start_time

    print(end_time)
    
output = 'COVID_vaccine_resistance_review/two_gene_vaccine/'
#output_argv = 'p'+str(params)
output_argv = 'p'+str(params)+'_r'+str(run_number)

with open(output + 'results_gillespie_paramgrid_' + output_argv + '.txt', 'wb') as f:
    pickle.dump({'params':params, 'times':timing_VR, 'presence_VR':appearance_VR}, f)