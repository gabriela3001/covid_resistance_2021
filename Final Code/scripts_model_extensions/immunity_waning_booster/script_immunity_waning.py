import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import time as tm

params = int(sys.argv[1])
run_number = int(sys.argv[2])

with open('paramgrid_appearance_barchart_immunitywaning.txt', 'rb') as f: 
    param_file = pickle.load(f)
    
param_dict = param_file[params]
print(param_dict)
s_threshold = param_dict['s_threshold']
i_threshold = param_dict['i_threshold']
s = 1/3
a = 0.25
N = param_dict['N']
beta1 = 3*a/N
beta2 = param_dict['delta']*3*a/N
mu = param_dict['mu']
d = param_dict['d']*a

c1 = param_dict['c']
c2 = param_dict['c']
q = param_dict['q']
theta = 180

period = -1
period_length = 1.

def calculate_reaction_rates(population, s, beta1, beta2, mu, q, c1, c2, a, d, theta):
    
    # 0:x; 1:y1; 2:y2x; 3:y2w; 4:z1; 5:z2; 6:w1; 7:w2; 8:v1; 9:v2; 10:v3; 11:v4; 12:v5: 13:v6; 14:v7; 15:v8
    
    x, y1, y2A, y2B, z1, z2, w1, w2, vx, vy1, vy2A, vy2B, vz1, vz2, vw1, vw2 = population
    
    total_WT_infected = y1+vy1
    total_MT_infected = y2A+y2B+vy2A+vy2B
    
    total_to_vaccinate = x+z1+z2
    total_to_boost = vx+vz1+vz2
    
    all_are_vaccinated = int(total_to_vaccinate==0)
    all_are_boosted = int(total_to_boost==0)
    
    rates = np.array([beta1*s*x*total_WT_infected, # susceptible gets infected by WT
                      beta2*s*x*total_MT_infected, # susceptible gets infected by MT
                      q*beta2*s*z1*total_MT_infected, # recovered from WT gets infected by MT
                      q*beta2*s*w1*total_MT_infected, # vaccinated to WT gets infected by MT
                      
                      beta1*s*vx*total_WT_infected, # susceptible gets infected by WT
                      beta2*s*vx*total_MT_infected, # susceptible gets infected by MT
                      q*beta2*s*vz1*total_MT_infected, # recovered from WT gets infected by MT
                      q*beta2*s*vw1*total_MT_infected, # vaccinated to WT gets infected by MT
                      
                      beta1*s*x*total_WT_infected*mu, # mutation occurs in susceptible getting infected
                      beta1*s*vx*total_WT_infected*mu, # mutation occurs in susceptible getting infected
                      
                      a*y1,  # recovered from WT
                      a*y2A, # recovered from MT, was susceptible
                      a*y2B, # recovered from MT, was vaccinated
                      a*vy1,  # recovered from WT
                      a*vy2A, # recovered from MT, was susceptible
                      a*vy2B, # recovered from MT, was vaccinated
                      
                      d*y1,  # died from WT
                      d*y2A, # died from MT, was susceptible
                      d*y2B, # died from MT, was vaccinated
                      d*vy1,  # died from WT
                      d*vy2A, # died from MT, was susceptible
                      d*vy2B, # died from MT, was vaccinated
                      
                      [(c1*x)/(total_to_vaccinate+[0,1e-50][all_are_vaccinated]),0][all_are_vaccinated],  # first vaccination of susceptible
                      [(c1*z1)/(total_to_vaccinate+[0,1e-50][all_are_vaccinated]),0][all_are_vaccinated], # first vaccination of WT recovered
                      [(c1*z2)/(total_to_vaccinate+[0,1e-50][all_are_vaccinated]),0][all_are_vaccinated], # first vaccination of MT recovered
                      
                      [(c2*vx)/(total_to_boost+[0,1e-50][all_are_boosted]),0][all_are_boosted],  # first vaccination of susceptible
                      [(c2*vz1)/(total_to_boost+[0,1e-50][all_are_boosted]),0][all_are_boosted], # first vaccination of WT recovered
                      [(c2*vz2)/(total_to_boost+[0,1e-50][all_are_boosted]),0][all_are_boosted], # first vaccination of MT recovered     
    
                      z1/theta,
                      w1/theta
                      
                     ])  
    return(rates)
    
                # infections by WT
update_vectors = [[-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                  [-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                  [0,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0],
                  [0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0],
                  [0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0],
                  
                  # mutations
                  [-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0],
                  
                  # recovery
                  [0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1],
                  
                  # death
                  [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],

                  # vaccination first doses
                  [-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1],
                  
                  # waning
                  [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
                  ]

appearance_VR = []
timing_VR = []
deaths = []

for nrun in range(10):

    start_time = tm.time()

    population =  [N-i_threshold, i_threshold, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    pop_evol = [np.array(population)]

    time = 0.
    times = [time]
    counter = 0

    reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, q, c1, c2, a, d, theta)

    vaccination_period = 0.

    VR_infected = False
    year_passed = False

    total_infections = [0.]
    infections_per_day = [0.]
    s = 1/3
    sevol = [s]

    infection_counter = 0

    while not all([population[cat] <= 0 for cat in [1,2,3,9,10,11]]):

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

        if population[2]+population[3]>population[1] and not VR_infected:
            timing_VR.append(time)
            print(time, population[2]+population[3], population[1])
            VR_infected = True

        if index in [0,1,2,3,4,5,6,7]:
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
        if time < theta:
            reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, q, c1, 0, a, d, theta)
        else:
            reaction_rates = calculate_reaction_rates(population, s, beta1, beta2, mu, q, c1, c2, a, d, theta)
            
        if time > 365*3:
          break
          
        if time > 365 and not year_passed:
          year_passed = True
          deaths.append(N-np.sum(population))

    appearance_VR.append(VR_infected)
    
    if year_passed == False:
      deaths.append(N-np.sum(population))
    end_time = tm.time() - start_time

    print(end_time)

output = 'COVID_vaccine_resistance_review/bar_chart/immunity_waning_booster/'
output_argv = 'p'+str(params)+'_r'+str(run_number)

with open(output + 'results_gillespie_paramgrid_' + output_argv + '.txt', 'wb') as f:
    pickle.dump({'params':params, 'times':timing_VR, 'presence_VR':appearance_VR, 'death_rates':deaths}, f)  

