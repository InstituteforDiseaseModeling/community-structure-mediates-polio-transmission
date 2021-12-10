import os
import sys
from collections import defaultdict
from Demographics import Demographics, write_dill #this is necessary!
import dill
import config
import json
import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "..", "lhs_search", "mass_action") #dump files in lhs_sampling folder
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

def search_files(directory='.', extension='.pkl'):
    found_files  = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                found_files.append(os.path.join(dirpath, name))
       
    return found_files

beta_global = int(sys.argv[1]) # the value of beta_global to be examined; allows for paramter sweep
iteration = int(sys.argv[2])
Matlab_village_files = search_files(r'/n/home04/weswong/multiscale_polio/MultiscaleModeling/PopSim/Demographics/Matlab_Villages', '.pkl')

fecal_oral_dose = 1.25e-06
c = {"fecal_oral_dose": fecal_oral_dose,
     "hh_transmission": 0,
     "beta_hh": 0,
     "bari_transmission": 0,
     "beta_bari": 0,
     "village_transmission": 0,
     "beta_village": 0,
     "age_assortivity": 0,
     "beta_inter_village": 0,

     "global_transmission": 1,
     "beta_global": beta_global}

config.params.overwrite(c)

prevalences = defaultdict(list)
r0, transmissions = [], []
sim_cohort_individuals = defaultdict(list)
sim_cohort_incidences = defaultdict(list)
transmission_levels = defaultdict(list)

for _ in range(10):#(30):
    random_idx = np.random.randint(0, len(Matlab_village_files) -1 )
    village_object_file = Matlab_village_files[random_idx]
    cohort_individuals = defaultdict(lambda: 0)
    
    
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)

    D.Taniuchi_study(config.params)
    for key in D.sim_prevalences:
        prevalences[key].append(D.sim_prevalences[key])
        
    for key in D.cohort_incidences:
        sim_cohort_incidences[key].append(D.cohort_incidences[key])
    
    for V in D.villages:
        for key in V.infants:
            cohort_individuals['infants_' + key] += len([i.id for i in V.infants[key]])
        for key in V.hh_contacts:
            cohort_individuals['hh_' + key] += len([i.id for i in V.hh_contacts[key]])

    for key in D.n_events_per_level:
        transmission_levels[key].append(D.n_events_per_level[key])

    for key in cohort_individuals:
        sim_cohort_individuals[key].append(cohort_individuals[key])

fout = 'mass_action_{x}.{iteration}'.format(x=beta_global, iteration = iteration)
with open(os.path.join(output_dir, '') + '{x}.json'.format(x=fout), 'w') as fp:
    json.dump([prevalences, sim_cohort_individuals, sim_cohort_incidences, transmission_levels], fp)
