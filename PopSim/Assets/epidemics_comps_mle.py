import os
import sys
from collections import defaultdict
from Demographics import Demographics, write_dill
import dill
import config
import json
import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))


def search_files(directory='.', extension='.pkl'):
    found_files = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                found_files.append(os.path.join(dirpath, name))

    return found_files


type = sys.argv[1] #should be either multi or ma (mass action); will automatically run the correct parameters
idx = str(sys.argv[2]) #an idx appended to the end of the file to allow multiple runs to occur in parallel

#this should be a directory path to where the demographic files that we will run the taniuchi study on
Matlab_village_files = search_files(r'/n/home04/weswong/multiscale_polio/MultiscaleModeling/PopSim/Demographics/Matlab_Villages', '.pkl')
fecal_oral_dose = 1.25e-06

#multiscale
if type == 'multi':
    transmission_configs = {"fecal_oral_dose" : 1.25e-06,
                        "hh_transmission" : 1,
                    "beta_hh" : 1,
                    "bari_transmission" : 1,
                    "beta_bari" : 15,
                    "village_transmission" : 1,
                    "beta_village" : 4,
                    "age_assortivity": 0,
                    "beta_inter_village":2,
                    "global_transmission" : 0}

else:
    transmission_configs = {"fecal_oral_dose" : 1.25e-06,
                    "hh_transmission" : 0,
                    "bari_transmission" : 0,
                    "village_transmission" : 0,
                    "global_transmission" : 1,
                    "beta_global": 2,#19,
                    "age_assortivity": 0}


print(transmission_configs)
config.params.overwrite(transmission_configs)   
prevalences = defaultdict(list)
r0, transmissions = [], []
sim_cohort_individuals = defaultdict(list)
sim_cohort_incidences = defaultdict(list)
transmission_levels = defaultdict(list)

for _ in range(10):
    random_idx = np.random.randint(0, len(Matlab_village_files) -1 )
    village_object_file = Matlab_village_files[random_idx]
    cohort_individuals = defaultdict(list)
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)

    D.Taniuchi_study(config.params)
    for key in D.sim_prevalences:
        prevalences[key].append(D.sim_prevalences[key])
        
    for key in D.cohort_incidences:
        sim_cohort_incidences[key].append(D.cohort_incidences[key])
    
    for V in D.villages:
        for key in V.infants:
            cohort_individuals['infants_' + key] += [i.id for i in V.infants[key]]
        for key in V.hh_contacts:
            cohort_individuals['hh_' + key] += [i.id for i in V.hh_contacts[key]]
    
    for key in D.n_events_per_level:
        transmission_levels[key].append(D.n_events_per_level[key])
    
    #transmissions.append(D.transmission_tracker)
    
    for key in cohort_individuals:
        sim_cohort_individuals[key].append(cohort_individuals[key])
        



if type == 'multi':
    fout = '{f}_{hh}_{bari}_{village}_{inter}_{idx}'.format(f=transmission_configs["fecal_oral_dose"],
                                                            hh=transmission_configs["beta_hh"],
                                                            bari=transmission_configs["beta_bari"],
                                                            village=transmission_configs["beta_village"],
                                                            inter=transmission_configs["beta_inter_village"], idx=idx)

    output_dir = os.path.join(dir_path, "..", "Taniuchi_epidemics", "multi", '')  # dump files in lhs_sampling folder
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(output_dir + 'epidemic_' + fout + '.json', 'w') as fp:
       json.dump([prevalences, sim_cohort_individuals, sim_cohort_incidences, transmission_levels], fp)
else:
    fout = '{f}_{g}_{idx}'.format(f=transmission_configs["fecal_oral_dose"],
                                                            g=transmission_configs["beta_global"], idx=idx)
    output_dir = os.path.join(dir_path, "..", "Taniuchi_epidemics", "global2", '')  # dump files in lhs_sampling folder
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    with open(output_dir + 'epidemic_' + fout + '.json', 'w') as fp:
       json.dump([prevalences, sim_cohort_individuals, sim_cohort_incidences, transmission_levels], fp)


#with open(os.path.join(dir_path, "..", "output", '') + 'transmissions.pkl', 'wb') as fp:
#    dill.dump(transmissions, fp)
#write_dill(D, basename='Matlab_Villages', output_dir=output_dir)
