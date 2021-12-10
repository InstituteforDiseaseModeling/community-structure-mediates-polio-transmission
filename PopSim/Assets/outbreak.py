import os
import sys
from collections import defaultdict
from Demographics import Demographics, write_dill
import dill
import config
import json
import copy
import numpy as np
from Individual import Individual
from Household import Household

dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "..", "outbreak_paper_sims")
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

year = sys.argv[1]
param_idx = int(sys.argv[2])
param_file = sys.argv[3]
attempt = sys.argv[4]
with open(param_file) as fin:
    parameters = json.load(fin)

print(parameters[year][param_idx])
p, strain_type, year, type, iteration = parameters[year][param_idx]

fout = '{strain_type}_{p}_{year}_{type}_{iter}.{attempt}'.format(strain_type = strain_type, p = str(p), year = str(year),
                                                                 type = type, iter = str(iteration), attempt = str(attempt))

if float(p) == 0.0:
    initial_conditions = 0
else:
    initial_conditions = 1

Matlab_village_files = search_files("/n/home04/weswong/multiscale_polio/MultiscaleModeling/PopSim/Demographics/Outbreak_Villages",'_{year}.pkl'.format(year = year))
#print(Matlab_village_files)
multiscale_transmission_configs = {"fecal_oral_dose" : 1.25e-06,
                        "hh_transmission" : 1,
                    "beta_hh" : 1,
                    "bari_transmission" : 1,
                    "beta_bari" : 15,
                    "village_transmission" : 1,
                    "beta_village" : 4,
                    "age_assortivity": 0,
                    "beta_inter_village":2,
                    "global_transmission" : 0,
                    'strain_type': strain_type,
                    }
global_transmission_configs ={"fecal_oral_dose" : 1.25e-06,
                    "hh_transmission" : 0,
                    "bari_transmission" : 0,
                    "village_transmission" : 0,
                    "global_transmission" : 1,
                    "beta_global": 19,
                    "age_assortivity": 0,
                    'strain_type': strain_type,}
                
if type == 'global':
    transmission_configs = global_transmission_configs
else:
    transmission_configs = multiscale_transmission_configs
                    
                    
config.params.overwrite(transmission_configs)   

                              
prevalences = defaultdict(list)
r0, transmissions = [], []
sim_cohort_individuals = defaultdict(list)
sim_cohort_incidences = defaultdict(list)
transmission_levels = defaultdict(list)
immunity = defaultdict(list)

for _ in range(1):
    print('Starting iteration: {i}'.format(i=_))
    random_idx = np.random.randint(0, len(Matlab_village_files) -1 )
    village_object_file = Matlab_village_files[random_idx] #Matlab_village_files[random_idx]
    cohort_individuals = defaultdict(list)
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)
    individuals = []
    n_hh = 0
            
    individuals = []
    hh_ids = []
    for V in D.villages:
        for bari in V.baris.values():
            for hhid in bari.households:
                hh_ids.append(hhid)
        individuals += V.return_individuals()
    start_iid = np.max([i.id for i in individuals])
    start_hid = np.max(hh_ids)
    Individual.override_id_generator(start_iid+10)
    Household.override_hid_generator(start_hid+10)

            
            
    if strain_type == 'WPV':
        for V in D.villages: #hack to make it WPV, WPV immunity assumed to be X-reactive bOPV2 induced S2 immunity
            for i in V.return_individuals():
                i.infection.strain_type = 'WPV'
            
    if str(initial_conditions) == str(0):
        print('point')
        shed_duration = D.outbreak(p = 0.0, children = False, single_infant = True, params = config.params, strain_type = strain_type, second_campaign = False)
    else:
        D.outbreak(p = p, children = True, single_infant = False, params = config.params, strain_type = strain_type, second_campaign = False)
        shed_duration = None
            
            
    for key in D.sim_prevalences:
        prevalences[key].append(D.sim_prevalences[key])
            
    for key in D.n_events_per_level:
        transmission_levels[key].append(D.n_events_per_level[key])
    for key in D.sim_immunity:
        immunity[key].append(D.sim_immunity[key])
            #transmissions.append(D.transmission_tracker)


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
            np.int16, np.int32, np.int64, np.uint8,
            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
            np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)): #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

with open(os.path.join(output_dir, '') + 'outbreak_' + fout + '.json', 'w') as fp:
    dumped = json.dumps(prevalences, cls = NumpyEncoder)
    json.dump(dumped, fp)

#with open(os.path.join(output_dir, '') + 'immunity_' + fout + '.json', 'w') as fp:
#    dumped = json.dumps(immunity, cls=NumpyEncoder)
#    json.dump(dumped, fp)

#if len(D.sim_immunity[key]) >= 360:
#    with open(os.path.join(output_dir, '') + 'demographics_' + fout + '.json', 'w') as fp:
#        dumped = json.dumps(D.demographics, cls=NumpyEncoder)
#        json.dump(dumped, fp)