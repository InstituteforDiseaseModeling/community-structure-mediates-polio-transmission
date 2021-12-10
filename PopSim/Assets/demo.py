import os
import sys
import json
import config
from collections import defaultdict
from Demographics import Demographics, write_dill
import dill
import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "test/")
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
            
            
            
def set_configs(type = 'multiscale', ):
    if type == 'multiscale':
        c  =  {'n_initial_villages' : 5, 'n_initial_baris' : 200,
        "fecal_oral_dose" : 3.5e-6,
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
        c = {'n_initial_villages': 5, 'n_initial_baris' : 200,
        "fecal_oral_dose" : 3.5e-6,
        "hh_transmission" : 0,
        "bari_transmission" : 0,
        "village_transmission" : 0,
        "global_transmission" : 1,
        "beta_global": 19, #1
        "age_assortivity": 0}
                       
    #user specified method to overwrite default params                   
    config.params.overwrite(c)  
    for key in c:
        print(key + ' changed to ' + str(c[key]))


def read_D_file(sim_files):
    #sim_files is an array specifying the path of all simulation files
    random_idx = np.random.randint(0, len(sim_files) -1 )
    sim_file = sim_files[random_idx]
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)
    individuals = []
    hh_ids = []
    for V in D.villages:
        for bari in V.baris.values():
            for hhid in bari.households:
                hh_ids.append(hhid)
        individuals += V.return_individuals()    
        
    #force id counts to be reset
    start_iid = np.max([i.id for i in individuals])
    start_hid = np.max(hh_ids)
    Individual.override_id_generator(start_iid+10)
    Household.override_hid_generator(start_hid+10)
    return D
    
def run_demographics(basename = 'demo_villages'):
    D = Demographics()
    #create a random set of 20, equally sized villages
    D.initialize_default_villages()
    #output Demographics object (containing all villages, households, etc etc) to file
    write_dill(D, basename='demo_villages', output_dir=output_dir)
    return D

def Taniuchi_study_demo(params = config.params, sim_files = None, n_iterations = 1, output_folder_basename = 'test'):   
    #various recording stats -- each row will represent the results from a single iteration
    prevalences = defaultdict(list)
    transmissions = [] #transmission tree
    sim_cohort_incidences = defaultdict(list)        
    #sim_files should be an array specifying the path of pregenerated simulation files
    for _ in range(n_iterations): #number of iterations to do
        if sim_files:
            D = read_D_file(sim_files)
        else:
            D = run_demographics()

        #run the mOPV2 challenge study described by the Taniuchi study
        D.Taniuchi_study(params)
        
        for key in D.sim_prevalences:#prevalence is split into one of the 8 cohorts defined by the Taniuchi study
            prevalences[key].append(D.sim_prevalences[key]) #D.sim_prevalences records time series data
        for key in D.cohort_incidences:
            sim_cohort_incidences[key].append(D.cohort_incidences[key])
            
        transmissions.append(D.transmission_tracker)
        
    with open(output_dir + 'epidemic.json', 'wb') as fp:
        dill.dump(prevalences, fp)

    with open(output_dir + 'incidences.json', 'wb') as fp:
        dill.dump(sim_cohort_incidences, fp)  
    
    with open(output_dir + 'transmissions.json', 'wb') as fp:
        dill.dump(transmissions, fp)  
              


def outbreak_demo(params = config.params, sim_files = None, p = 0, n_iterations = 1, output_folder_basename = 'test'): 
    #initial conditions = 0 or 1, 0 to start single infant outbreak, 1 to do percentage based outbreak
    prevalences = defaultdict(list)
    transmissions = []
    sim_cohort_incidences = defaultdict(list)
    for _ in range(1): 
        if sim_files:
            D = read_D_file(sim_files)
        else:
            D = run_demographics()
        D.set_outbreak_conditions(strain_type='S2')
        print('Starting iteration: {i}'.format(i=_))       
        if p == 0:
            D.outbreak(p = 0.0, children = False, single_infant = True, params = params, strain_type = 'S2')
        else:
            D.outbreak(p = p, children = True, single_infant = False, params = params, strain_type = 'S2')                
                
        for key in D.sim_prevalences:
            prevalences[key].append(D.sim_prevalences[key])
        for key in D.incidence:
            sim_cohort_incidences[key].append(D.incidence[key])
                
        transmissions.append(D.transmission_tracker)
                
    with open(output_dir + 'epidemic.json', 'wb') as fp:
        dill.dump(prevalences, fp)

    with open(output_dir + 'incidences.json', 'wb') as fp:
        dill.dump(sim_cohort_incidences, fp)  

    with open(output_dir + 'transmissions.json', 'wb') as fp:
        dill.dump(transmissions, fp) 
    
if __name__ == '__main__':
    transmission_type = sys.argv[1]
    run_type = sys.argv[2]
    set_configs(transmission_type)
    
    if run_type == 'outbreak':
        outbreak_demo(params = config.params)
    elif run_type == 'Taniuchi':
        Taniuchi_study_demo(params = config.params)
    