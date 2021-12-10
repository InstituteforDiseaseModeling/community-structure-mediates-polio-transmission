import itertools
import numpy as np
from collections import defaultdict, Counter
from Individual import Individual, Female, Male
from Infection import ImmunoInfection
from Household import Household
import Utility
from Bari import Bari
import config 
import copy
import math
import multiprocessing as mp

class Village:
    vid_generator = itertools.count()
    def __init__(self):
        self.baris = {} #either household or bari baris
        self.id = next(Village.vid_generator)
        self.singles = defaultdict(list)
        self.f_marriage_ages = []
        self.m_marriage_ages = []
        self.time = 0
        self.bari_age_structure  = {}
        self.infants = defaultdict(list) # shallow copies...python will not recursively copy objects unless specifically told to
        self.hh_contacts = defaultdict(list)
        self.community_cases = [] #5% random sampling
        
    @classmethod
    def initialize_village(cls, n_baris = 5, n_households = 1):
        V = Village()
        for _ in range(n_baris): 
            new_bari = Bari.initialize_bari(n_households)
            V.baris[new_bari.id]=new_bari
        return V
    
    def age_matrix(self):
        age_matrix = {}
        for bari in self.baris.values():
            age_matrix[bari.id] = bari.age_structure()
        return age_matrix

    def find_singles(self):
        self.singles = {}
        single_females, single_males = [],[]
        for bari in self.baris.values():
            s_females, s_males = bari.find_singles()
            single_females += s_females
            single_males += s_males
        self.singles[1] = single_males
        self.singles[0] = single_females
        
    def create_unions(self):
        if len(self.singles[0]) > 0:
            for f,m in itertools.zip_longest(self.singles[0], self.singles[1]):
                if f != None and m != None:
                    f_bari_id, f_hhid, f_nid, f_iid = f[0], f[1], f[2], f[3]
                    m_bari_id, m_hhid, m_nid, m_iid = m[0], m[1], m[2], m[3] #assortative mating
                   
                    if f_bari_id != m_bari_id:
                        #Remove Female from her former household
                        F_individual = self.baris[f_bari_id].households[f_hhid].individuals.pop(f_iid)
                        self.f_marriage_ages.append((self.time, F_individual.age))
                            
                        self.baris[f_bari_id].households[f_hhid].reverse_tree_map.pop(f_iid)
                        self.baris[f_bari_id].households[f_hhid].tree_map[f_nid].remove_individual(f_iid)
                            
                        #add female to her new household
                        self.baris[m_bari_id].households[m_hhid].new_union(m_nid, F_individual)
                        self.m_marriage_ages.append((self.time, self.baris[m_bari_id].households[m_hhid].individuals[m_iid].age))
            
    def calculate_marriage_differential(self):
        marriage_differential = []
        for bari in self.baris.values():
            for hh in bari.households.values():
                marriage_differential += hh.find_marriage_differential()
        return marriage_differential
    
    def hh_sizes(self):
        hh_sizes = []
        for bari in self.baris.values():
            for hh in bari.households.values():
                hh_sizes.append(hh.n_individuals())
        return hh_sizes
    
    def bari_n_households(self):
        bari_sizes = []
        for bari in self.baris.values():
            bari_sizes.append(bari.n_households())
        return bari_sizes
    
    def bari_n_individuals(self):
        individuals_per_bari = []
        for bari in self.baris.values():
            n_individuals = 0
            for hh in bari.households.values():
                n_individuals += hh.n_individuals()
            individuals_per_bari.append(n_individuals)
        return individuals_per_bari
                    
    
    def infant_bari_count(self):
        infant_baris = []
        for bari in self.baris.values():
            for hh in bari.households.values():
                for individual in hh.individuals.values():
                    if individual.age <= config.params.infant_age:
                        infant_baris.append(bari.id)
        infant_bari_count = list(Counter(infant_baris).values())
        return infant_bari_count
        
    def noinfant_bari_prop(self):
        noinfant_bari_count = 0
        total_bari_count = 0
        for bari in self.baris.values():
            infant_count = 0
            for hh in bari.households.values():
                for individual in hh.individuals.values():
                    if individual.age <= 1:
                        infant_count += 1
            if infant_count == 0:
                noinfant_bari_count += 1
            total_bari_count += 1
        if total_bari_count:
            return noinfant_bari_count/total_bari_count
        else:
            return 0
    
    def n_baris(self):
        return len(self.baris)
    
    def update(self, dt, params = config.params):
        self.time += dt
        dead_uids = []
        self.transmissions = []
        self.demographics_tracker = []
        for bari in self.baris.values():
            transmissions, demographics = bari.update(dt, params)
            self.transmissions += transmissions
            for D in demographics:
                D.village_id = self.id
            self.demographics_tracker += demographics
            if bari.isalive == False:
                dead_uids.append(bari.id)
                
        for uid in dead_uids:
            self.baris.pop(uid)
        
        if len(self.transmissions) > 0:
            for Tx in self.transmissions:
                Tx.village_id = self.id
            if params.village_transmission:
                self.find_individuals() #function used for updating inter village stats
                self.inter_bari_transmission()
            

        self.find_singles()
        self.create_unions()
        self.check_status()
    
    def find_individuals(self):
        '''used as an individual tracker for transmission'''
        self.individuals, self.bari_tracker, self.individual_ages = [], [], []
        for bari in self.baris.values():
            bari_individuals = bari.return_individuals()
            self.individuals += bari_individuals
            self.bari_tracker += len(bari_individuals) * [bari.id]
            self.individual_ages += [int(round(I.age,0)) if I.age < 80 else 80 for I in bari_individuals]
        self.individuals = np.asarray(self.individuals)
        self.individual_ages = np.asarray(self.individual_ages)
        self.bari_tracker = np.asarray(self.bari_tracker)
    
    def update_individuals(self):
        birth_individuals,birth_ages, birth_baris = [], [], []

        death_iids = []
        if len(self.demographics_tracker) > 0:
            for D in self.demographics_tracker:
                if D.type: #birth event occured
                    birth_individuals.append(D.individual)
                    birth_ages.append(int(max(round(D.individual.age, 0), 80)))
                    birth_baris.append(D.bari_id) 
                else: #death event occured
                    death_iids.append(D.individual.id)
                    
            if len(death_iids) > 0:
                mask = np.asarray([individual.id not in death_iids for individual in self.individuals])
                self.individuals = self.individuals[mask]
                self.individual_ages = self.individual_ages[mask]
                self.bari_tracker = self.bari_tracker[mask]
                
            np.append(self.individuals, birth_individuals)
            np.append(self.individual_ages, birth_ages)
            np.append(self.bari_tracker, birth_baris)
        
    def inter_bari_transmission(self, params = config.params):
        #output = mp.Queue()
        for Tx in self.transmissions: 
            source_bari = Tx.bari_id
            mask = self.bari_tracker != source_bari
            
            infectees = self.individuals[mask]
            infectee_ages = self.individual_ages[mask]
            ego_age = Tx.source_age
            contact_probabilities = config.params.village_age_assortivity_matrix[ego_age][infectee_ages]
            cum_contact_probabilities = np.sum(contact_probabilities)
            r_prob = contact_probabilities / cum_contact_probabilities
            potential_infectee_indexes = range(len(infectees))

            viral_dose = Tx.dose * params.fecal_oral_dose
            if len(infectees) > 0:
                infectee_idxes = np.random.choice(potential_infectee_indexes, params.beta_village, p = r_prob)
                for infectee_idx in infectee_idxes:
                    infected = infectees[infectee_idx]
                    infected_status, prechallenge_immunity = infected.receive_transmission_dose(viral_dose)
                    Tx.secondary_cases += infected_status
                    if infected_status:
                        Tx.tx_history['village'].append((infected.id, infected.age, prechallenge_immunity, infected.times_infected))          

                        
                        
    
    def track_individuals(self,iid):
        if iid in [i.id for i in V.return_individuals()]:
            for bari in V.baris.values():
                if iid in [i.id for i in bari.return_individuals()]:
                    for hh in bari.households.values():
                        if iid in [i.id for i in hh.individuals.values()]:
                            return [hh.id, bari.id, self.id]
            
    
    
    
    def check_status(self):
        if len(self.baris) > 1:
            self.isalive = True
        else:
            self.isalive = False
    
    def vaccination_campaign(self, trial, vaccine_dose = 10**6, initialize_pre_existing_immunity = True):
        '''Taniuchi et al design
        trial: bOPV vs tOPV'''
        self.community_cases = []
        self.infants = defaultdict(list) # shallow copies...python will not recursively copy objects unless specifically told to
        self.hh_contacts = defaultdict(list)
        
        for bari in self.baris.values():
            bari_vaccination_flag = False
            bari_infants = defaultdict(list)
            bari_hh_contacts = defaultdict(list)
            bari_community_cases = []
            for hh in bari.households.values():
                infants, hh_contacts, other_hh_members = hh.return_Taniuchi_individuals()
                n_infants = len(infants)
                infant_status = None
                infant_enrolled = False
                if n_infants > 0:
                    for infant in infants:
                        if initialize_pre_existing_immunity:
                            infant.infection = ImmunoInfection.initialize_infant_immunity(trial)
                        if np.random.random() < 1/2:
                            infant_enrolled = True
                            if np.random.random() < 0.33:
                                infant.receive_transmission_dose(10**6, index=True) #mOPV2
                                bari_infants['i+'].append(infant)
                                bari_vaccination_flag = True
                                infant_status = True
                            else:
                                bari_infants['i-'].append(infant)
                                infant_status = False
                        
                        if n_infants > 1:
                            infant_status = 'mult'
                    
                    for hh_contact in hh_contacts:
                        if initialize_pre_existing_immunity:
                            hh_contact.infection = ImmunoInfection.initialize_noninfant_immunity(hh_contact.age, trial)
                        if infant_enrolled:
                            if np.random.random() < 0.05:
                                hh_contact.receive_transmission_dose(10**6, index = True)
                                bari_vaccination_flag = True
                                if infant_status == True:
                                    bari_hh_contacts['hh+_i+'] += [hh_contact]
                                elif infant_status == False:
                                    bari_hh_contacts['hh+_i-'] += [hh_contact]
                                else:
                                    bari_hh_contacts['mult_infants'] += [hh_contact]
                            else:
                                if infant_status == True:
                                    bari_hh_contacts['hh-_i+'] += [hh_contact]
                                elif infant_status == False:
                                    bari_hh_contacts['hh-_i-'] += [hh_contact]
                                else:
                                    bari_hh_contacts['mult_infants'] += [hh_contact]
                        else:
                            if hh_contact.age <=5:
                                hh_contact.receive_transmission_dose(10**6, index = True)
                                bari_community_cases.append(hh_contact)
                                bari_vaccination_flag = True
                    
                    
                    for other_individual in other_hh_members:
                        if initialize_pre_existing_immunity:
                            other_individual.infection = ImmunoInfection.initialize_noninfant_immunity(other_individual.age, trial)
                        if other_individual.age <=5:
                            if np.random.random() < 0.40:
                                other_individual.receive_transmission_dose(10**6, index = True)
                                bari_community_cases.append(other_individual)
                                bari_vaccination_flag = True
                    

                else: # no infants in the household
                    for individual in hh.individuals.values():
                        if initialize_pre_existing_immunity:
                            individual.infection = ImmunoInfection.initialize_noninfant_immunity(individual.age, trial)
                        if individual.age <= 5:
                            if np.random.random() < 0.40:
                                individual.receive_transmission_dose(10**6, index = True)
                                bari_community_cases.append(individual)
                                bari_vaccination_flag = True
            self.community_cases += bari_community_cases
            
            for key in ['hh+_i+', 'hh+_i-', 'hh-_i+', 'hh-_i-']:
                if bari_vaccination_flag:
                    self.hh_contacts['b+_' + key] += bari_hh_contacts[key]
                else:
                    self.hh_contacts['b-_' + key] += bari_hh_contacts[key]

            for key in ['i+', 'i-']:
                if bari_vaccination_flag:
                    self.infants['b+_' + key] += bari_infants[key]
                else:
                    self.infants['b-_' + key] += bari_infants[key]
                        
        self.find_individuals()
    def n_individuals(self, infants = False):
        n_individuals, n_infants = 0,0
        if infants:
            for bari in self.baris.values():
                for hh in bari.households.values():
                    n_individuals += hh.n_individuals()
                    n_infants += hh.n_infants()
            return n_individuals, n_infants
        else:
            for bari in self.baris.values():
                for hh in bari.households.values():
                    n_individuals += hh.n_individuals()
            return n_individuals
    
    def n_households(self):
        n_households = 0
        for bari in self.baris.values():
            n_households += len(bari.households)
        return n_households
    
    def households(self):
        hh_list = []
        for bari in self.baris.values():
            for hh in bari.households.values():
                hh_list.append(hh)
                return hh_list
                
    def return_individuals(self):
        individual_list = []
        for bari in self.baris.values():
            for hh in bari.households.values():
                individual_list += list(hh.individuals.values())
        return individual_list
        
    
    
            
   