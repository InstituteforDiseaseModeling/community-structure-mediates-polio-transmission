import itertools
import numpy as np
from collections import defaultdict
import copy
from Individual import Individual, Female, Male
from Household import Household
from Node import Node
import config
import Utility

class Bari:
    bid_generator = itertools.count()
    def __init__(self):
        self.households = {}
        self.id = next(Bari.bid_generator)
        self.generation = 0
        self.unroot_events = []
        self.hh_age_structure = {} #age bin dictionary
        
    @classmethod    
    def initialize_bari(cls, n_households = 1):
        B = Bari()
        for _ in range(n_households):
            new_hh = Household.initialize_household([Female(config.params.f_marriage_age()), Male(config.params.m_marriage_age())])
            B.households[new_hh.id] = new_hh #start just shy of fertility age
            B.hh_age_structure[new_hh.id] = new_hh.age_structure()
        B.isalive = True        
        return B
    
    def find_singles(self):
        singles = defaultdict(list)
        for hh in self.households.values():
            single_males, single_females = hh.find_singles(self.id)
            singles[0] += single_females
            singles[1] += single_males                    
        return singles[0], singles[1]
    
    def n_individuals(self):
        count = 0
        for hh in self.households.values():
            count += hh.n_individuals()
        return count
    
    def return_individuals(self):
        individuals = []
        for hh in self.households.values():
            individuals += hh.individuals.values()
        return individuals
    
    def n_households(self):
        return len(self.households)
    
    def check_status(self):
        if len(self.households) >= 1:
            self.isalive = True
        else:
            self.isalive = False
            
    def split_bari(self, hh_id):
        """dummy code, doesn't do anything right now"""
        removed_hh = self.households.pop(hh_id)
        return removed_hh
    
    def unroot_hh(self, hh_id):
        self.unroot_events.append(self.generation)
        hh = self.households[hh_id]
        if hh.initiate_unroot == True: #todo: parameter of splitting after patriarch dies, would have to to be modified to include an if root empty condition
            hh.isactive = False #the original hh is being split and destroyed
            new_hhs, true_orphans = hh.unroot()
            for new_hh in new_hhs:
                self.households[new_hh.id] = new_hh
            
            for orphan in true_orphans:
                orphan_object, status = orphan
                adopter_id = np.random.choice(list(self.households.keys()))
                self.households[adopter_id].external_adoption('0', orphan_object, status) #assign orphan to the top level node...?
    
    def elder_adopts(self):
        elders = []
        potential_adopters = [hhid for hhid in self.households.keys() if len(self.households[hhid].individuals) >=2]
        if len(potential_adopters) > 0:
            for hh in self.households.values():
                if len(hh.individuals) == 1:
                    iid = list(hh.individuals.keys())[0]
                    if hh.individuals[iid].age > config.params.elder_age:
                        nid = hh.reverse_tree_map.pop(iid)
                        node_status = copy.deepcopy(hh.tree_map[nid].status)

                        elders.append(hh.individuals.pop(iid))
                        hh.tree_map[nid].remove_individual(iid)
                        hh.isactive = False
        for elder in elders:
            adopter_hhid = np.random.choice([hhid for hhid in self.households.keys() if len(self.households[hhid].individuals) >=2])
            self.households[adopter_hhid].external_adoption('0', elder, 'adopted_elder')
                    
    def bari_transmission(self, transmissions, params = config.params):
        hh_individuals, hh_tracker = [],[]
        for hh in self.households.values():
            hh_individuals += hh.individuals.values()
            hh_tracker += len(hh.individuals) * [hh.id]
        hh_individuals = np.asarray(hh_individuals)
        hh_tracker = np.asarray(hh_tracker)
        hh_individual_ages = np.asarray([min(int(round(individual.age,0)), 80) for individual in hh_individuals])
        for Tx in transmissions:
            source_hh = Tx.hhid
            mask = hh_tracker != source_hh
            infectees = hh_individuals[mask]
            potential_infectee_indexes = range(len(infectees)) 
            infectee_ages = hh_individual_ages[mask]
            
            ego_age = Tx.source_age
            contact_probabilities = config.params.bari_hh_age_assortivity_matrix[ego_age][infectee_ages]
            cum_contact_probabilities = np.sum(contact_probabilities)
            r_prob = contact_probabilities / cum_contact_probabilities
            
            viral_dose = Tx.dose * params.fecal_oral_dose
            if len(infectees) > 0 :
                infectee_iidxes = np.random.choice(potential_infectee_indexes, params.beta_bari, p = r_prob)
                for infectee_idx in infectee_iidxes: 
                    infected = infectees[infectee_idx]
                    infected_status, prechallenge_immunity = infected.receive_transmission_dose(viral_dose)
                    Tx.secondary_cases += infected_status
                    if infected_status:
                        Tx.tx_history['bari'].append((infected.id, infected.age, prechallenge_immunity, infected.times_infected))

                    
    def update(self, dt=1, params = config.params):
        self.generation += dt
        dead_hhids = []
        unroot_hhids = []
        new_hhs = []
        self.hh_age_structure = {}
        self.elder_adopts()
        transmissions = []
        
        
        for hh in self.households.values():
            new_hh_list, hh_tx, demographics_tracker = hh.update(dt, params) 
            new_hhs += new_hh_list
                
            if hh.isactive == False: #unrooting trumps this
                dead_hhids.append(hh.id)
            elif hh.initiate_unroot == True:
                unroot_hhids.append(hh.id)
                transmissions += hh_tx
           
            else:
                transmissions += hh_tx
                
        for hhid in dead_hhids:
            self.households.pop(hhid)
            
        if len(transmissions) > 0:
                for Tx in transmissions:
                    Tx.bari_id = self.id
                    
        if params.bari_transmission:
            self.bari_transmission(transmissions)
        
        #unrooting the household is last to occur...mainly because it's a pain to track individuals after they've been split
        if len(unroot_hhids) > 0:
            for hhid in unroot_hhids:
                self.unroot_hh(hhid)
                self.households.pop(hhid)
            
        for hh in new_hhs:
            self.households[hh.id] = hh
            
        self.check_status()
        
        for D in demographics_tracker:
            D.bari_id = self.id
        return transmissions, demographics_tracker
        
    def render_tree(self, base_fp=None):
        count = 1
        for hh in self.households.values():
            t = hh.tree    
            for node in t.traverse():    
                nstyle = NodeStyle()
                nstyle["size"] = 20
                n =hh.tree_map[node.name]        

                if n.status == 'married':
                    nstyle["shape"] = "sphere"
                    nstyle["fgcolor"] = "DarkRed"

                elif n.status == 'widow': #widow is female
                    nstyle["shape"] = "circle"
                    nstyle["fgcolor"] = 'lightpink'

                elif n.status == 'widower':
                    nstyle["shape"] = "square"
                    nstyle["fgcolor"] = 'lightblue'

                elif n.status == 'single_male':
                    nstyle["shape"] = "square"
                    nstyle["fgcolor"] = 'DodgerBlue'

                elif n.status == 'single_female':
                    nstyle["shape"] = "circle"
                    nstyle['fgcolor'] = 'Violet'

                else:
                    nstyle["shape"] = "circle"
                    nstyle['fgcolor'] = 'grey'

                node.set_style(nstyle)

            ts = TreeStyle()
            ts.rotation = 90
            ts.branch_vertical_margin  = 10
            ts.show_leaf_name = False
            ts.show_branch_length = False
            ts.show_branch_support = False
            ts.show_scale = False
            
            if base_fp:
                t.render(base_fp + '_' +str(count) + '.svg', tree_style=ts)
                count += 1
            else:
                t.show(tree_style = ts)

