#Tree node for household model
from collections import Counter
import config

class Node:
    def __init__(self):
        self.individuals = {}
        self.status = 'dead'
        self.successor = False
    
    @classmethod
    def create_node(cls, iids, isexes): #iids = individual id list, isexes = individual sex list
        N= Node()
        for iid,sex in zip(iids, isexes):
            N.individuals[iid] = sex
        if len(N.individuals.values()) > 0:
            if 0 in N.individuals.values() and 1 in N.individuals.values():
                N.status = 'married'
            elif 0 in N.individuals.values():
                N.status = 'single_female'
            else:
                N.status = 'single_male'
        else:
            N.status = 'dead'
        return N
    
    def remove_individual(self, individual_id):
        self.individuals.pop(individual_id)
        self.update()
    
    def add_member(self, individual_id, individual_sex):
        self.individuals[individual_id] = individual_sex
        self.update()
        
    def n_preferred_sex(self):
        '''find number of favored sex in the node
        patriarchical = 1, matriarchical = 0'''
        n_preferred_sex = Counter(self.individuals.values())[config.params.preferred_sex]
        return n_preferred_sex                 
    
    def update(self):
        '''only updated when some aspect of the node is affected (ie adding a new member or deleting a new member)'''
        if len(self.individuals.values()) > 0:
            if 0 in self.individuals.values() and 1 in self.individuals.values():
                self.status = 'married'
            elif 0 in self.individuals.values():
                if self.status == 'married':
                    self.status = 'widow'
                elif self.status != 'widow':
                    self.status = 'single_female'
            else:
                if self.status == 'married':
                    self.status = 'widower'
                elif self.status != 'widow':
                    self.status = 'single_male'
        else:
            self.status = 'dead'