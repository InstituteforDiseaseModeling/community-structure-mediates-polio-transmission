import itertools
import numpy as np
from collections import defaultdict, Counter
from Individual import Individual, Female, Male
from Node import Node
from ete3 import Tree#, TreeStyle, NodeStyle
import copy
from Utility import remove_subtree, add_child, print_tree, retrieve_root_distance, relabel_tree,find_siblings
import Utility
import config
	
#Household Class	
class Household:
    '''Household defined by active couple units that generate children --
    Allow for multi-generational living (grandparents, children, and grandchildren cohabitating)
    Centered around the concept of a family tree: family trees are generated for each household
    explicitly keeps track of lineages
    loosely based on family systems defined by William Skinner in CH 3 of anthropological demography toward a new synthesis'''
    hid_generator = itertools.count()
    
    def __init__(self):
        '''4 layers of organization that must be continuously updated:
        1) family tree
        2) tree_map that relates nodes to individuals: more accurately: the Node Class
        3) reverse_tree_map that relates individuals to nodes: only stores iid
        4) individual dictionary that holds the individuals: only place where Individual Class resides'''
        self.individuals = {} # a container to hold the individual classes relevant to this household
        self.id = next(Household.hid_generator) 
        self.tree_map = {}
        self.nodeid_generator = itertools.count()
        self.reverse_tree_map = {}
        self.isactive=True
        self.initiate_unroot = False
        self.n_transmitting = 0
        self.index_cases = []

        #self.p_above15_moveout = 0.0
        
    @classmethod
    def override_hid_generator(cls, start):
        Household.hid_generator = itertools.count(start)
        
            
    def set_nodid_generator(self, start):
        self.nodeid_generator = itertools.count(start)
    
    @classmethod
    def initialize_household(cls, founder_list):
        '''initialize household from scratch -- ie only the founder members'''
        HH = Household()
        for individual in founder_list:
            HH.individuals[individual.id] = individual
            individual.ismarried=True
            individual.t_since_marriage = 0
            individual.t_since_marriage_moveout = None
        
        HH.male_root = True #is there a male alive in the tree root
        node_name = str(next(HH.nodeid_generator))
        HH.tree = Tree(name=node_name) #true lineage mapping: dead individuals not removed
        HH.tree_map[node_name] = Node.create_node(list(HH.individuals.keys()), [individual.sex for individual in HH.individuals.values()])
        for iid in HH.individuals.keys():
            HH.reverse_tree_map[iid] = node_name
        return HH
    
    def birth(self, params = config.params):
        recent_mothers = []
        demographics_tracker = []
        for iid in self.individuals:
            if self.individuals[iid].sex == 0:
                if self.individuals[iid].just_gave_birth:
                    recent_mothers.append(iid)
        for f_id in recent_mothers:
            parental_node = self.reverse_tree_map[f_id]
            child_node = str(next(self.nodeid_generator))
            self.tree = add_child(self.tree, child_node, parental_node=parental_node) #update the tree
            if np.random.random() < params.p_male_child:
                new_child = Male(0, strain_type = params.strain_type)
            else:
                new_child = Female(0,  strain_type = params.strain_type)
            self.individuals[new_child.id] = new_child #update individual dictionary
            self.tree_map[child_node] = Node.create_node([new_child.id],[new_child.sex]) #update tree map
            self.reverse_tree_map[new_child.id] = child_node #update reverse tree map
            
            self.individuals[f_id].children_stats[new_child.sex].append(new_child.id)#temp?
            demographics_tracker.append(Utility.Demographics(new_child, 1))
        return demographics_tracker
    
    def remove_dead(self, dead_ids):
        for dead_iid in dead_ids:
            individual = self.individuals.pop(dead_iid) #remove dead individuals from household
            
            node=self.reverse_tree_map.pop(dead_iid) #check what node that individual belonged to on the tree
            self.tree_map[node].remove_individual(dead_iid) #remove that individual from the node and update node status
            
            if node == '0': # if the dead individual occured in the root
                root_node = self.tree_map['0']
                if config.params.preferred_sex is not None:
                    if root_node.n_preferred_sex() == 0:
                        #this should only initiate if the patriarch dies, if a widow is acting as the root then it should still 
                        #trigger when she dies
                        self.initiate_unroot = True #unroot the tree
                else:
                    if root_node.status == 'dead':
                        self.initiate_unroot = True
            
            if self.tree_map[node].status != 'married': #if no longer married, reset individual status to not married
                for iid in self.tree_map[node].individuals:
                    self.individuals[iid].ismarried = False
                    self.individuals[iid].t_since_marriage = None
    
    #adding new members ------------------------------------------
    def new_union(self, node, individual):
        self.individuals[individual.id] = individual
        self.tree_map[node].add_member(individual.id, individual.sex)
        self.reverse_tree_map[individual.id] = node
        if self.tree_map[node].status == 'married': #the married flag is controlled by the node object, considers 1 of each sex as being married
            parental_node_id = (self.tree&str(node)).up.name
            if self.tree_map[parental_node_id].successor == False:
                #move_out_date = 1e8
                move_out_date = config.params.hh_move_out_fn()
                self.tree_map[parental_node_id].successor = True
            else:
                move_out_date = config.params.hh_move_out_fn()
            for iid in self.tree_map[node].individuals:
                self.individuals[iid].ismarried = True
                self.individuals[iid].t_since_marriage = 0
                self.individuals[iid].t_since_marriage_moveout = move_out_date
                self.individuals[iid].marriage_age.append(self.individuals[iid].age)
    
    def internal_adoption(self, parental_node, orphan_nodeid):
        '''within_hh adoption procedure: uses nodeid because it just rewires the relationship between nodes
        only the nodeid is required'''
        #print_tree(self.tree)
        self.tree = remove_subtree(self.tree, orphan_nodeid)
        #print_tree(self.tree)
        self.tree = add_child(self.tree, orphan_nodeid, parental_node = parental_node)
        #print_tree(self.tree)
        
    def external_adoption(self, parental_node, orphan_individual, status = None):
        '''external_hh adoption procedure: more involved because previous node relationships are meaningless
        orphan individual = individual class'''
        self.individuals[orphan_individual.id] = orphan_individual
        child_node = str(next(self.nodeid_generator))
        self.tree = add_child(self.tree, child_node, parental_node = parental_node)
        self.tree_map[child_node] = Node.create_node([orphan_individual.id], [orphan_individual.sex])
        self.reverse_tree_map[orphan_individual.id] = child_node
        
        if status: #used for widow
            self.tree_map[child_node].status = status
            
    def merge_households(self, hh2):
        for node in hh2.tree.traverse('levelorder'):
            new_nid = str(next(self.nodeid_generator))
            self.tree_map[new_nid] = hh2.tree_map.pop(node.name)
            for iid in self.tree_map[new_nid].individuals:
                hh2.reverse_tree_map.pop(iid)
                self.reverse_tree_map[iid] = new_nid
                self.individuals[iid] = hh2.individuals.pop(iid)
            node.name = new_nid
        subtree = self.tree&'0'
        subtree.add_child(hh2.tree)
        self.tree = subtree
        hh2.isactive = False               
                
            
    
    #splitting the tree -----------------------------
    def split_household(self,split_node):
        '''split the tree at the split_node. The split_node is treated as the new tree root
        all downstream individuals are carried with it'''
        self.tree, split_tree = remove_subtree(self.tree, split_node, return_subtree = True)
        split_HH = Household()
        
        count =0
        for node in split_tree.traverse('levelorder'):
            split_HH.tree_map[str(count)] = self.tree_map.pop(node.name) #transfer node to new split HH, reset so that root = 0
            for iid in split_HH.tree_map[str(count)].individuals:
                self.reverse_tree_map.pop(iid)
                split_HH.reverse_tree_map[iid] = str(count)
                split_HH.individuals[iid] = self.individuals.pop(iid) #transfer individual to new split HH
            node.name = str(count)
            count +=1
        split_tree.name='0'
        split_HH.tree = split_tree
        split_HH.set_nodid_generator(count)
        return split_HH
        
    def prepare_unroot(self):
        '''In the event that the root is unable to lead the household, the lower generations will split
        currently triggers when no male is present in the root. Old widows will not trigger an unrooting event
        Does not actually unroot the tree: identifies the viable (age>15) nodes in the next sublevel and splits them off into new household
        Orphans (those without living parental units) will be adopted by other parental units within the household'''
        
        root_tree = self.tree&'0'
        true_orphans = [] #non-independent individuals (elder nodes, children) with no living person in the hierarchy
        children_nodes = defaultdict(list)
        
        for n in root_tree.children:
            children_nodes[self.tree_map[n.name].status].append(n.name)
        
        #find innermost living nodes
        while len(children_nodes['dead']) > 0:
            dead_nodes = copy.deepcopy(children_nodes['dead'])
            children_nodes['dead'] = []
            for dead_node in dead_nodes:
                subroot = self.tree.search_nodes(name=dead_node)[0]
                for node in subroot.children:
                    children_nodes[self.tree_map[node.name].status].append(node.name)
        children_nodes.pop('dead', 0)
        split_nodes = []
        orphaned = []
        widows = []

        #identify independent CU and and orphans
        for key in children_nodes:
            if key in ['married']:#, 'widower']:
                split_nodes += children_nodes[key]
            elif 'single' in key:
                #self.individuals are single from this point on
                for node in children_nodes[key]:
                    for iid in self.tree_map[node].individuals.keys():
                        orphaned.append(node)
            else: #widows & widowers
                widows += children_nodes[key]
        
        #find a suitable adopter -- only married individuals are the adopter
        if 'married' in children_nodes:
            adopter = lambda: np.random.choice(children_nodes['married'])  
        else:
            adopter = None
            
        if len(orphaned) > 0:
            if adopter:
                for orphan_node in orphaned:
                    self.internal_adoption(adopter(), orphan_node)
            else:
                for orphan_node in orphaned:
                    orphan_node_object = self.tree_map.pop(orphan_node)
                    orphan_iid = [k for k in orphan_node_object.individuals.keys()][0]
                    self.reverse_tree_map.pop(orphan_iid)
                    true_orphans.append((self.individuals.pop(orphan_iid), None))
        
        if len(widows) > 0: #any lower level widows
            if adopter: #attempt to adopt into the family
                for widow_node in widows:
                    self.internal_adoption(adopter(), widow_node)
            else: #return to true orphan pool
                for widow_node in widows:
                    widow_object = self.tree_map.pop(widow_node)
                    widow_iid = [k for k in widow_object.individuals.keys()][0]
                    self.reverse_tree_map.pop(widow_iid)
                    true_orphans.append((self.individuals.pop(widow_iid), 'widow'))
        
        #assign widow to descendents 
        if len(self.tree_map['0'].individuals) > 0:
            widow_node = self.tree_map.pop('0')
            widow_iid = [k for k in widow_node.individuals.keys()][0]
            self.reverse_tree_map.pop(widow_iid) #remove all references to the widow
            
            widow_object = self.individuals.pop(widow_iid)
            adopter_flag = False
            for son_id in widow_object.children_stats[1]:
                if adopter_flag == False:
                    if son_id in self.individuals:
                        adopter_node = self.reverse_tree_map[son_id]
                        adopter_flag == True
            if adopter_flag == True:
                self.external_adoption(adopter_node, widow_object, 'widow')
            else:
                true_orphans.append((widow_object, 'widow'))          
            
        return split_nodes, true_orphans
        
    def unroot(self, inherit=False):
        new_hh_statuses = defaultdict(list)
        split_node_ids, true_orphans = self.prepare_unroot()
        new_hhs = []
        for split_node_id in sorted(split_node_ids, reverse = True):
            new_hhs.append(self.split_household(split_node_id))
        
        #in the event that, for whatever reason, you want one child to inherit the ENTIRE HOUSEHOLD: INCLUDING BROTHERS AND SISTERS
        #not currently used...use at own risk?
        if inherit and len(new_hhs) > 1:
            inheritor_idx = np.random.choice([idx for idx, x in enumerate(new_hhs)])
            inheritor_hh = new_hhs[inheritor_idx]
            for idx,hh in enumerate(new_hhs):
                if idx != inheritor_idx:
                    inheritor_hh.merge_households(hh)
            for orphan in true_orphans:
                orphan_object, status = orphan
                inheritor_hh.external_adoption('0', orphan_object, status)
            true_orphans = []
            new_hhs = [inheritor_hh]
                
        return new_hhs, true_orphans
    
    
    def hh_transmission(self, transmissions, params = config.params):
        individuals = np.asarray(list(self.individuals.values()))
        for Tx in transmissions:
            source = Tx.source_iid
            mask = np.asarray([individual.id for individual in individuals]) != source
            infectees = individuals[mask]
            if len(infectees) > 0:
                infectee_ages = np.asarray([min(int(round(individual.age,0)),80) for individual in infectees])
                ego_age = Tx.source_age
                contact_probabilities = config.params.bari_hh_age_assortivity_matrix[ego_age][infectee_ages]
                cum_contact_probabilities = np.sum(contact_probabilities)
                r_prob = contact_probabilities / cum_contact_probabilities
                potential_infectee_indexes = range(len(infectees))
                viral_dose = Tx.dose * params.fecal_oral_dose
                infectee_idxes = np.random.choice(potential_infectee_indexes, params.beta_hh, p = r_prob)
                for infectee_idx in infectee_idxes:
                    infected = infectees[infectee_idx]
                    infected_status, prechallenge_immunity = infected.receive_transmission_dose(viral_dose)
                    Tx.secondary_cases += infected_status
                    if infected_status:
                        Tx.tx_history['hh'].append((infected.id, infected.age, prechallenge_immunity, infected.times_infected))
                #print('hh tx {param}: infecting {id} with {viral_dose} from {source}'.format(param = params.beta_hh, id = self.individuals[infectee_id].id, viral_dose = viral_dose, source = Tx.source_iid))
    
    # update------------------------
    def update(self, dt = 1, params = config.params): #years
        move_out_nodes = []
        dead_ids = []
        transmissions = []
        age_structure = np.zeros((1, len(params.age_bins)))
        self.n_transmitting = 0
        self.find_active_index_cases()
        demographics_tracker = []
        for individual in self.individuals.values():
            individual.update(dt)
            if individual.isalive == False:
                dead_ids.append(individual.id)
                demographics_tracker.append(Utility.Demographics(individual, 0))
            else:
                if individual.sex == 1 and individual.check_move_out():
                    move_out_node = self.reverse_tree_map[individual.id]
                    if move_out_node != '0':
                        move_out_nodes.append(move_out_node)
                
                if individual.is_transmitting:
                    self.n_transmitting += 1
                    Tx = individual.return_transmission()
                    Tx.hhid = self.id 
                    transmissions.append(Tx)
                    
        if params.hh_transmission:
            self.hh_transmission(transmissions)
                
        self.remove_dead(dead_ids)
        demographics_tracker += self.birth(params = config.params)
        #infections occur before demographic changes
        #this choice is unjustifiable (but also not a major sticking point)...
        
        
        split_hhs = []
        try:
            for nodeid in sorted(move_out_nodes, reverse = True):
                split_hhs.append(self.split_household(nodeid))
        except: #in the event of a very very rare bug....
            print(nodeid, move_out_nodes, print_tree(self.tree))
            
        if len(self.individuals) == 0:
            self.isactive = False
        return split_hhs, transmissions, demographics_tracker
        
    #useful attribute functions--------------------------
    def find_active_index_cases(self):
        self.active_index_cases = [individual.id for individual in self.individuals.values() if (individual.index == True) and individual.infection.is_shed]
    
    def n_individuals(self):
        return len(self.individuals)
        
    def n_infants(self):
        return len([individual.id for individual in self.individuals.values() if individual.age <=config.params.infant_age] )
    
    def return_Taniuchi_individuals(self):
        infants,hh_contacts, other = [],[],[]
        for individual in sorted([individual for individual in self.individuals.values()], key = lambda x: x.age, reverse = True):
            if individual.age < config.params.infant_age:
                infants.append(individual)
            elif individual.age <= 14 or individual.sex == 0:
                if len(hh_contacts) < 2:
                    hh_contacts.append(individual)
                else:
                    other.append(individual)
            else:
                other.append(individual)
        return (infants, hh_contacts, other)
        
        
        if len(hh_contacts) >= 2:
            youngest_contacts = hh_contacts[0:2]
            other_contacts = hh_contacts[2:]
        else:
            youngest_contacts = hh_contacts
            other_contacts = []
        return (infants, youngest_contacts, other_contacts)
    
    def age_structure(self):
        '''returns a dictionary whose keys are the binned age indexes of the household. Refer to config.params.age_bins for the bin ranges'''
        age_struct = [0 for _ in range(len(config.params.age_bins))]
        for iid in self.individuals:
            age_idx = Utility.find_bin_idx(self.individuals[iid].age, config.params.age_bins)
            age_struct[age_idx] += 1
        return age_struct
            
    
    def find_eldest_node(self):
        eldest_idx = np.argmax([i.age for i in self.individuals.values()])
        eldest_node = self.reverse_tree_map[list(self.individuals.keys())[eldest_idx]]
        return eldest_node
    
    def family_type(self):
        #find the max distance between LIVING nodes
        '''find the max distance between the eldest node and the terminal living nodes (leaf)
        nuclear: distance <= 1
        joint/multi-generational: distance > 1'''
        leaves = self.tree.get_leaves()
        eldest_node = self.find_eldest_node()
        distances = [self.tree.get_distance(eldest_node, leaf.name) for leaf in leaves if self.tree_map[leaf.name].status != 'dead' and eldest_node != leaf.name]
        
        if len(self.tree_map) > 1 and len(distances) > 0:
            if np.max(distances) == 1:
                if len([n.status for n in self.tree_map.values() if n.status == 'married']) == 1:
                    return 'strict-nuclear'
                else:
                    return 'extended-nuclear'
            elif np.max(distances) > 1:
                return 'multi-generational'
        elif len(self.tree_map) == 1:
            return self.tree_map[eldest_node].status
        else:
            return 'misc'
    
    def find_singles(self, bari_id = 1): #if no bari id given use an arbitrary number
        single_males, single_females = [], []
        for individual in self.individuals.values():
            node = self.reverse_tree_map[individual.id]
            if individual.ismarried == False and individual.marriage_eligible == True and 'single' in self.tree_map[node].status:
                #individual is not currently married, wants to get married, and is not a widow(er)
                if individual.sex ==1:
                    single_males.append((bari_id, self.id, node, individual.id)) #hhid -> node -> iid
                else:
                    single_females.append((bari_id, self.id, node, individual.id))
        return single_males, single_females
                
    def find_marriage_differential(self):
        marriage_differential  = []
        for n in self.tree.traverse():
            node = self.tree_map[n.name]
            if node.status == 'married':
                i1,i2 = node.individuals
                if self.individuals[i1].sex == 1:
                    age_differential = self.individuals[i1].age - self.individuals[i2].age
                else:
                    age_differential = self.individuals[i2].age - self.individuals[i1].age
                marriage_differential.append(age_differential)
        return marriage_differential   
    
    def rename_tree(self):
        renamed_tree = copy.deepcopy(self.tree)
        for n in renamed_tree.traverse():
            if len(self.tree_map[n.name].individuals) > 0:
                rename = (':').join([str(iid) for iid in self.tree_map[n.name].individuals])
            else:
                rename = 'x'
            n.name = rename
        return renamed_tree
    
    def print_tree(self, node_idxes = False, status = False):
        if node_idxes == False:
            renamed_tree = self.rename_tree()
            print(renamed_tree.get_ascii(show_internal=True))
        
        else:
            print_tree(self.tree)
    
    def print_stats(self):
        self.print_tree()
        print('Household Active: ' + str(self.isactive))
        print('n_individuals: ' + str(self.n_individuals()))
        for node_id in self.tree_map:
            print(node_id, self.tree_map[node_id].__dict__, [self.individuals[iid].age for iid in self.tree_map[node_id].individuals.keys()])