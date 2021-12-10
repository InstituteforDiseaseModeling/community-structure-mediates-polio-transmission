from ete3 import Tree#, TreeStyle, NodeStyle
import scipy.stats
import numpy as np
import config
from collections import defaultdict
#Tree manipulation Functions
def remove_subtree(tree, internal_node, return_subtree = False):
    subtree = tree.search_nodes(name = internal_node)[0]

    removed = subtree.detach()
    
    if return_subtree:
        return tree, removed
    else:
        return tree

def add_child(tree, child_name, parental_node = None):
    if parental_node:
        subtree = tree.search_nodes(name=parental_node)[0]
    else:
        subtree = tree
    subtree.add_child(name=child_name)
    return tree

def print_tree(tree):
    print(tree.get_ascii(show_internal=True))

def retrieve_root_distance(tree, node):
    if node == tree.name:
        distance = 0
    else:
        distance = tree.get_distance(tree.name, node)
    return distance

def relabel_tree(tree):
    '''relabels the node names for subtrees so that 0 is always the root'''
    count =0
    for node in tree.traverse('levelorder'): #traverse down the tree row by row
        node.name = count
        count +=1
    tree.name = 0 #set the root to have id zero
    return tree

def find_siblings(tree, ref_node):
    node_generation = retrieve_root_distance(tree, ref_node)
    def conditional_function(target_node):
        tnode_generation = retrieve_root_distance(tree, target_node)
        if tnode_generation == node_generation:
            return True
        else:
            return False
    
    matches = [n.name for n in filter(conditional_function, tree.traverse())]

    return matches
    
def modified_counts(counter, max_range = 9):
    modified_counts = []
    for number in range(max_range):
        modified_counts.append(counter.pop(number,0))
    remainder = np.sum(list(counter.values()))
    modified_counts.append(remainder)
    return modified_counts   

def find_bin_idx(x, bins):
    '''find the bin index for the number x'''
    return np.digitize(x, bins)

class Distribution:
    '''wrapper class for various distributions to optimize'''
    def __init__(self, pdf, cdf, integral=None):
        self.pdf = pdf #or pmf
        self.cdf = cdf
        self.integral = integral  
        
def weibull_pdf(x, params_array):
    rate, k = params_array[0], params_array[1]
    pdf = k/rate * (x/rate)**(k-1) * np.e**-(x/rate)**k
    return pdf

def weibull_cdf(x, params_array):
    rate, k = params_array[0], params_array[1]
    cdf = 1 - np.e**-(x/rate) **k
    return cdf
#exponentiated weibull distributions -- 3 parameter distribution, when alpha =1 the same as a weibull
def weibull_exponentiated_distribution(params_array,x):
    rate, k, alpha = params_array[0], params_array[1], params_array[2]
    return alpha * (k/rate)*(x/rate)**(k-1)*(1- np.e**-(x/rate)**k)**(alpha -1) * np.e**(-(x/rate)**k)

def weibull_exponentiated_cdf(params_array, x):
    rate, k, alpha = params_array[0], params_array[1], params_array[2]
    return (1 - np.e **-(x/rate) **k)**alpha

def weibull_exponentiated_integral(params_array, lower, upper):
    integral = weibull_exponentiated_cdf(params_array, upper) - weibull_exponentiated_cdf(params_array, lower)
    return integral

#x is number of events/successes that occur

def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y
    
Weibull_exponentiated = Distribution(weibull_exponentiated_distribution, weibull_exponentiated_cdf, weibull_exponentiated_integral) #params_array = rate, k, alpha

        
class Transmission:
    '''helper class to organize transmissions: serves as a more readable dictionary'''
    def __init__(self, source_iid, source_n_infected, source_age, dose):
        self.dose = dose
        self.source_iid = source_iid
        self.source_age = int(min(round(source_age,0), 80))
        self.hhid = None
        self.bari_id = None
        self.tx_history = defaultdict(list)
        self.secondary_cases = 0
        self.source_n_infected = source_n_infected #times the source infection was infected
    
    def return_source_stats(self):
        return [self.bari_id, self.hhid, self.source_iid, self.source_age]
        
    def optimal_contact_age_idx(self, n_age_bins, p, mask):
        nonzero_age_idxes = np.asarray([_ for _ in range(n_age_bins+1)])[mask]
        relative_probabilities = np.asarray(p)[mask] / np.sum(np.asarray(p)[mask])
        return np.random.choice(nonzero_age_idxes, p= relative_probabilities)
        
class Demographics:
    '''helper class to keep track of births and deaths -- only used at inter-Village level because we can safely assume that only births and deaths affect village-level census'''
    def __init__(self,Individual, type): #type = 0 = 'death', 1 = 'birth'
        self.id = Individual.id
        self.source_age = Individual.age
        self.individual = Individual
        self.type = type
        self.village_id = None
        
        