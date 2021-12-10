import numpy as np
from collections import defaultdict
import scipy
import scipy.interpolate as interpolate
import Utility
import json

class Params:
    def __init__(self, **kwargs):
        #self.set_start_conditions()
        
        self.n_initial_villages = 20 #used only if not using the preloaded taniuchi stats
        self.n_initial_baris = 2000 #used only if not using the preloaded taniuchi stats
        self.sim_expansion_years = 140
        self.sim_contraction_years = 40
        self.n_initial_hh_per_bari = 1
    
        #only affects the unrooting procedure: if none then unrooting occurs when root completely empties
        self.preferred_sex = 1 #1 for male, 0 for female:modify at your own risk... if you don't want sex bias to affect unrooting procedure set to None
        
        #INDIVIDUAL CLASS
        self.minimum_fertility_age = 15
        self.maximum_fertility_age = 50
        self.elder_age = 60
        self.infant_age = 1
        self.p_male_child = 0.5 #allow for biased birth probabilities
        self.f_marriage_age = lambda: np.random.normal(19.3, 3.7) #Matlab data --Table 6.1
        self.m_marriage_age = lambda: np.random.normal(27.3, 4.5) #Matlab data --Table 6.1
        self.historical_fertility_scale_fn = lambda x: 1
        self.historical_fertility_multiplier = 3
        self.strain_type = 'S2'
        
        #sample death ages: values are the cumsum and bin edges of an empirically derived simulated death age distribution
        #simulated 10000 Males and Females and recorded their death ages
        self.historical_m_death_age_sampler = interpolate.interp1d([0.0, 0.068, 0.13, 0.17300000000000001, 0.22000000000000003, 0.276, 0.319, 0.375, 0.431, 0.458, 0.495, 0.523, 0.556, 0.6080000000000001, 0.6990000000000001, 0.769, 0.861, 0.94, 0.978, 0.996, 1.0] , [0.0, 4.65, 9.3, 13.950000000000001, 18.6, 23.25, 27.900000000000002, 32.550000000000004, 37.2, 41.85, 46.5, 51.150000000000006, 55.800000000000004, 60.45, 65.10000000000001, 69.75, 74.4, 79.05000000000001, 83.7, 88.35000000000001, 93.0])
                
        self.historical_f_death_age_sampler = interpolate.interp1d([0.0, 0.07400000000000001, 0.14800000000000002, 0.199, 0.25, 0.301, 0.356, 0.39199999999999996, 0.44399999999999995, 0.477, 0.514, 0.552, 0.5910000000000001, 0.6350000000000001, 0.6860000000000002, 0.7620000000000001, 0.8250000000000002, 0.9030000000000001, 0.9620000000000002, 0.9880000000000002, 1.0000000000000002] , [0.0, 4.7, 9.4, 14.100000000000001, 18.8, 23.5, 28.200000000000003, 32.9, 37.6, 42.300000000000004, 47.0, 51.7, 56.400000000000006, 61.1, 65.8, 70.5, 75.2, 79.9, 84.60000000000001, 89.3, 94.0])
        
        self.m_death_age_sampler_2015 = interpolate.interp1d([0.    , 0.0358, 0.0452, 0.0469, 0.0497, 0.0521, 0.0546, 0.0576,
            0.0599, 0.0626, 0.0647, 0.0682, 0.0705, 0.0726, 0.0775, 0.0814,
            0.0889, 0.0971, 0.1039, 0.1192, 0.1322, 0.1509, 0.1764, 0.1975,
            0.2388, 0.2707, 0.3294, 0.392 , 0.4402, 0.526 , 0.5845, 0.6772,
            0.771 , 0.8264, 0.8964, 0.9325, 0.9685, 0.9884, 0.9948, 0.9989,
            1.    ], 
            [  0. ,   2.6,   5.2,   7.8,  10.4,  13. ,  15.6,  18.2,  20.8,
            23.4,  26. ,  28.6,  31.2,  33.8,  36.4,  39. ,  41.6,  44.2,
            46.8,  49.4,  52. ,  54.6,  57.2,  59.8,  62.4,  65. ,  67.6,
            70.2,  72.8,  75.4,  78. ,  80.6,  83.2,  85.8,  88.4,  91. ,
            93.6,  96.2,  98.8, 101.4, 104. ])

        self.f_death_age_sampler_2015 = interpolate.interp1d(
            [0.    , 0.033 , 0.0416, 0.0433, 0.0457, 0.0487, 0.0511, 0.0553,
            0.0583, 0.0603, 0.0631, 0.0655, 0.0677, 0.0712, 0.0747, 0.0774,
            0.0809, 0.0872, 0.0941, 0.1052, 0.1148, 0.1335, 0.1545, 0.1724,
            0.2076, 0.2499, 0.2861, 0.3509, 0.4291, 0.487 , 0.5825, 0.6872,
            0.7596, 0.8549, 0.9253, 0.9584, 0.9862, 0.9969, 0.9989, 0.9997,
            1. ],
            [  0.  ,   2.65,   5.3 ,   7.95,  10.6 ,  13.25,  15.9 ,  18.55,
             21.2 ,  23.85,  26.5 ,  29.15,  31.8 ,  34.45,  37.1 ,  39.75,
             42.4 ,  45.05,  47.7 ,  50.35,  53.  ,  55.65,  58.3 ,  60.95,
             63.6 ,  66.25,  68.9 ,  71.55,  74.2 ,  76.85,  79.5 ,  82.15,
             84.8 ,  87.45,  90.1 ,  92.75,  95.4 ,  98.05, 100.7 , 103.35,
            106.  ])
        
        self.m_death_age_sampler = self.m_death_age_sampler_2015
        self.f_death_age_sampler = self.f_death_age_sampler_2015
 
        
        
        #exponentiated weibull parameters for birth intervals
        self.f_birth_interval_params =  { '15-19':  [7.008200010283388e-07, 0.18519982325420162, 4955315.919483488],
                                         '20-29': [0.617506311482934, 0.6131814950322227, 15.055481013015964],
                                         '30-39': [0.004432998576956875, 0.25617711748578326, 327.0658166547044],
                                         '40-49': [8.433674445639343e-10, 0.10913540296074727, 110510.2718509799]}
        
        #coefficients for the 4-D age specific base fertility function
        self.poly_4_fertility_coefficients = [-1.25493395e-06, 1.79042907e-04,-9.63000029e-03, 2.09000305e-01, -1.19246824e+00]
        
        #birth desire coefficients
        self.birth_desire_params = [ 1.36874821,  4.80536379e-01,  1-4.80536379e-01, -4.101415]        
        
        
        #HOUSEHOLD CLASS
        self.p_stem_hh = 1.0
        
        self.post_marriage_move_out= 0# [1.5773640885, 1.8171503040000005] #rate, scale parameters of a weibull distribution
        #if set as None, then assume never move out, if 0 assume immediate moveout
        if self.post_marriage_move_out != None and self.post_marriage_move_out != 0:
            rate, scale = self.post_marriage_move_out
            self.hh_move_out_fn = lambda: np.random.weibull(scale) * rate
        elif self.post_marriage_move_out == None:
            self.hh_move_out_fn = lambda: 1e8 #arbitrarily big number to prevent moveout
        else:
            self.hh_move_out_fn = lambda: 0 #immediately moveout after marriage
        
        self.__dict__.update(kwargs) #overwrite  
        
        self.file_path = 'sims/simulation_1.pkl'
        
        #age assortivity matrix stats
        self.n_age_bins = 15
        self.n_age_bin_step = 5
        self.age_bins = [_ for _ in range(self.n_age_bin_step, self.n_age_bin_step * (self.n_age_bins +1), 5)] + [np.inf]
        
        self.age_assortivity = True  
        self.create_assortivity_matrix()
        
        self.fecal_oral_dose = 5e-7
        self.beta_hh = 10
        self.beta_bari = 5
        self.beta_village = 5
        self.beta_inter_village = 5
        
        
        self.single_parameter = False
        self.beta_global = 0
        
        self.vaccination_campaign = 'bOPV'
        self.hh_transmission = False
        self.bari_transmission = False
        self.village_transmission = False
        
        
        self.infant_initial_Nab = {'tOPV': 200, 'bOPV': 31}
        self.initial_Nab_under_5 = 580
     
    def create_assortivity_matrix(self):
        self.village_age_assortivity_matrix = np.zeros((81,81))
        self.bari_hh_age_assortivity_matrix = np.zeros((81,81))
        for c in range(0,81):
            for e in range(0,81):
                self.village_age_assortivity_matrix[e][c] = Age_assortivity_MM.mixture_pdf_fn(c, e, self.age_assortivity) #rows = ego age, columns = contact age
                self.bari_hh_age_assortivity_matrix[e][c] = HH_Age_Assortivity_MM.mixture_pdf_fn(c,e,self.age_assortivity)  
        
        
    def overwrite(self, *json_dict, **kwargs):
        for params in json_dict:
            for key in params:
                setattr(self, key, params[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])
        self.create_assortivity_matrix()

    def set_start_conditions(self, json_file='Taniuchi_village_start.json'):
        self.start_conditions = json.load(open(json_file))
        self.n_initial_villages = np.sum([x[1] for x in self.start_conditions.values()])
        
class Age_assortivity_MM:
    '''between bari age assortivity matrix'''
    def __init__(self):
        pass
    
    @classmethod
    def stdev(cls, x, v=3.28,d=0.45,b=-4.98):
        return v*(x<20) + (x>=20)*(d*x + b)
    
    @classmethod
    def weight(cls, x, e=-0.0017,f=0.066,g=-0.021):
        return max((x<=40) * (e*x**2 + f*x + g) + (x>40)*0, 0)
    
    @classmethod
    def mixture_pdf_fn(cls, contact_age, ego_age, age_assortivity_flag = True):
        #mu = GMM.mu(ego_age, *coefficients_mu)
        if age_assortivity_flag:
            stdev = cls.stdev(ego_age)
            w = cls.weight(ego_age)
            probability =  w * scipy.stats.norm.pdf(contact_age, loc=ego_age,scale=stdev) + (1-w)/80
        else:
            probability = 1/80.
        return probability
    
    @classmethod
    def bin_pdf(cls, age1, age2, ego_age):
        p = cls.mixture_cdf_fn(age2, ego_age) - cls.mixture_cdf_fn(age1, ego_age)
        return p
    
    @classmethod
    def mixture_cdf_fn(cls, contact_age, ego_age):
        stdev = cls.stdev(ego_age)
        w = cls.weight(ego_age)
        p = (w) * scipy.stats.norm.cdf(contact_age, loc = ego_age, scale = stdev) + (1-w)/80 * contact_age
        return p
    
    @classmethod
    def age_bin_pdf(cls, ego_age, age_bins):
        age_weights = []
        for bin in age_bins:
            if bin != np.inf:
                right = bin
                left = bin - 5
                age_weights.append(cls.bin_pdf(left, right, ego_age))
            else:
                age_weights.append(1-cls.mixture_cdf_fn(age_bins[-2], ego_age))
                
        return np.asarray(age_weights) / np.sum(age_weights) #empirical rounding errors
        
class HH_Age_Assortivity_MM:
    '''within hh age assortivity matrix-- Three generation Gaussian Mixture Model
    Each gaussian corresponds to a different generation (ages 0-20, 20-50, 50+)'''
    def __init__(self):
        pass
    
    @classmethod
    def concave_fn(cls,x, a,b,c):
        return a*(x - b)**2 + c 
    
    @classmethod
    def stdevs(cls, x):
        if x <=20:
            m_i, b_i = [0.26960995, 2.70785968] #infant gaussian
            m_a, b_a = [0.25812943, 7.4845975 ] #adult gaussian
            m_e, b_e = [-0.01787006,  9.7365903] #elder gaussian
            
        elif x <=50:
            m_i, b_i = [ 0.18893171, -2.31804373]
            m_a, b_a = [0.25812943, 7.4845975 ]
            m_e, b_e = [-0.37965034, 22.96636377]
        else:
            m_i, b_i = [-0.01596075,  7.73650253]
            m_a, b_a = [-1.08231331e-02,  1.75959036e+01]
            m_e, b_e = [ 0.17611273, -0.8546759 ]
        
        y_i = m_i*x + b_i  #stdev of the infant gaussian
        y_a =  m_a*x + b_a #stdev of the adult gaussian
        y_e = m_e*x + b_e #stdev of the elder gaussian
        
        return np.asarray([y_i, y_a, y_e])
        
    @classmethod
    def means(cls, x):
        if x <= 20:
            m_infant, m_older, b_i, b_a, b_e = [ 0.66544271,  0.66544271, 3.24819905, 30.8177947 , 58.36860468]
            y_i = m_infant*x + b_i
            y_a = m_older*x + b_a
            y_e = m_older*x + b_e
        elif x <=50:
            m_infant, m_older, b_i, b_a, b_e = [ 0.48622295,  1.00557581, -9.22902505, -1.46407176, 28.30562542]
            y_i = m_infant*x + b_i
            y_a = m_older*x + b_a
            y_e = m_older*x + b_e
        else:
            infant_concave_params = [3.59379721e-02, 6.70357768e+01, 8.23084115e+00]
            adult_concave_params = [3.13007006e-02, 6.10738322e+01, 3.55015021e+01]
            elder_concave_params = [-4.46633477e-05, -1.60549161e+03,  1.92760260e+02]

            y_i = cls.concave_fn(x, *infant_concave_params)
            y_a = cls.concave_fn(x, *adult_concave_params)
            y_e = cls.concave_fn(x, *elder_concave_params)
        return np.asarray([y_i, y_a, y_e])
        
    @classmethod
    def weights(cls, x):
        if x <= 20:
            m_i, b_i = [-0.01291969,  0.51083052]
        else:
            m_i, b_i = [-0.00200305,  0.4674363 ]

        y_e = m_i*x + b_i
        y_a = cls.concave_fn(x, 7.37508443e-05, 4.46613586e+01, 3.05877495e-01)
        y_i = 1 - y_a - y_e
        return np.asarray([y_i, y_a, y_e])
        
    @classmethod
    def mixture_pdf_fn(cls, contact_age, ego_age, age_assortivity_flag = True):
        
        if age_assortivity_flag:
            gaussian_means = cls.means(ego_age)
            gaussian_stdevs = cls.stdevs(ego_age)
            gaussian_weights = cls.weights(ego_age)
            probability = 0
            for i in range(3):
                probability += gaussian_weights[i] * scipy.stats.norm.pdf(contact_age, loc = gaussian_means[i], scale = gaussian_stdevs[i])
        else: #uniform
            probability = 1/80.
        return probability
    
    @classmethod
    def mixture_cdf_fn(cls, contact_age, ego_age):
        gaussian_means = cls.means(ego_age)
        gaussian_stdevs = cls.stdevs(ego_age) 
        gaussian_weights = cls.weights(ego_age)
        probability = 0
        for i in range(3):
            probability += gaussian_weights[i] * scipy.stats.norm.cdf(contact_age, loc = gaussian_means[i], scale = gaussian_stdevs[i])
        return probability
        
    @classmethod
    def bin_pdf(cls, age1, age2, ego_age):
        p = cls.mixture_cdf_fn(age2, ego_age) - cls.mixture_cdf_fn(age1, ego_age)
        return p
        
    @classmethod
    def age_bin_pdf(cls, ego_age, age_bins):
        age_weights = []
        for bin in age_bins:
            if bin != np.inf:
                right = bin
                left = bin - 5
                age_weights.append(cls.bin_pdf(left, right, ego_age))
            else:
                age_weights.append(1-cls.mixture_cdf_fn(age_bins[-2], ego_age))
                
        return np.asarray(age_weights) / np.sum(age_weights) #empirical rounding errors
                
params = Params()



