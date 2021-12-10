import numpy as np
from collections import defaultdict
import itertools
import config
from Utility import Weibull_exponentiated, sigmoid, Transmission
from Infection import ImmunoInfection
import copy

class Individual:
    id_generator = itertools.count()
    def __init__(self, age, strain_type = 'S2'):
        self.id = next(Individual.id_generator) #iterator that generates a new unique id everytime it is called
        self.age = age
        self.isalive =True
        self.check_fertility()
        self.ismarried = False
        self.t_since_marriage = None
        self.marriage_eligible = False
        self.t_since_marriage_moveout = None
        self.marriage_age = []
        self.importation = False
        self.infection = ImmunoInfection(strain_type)
        self.is_transmitting = False
        self.index = False
        self.times_infected = 0
    
    @classmethod
    def override_id_generator(cls, start):
        Individual.id_generator = itertools.count(start)
    
    def check_alive(self): #check if person is still alive after specified time period
        if self.age > self.death_age:
            self.isalive = False
            
    def check_fertility(self):
        if self.age >=config.params.minimum_fertility_age and self.age < config.params.maximum_fertility_age: #assume males virility follow fertile rules in females
            self.isfertile = True
        else:
            self.isfertile = False

    def poly_4(self):
        a,b,c,d,e = config.params.poly_4_fertility_coefficients
        if self.age < config.params.minimum_fertility_age or self.age > config.params.maximum_fertility_age:
            return 0
        else:
            x = self.age
            return max(0, a *x**4 + b*x**3 + c*x**2 + d*x + e)
    
    def base_p_child(self):
        '''baseline probabilities of having a child-- P(Child | no_children)'''
        p_child = self.poly_4()*1.6 * config.params.historical_fertility_scale_fn(self.age) #1.375 is an additional scaling factor to match age pyramid stats
        return p_child 
    
    @classmethod
    def create_child(cls):
        if np.random.random() > config.params.p_male_child:
            return Male(0)
        else: 
            return Female(0)
        
    def check_move_out(self): 
        #based on Mohandpur Village A in 1991 "family structure and change in bangladesh" Amin 1998
        if self.t_since_marriage is not None and self.t_since_marriage_moveout is not None: #the latter occurs when they are the root
            if self.t_since_marriage >= self.t_since_marriage_moveout:
                return True
            else:
                return False
        else:
            return False
            
    def check_transmitting(self):
        self.is_transmitting = self.infection.is_shed
    
    
    def return_transmission(self):
        '''virus shed in a typical fecal-oral dose, returned as an ordered list. first idx = S1, 2nd = S2, 3rd = S3, 4th = WPV'''
        Tx = Transmission(self.id, self.times_infected, round(self.age,0), self.infection.viral_shed(self.age)) 
        return Tx
    
    def receive_transmission_dose(self, dose, modifier = 1, index=False):
        '''last step after all other updates are performed. This prevents the infection from being immediately infectious, tx =Transmission instance'''
        prechallenge_immunity = copy.deepcopy(self.infection.current_immunity)
        infection_success_flag = self.infection.attempt_infection(dose, modifier)   
        if infection_success_flag:
            self.times_infected += 1
        if index == True:
            self.index = True
        return (infection_success_flag, prechallenge_immunity)
            
    def update_transmissions(self, dt):
        if not self.infection.naiive:
            self.infection.update(dt) #dt = years
            self.check_transmitting()
        
class Male(Individual):
    def __init__(self, *args, **kwargs):
        Individual.__init__(self, *args, **kwargs)
        self.sex = 1
        self.minimum_marriage_age = config.params.m_marriage_age() 
        self.death_age = min(float(config.params.m_death_age_sampler(np.random.uniform(0,1))), 85)       
    
    def mortality_rate(self, minimum_mortality = 0.001):
        '''left over code, not removing becuase it's useful to keep
        used to determine the probability of death at a given age
        has been deprecated in favor of a death age sampler interpolated from simulated data using this function'''
        #the age thresholds for child_fn and elder_fn are guestimated from the data -- will need to check this when applying to other countries
        #minimum mortality is also guestimated from the data -- this is mostly to do away with a 0 mortality rate
        if self.age <=12:
            m = max(minimum_mortality, np.e**(-4.038012826363843 + -0.4620757529135555*self.age))
        elif self.age > 30: 
            m = max(minimum_mortality, np.e**(-9.57843000677779 + 0.08987439511446728*self.age))
        else:
            m = minimum_mortality
        return m
    
    def check_eligibility(self):
        #table 6.1... Data is from Matlab!
        #assuming normal distribution...
        if self.age >= self.minimum_marriage_age:
            self.marriage_eligible = True
            
    def update(self, dt):
        self.check_alive()
        if self.isalive:
            self.age += dt
            self.check_fertility()
            self.update_transmissions(dt)
            if self.ismarried == False and self.marriage_eligible == False:
                self.check_eligibility()
            elif self.ismarried and self.t_since_marriage is not None:
                self.t_since_marriage += 1
    
class Female(Individual):
    def __init__(self, *args, **kwargs):
        Individual.__init__(self, *args, **kwargs)
        self.n_children = 0 #keep track of n_children only on female side
        self.sex =0
        self.t_since_birth = None #period of time where she does not give birth again -- rounding it to 1 year-- this could probably be turned off tbh
        self.t_secondary_births =[]
        self.birth_age = []
        self.just_gave_birth = False #hack-y tracker to keep track of number of births per time unit for calibrating to Births/individual/year
        self.death_age = min(float(config.params.f_death_age_sampler(np.random.uniform(0,1))),85)
        self.minimum_marriage_age = config.params.f_marriage_age()
        self.children_stats = defaultdict(list)

    def mortality_rate(self, minimum_mortality = 0.001):
        '''left over code, not removing becuase it's useful to keep, see male version above'''
        if self.age <=12:
            m = max(minimum_mortality, np.e**(-4.154414473758413 + -0.4420609900131315*self.age)) #m = age specific mortality rate (deaths/individual/year)
        elif self.age > 30: 
            m = max(minimum_mortality, np.e**(-10.527925730119318 + 0.10104736306154925*self.age))
        else: 
            m= minimum_mortality
        return m    

    def birth(self, dt):
        '''determines whether a woman gives birth. 
        other than imposing a 1 year refractory period births are 
        births per year! To scale down to day, assume that Uniformly distributed across that time period'''
        if self.isfertile and self.ismarried:
            if self.t_since_birth != None:
                if self.t_since_birth < 1:
                    p_birth = 0
                else:
                    p_birth = np.max([0, self.base_p_child() * self.birth_interval() * self.birth_desire()]) # prevent it from going negative
            else:
                p_birth = self.base_p_child() #*self.birth_desire()
            
            if np.random.random() <= p_birth * dt:
                self.n_children +=1
                self.birth_age.append(self.age)
                if self.t_since_birth != None:
                    self.t_secondary_births.append((self.age, self.t_since_birth))
                self.t_since_birth = 0 
                self.just_gave_birth = True #flips on if birth happened
                return True
            else:
                self.just_gave_birth = False
                return False
        else: 
            self.just_gave_birth = False
            return False
    
    def birth_interval(self):
        '''scaling factor that modifies p_child given t years after the previous birth'''

        if self.age >= 15 and self.age <20:
            p_birth_given_time = Weibull_exponentiated.cdf(config.params.f_birth_interval_params['15-19'], self.t_since_birth)
        elif self.age < 30:
            p_birth_given_time = Weibull_exponentiated.cdf(config.params.f_birth_interval_params['20-29'], self.t_since_birth)
        elif self.age < 40:
            p_birth_given_time = Weibull_exponentiated.cdf(config.params.f_birth_interval_params['30-39'], self.t_since_birth)
        elif self.age < 50:
            p_birth_given_time = Weibull_exponentiated.cdf(config.params.f_birth_interval_params['40-49'], self.t_since_birth)
        else:
            p_birth_given_time = 0
        
        return p_birth_given_time
            
    def birth_desire(self):#obtained from Bangladesh DHSS 2014 Table 6.1
        '''scaling factor that modifies p_child given they have n living children already
        Based on the percent distribution of women aged 15-49 who still desire children given n living children'''
        p_desire = sigmoid(config.params.birth_desire_params, self.n_children) #this data is on still LIVING children...
        return p_desire
    
    def check_eligibility(self):
        if self.age >= self.minimum_marriage_age:
            self.marriage_eligible = True
    
    def update(self, dt):
        self.check_alive()
        if self.isalive:            
            if self.t_since_birth != None:
                self.t_since_birth += dt
            self.birth(dt)
            self.age += dt #advance time
            self.check_fertility() #advance time
            self.update_transmissions(dt)
            if self.ismarried == False and self.marriage_eligible == False: #once marriage eligible always marriage eligible
                self.check_eligibility()
            elif self.ismarried and self.t_since_marriage is not None:
                self.t_since_marriage += 1
                
def cdf_married(age, coefficients = [0.12459153286506645, 0.48875807794081044, 0.9875643033512792]):
    '''not used anymore, but left in because its nice to keep'''
    def growth_rate(x, p0,  r, k):
        '''population growth logistic equation: Verhulst equation'''
        return (k*p0*np.exp(r*x)) / (k + p0*(np.exp(r*x) - 1))

    x = age - 15 #shift the growth_rate equation for better calculation
    y = growth_rate(x, *coefficients)
    return y