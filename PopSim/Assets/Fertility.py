import numpy as np
from collections import Counter, defaultdict
import itertools
import scipy.stats
from scipy.special import gamma as gamma_fn
from scipy.special import loggamma as loggamma_fn
from scipy.special import gammainc as lower_inc_gamm_fn
from scipy.special import erf

class Distribution:
    '''wrapper class for various distributions to optimize'''
    def __init__(self, pdf, cdf, integral=None):
        self.pdf = pdf #or pmf
        self.cdf = cdf
        self.integral = integral  
		
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
Weibull_exponentiated = Distribution(weibull_exponentiated_distribution, weibull_exponentiated_cdf, weibull_exponentiated_integral) #params_array = rate, k, alpha

def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def multinomial(p_array, x_array):
    n_factorial = gamma_fn(np.sum(x_array) + 1)
    x_i_factorial_product_sum = 1
    for x_i in x_array:
        x_i_factorial = gamma_fn(x_i + 1)
        x_i_factorial_product_sum *= x_i_factorial
    p_i_product_sum = 1
    for x_i, p_i in zip(x_array, p_array):
        p_i_product_sum *= (p_i**x_i)
    pmf = (n_factorial / x_i_factorial_product_sum) * p_i_product_sum
    return pmf
    
class Individual:
    id_generator = itertools.count()
    scale = 1.0
    def __init__(self, age, base_fertility_params):
        self.id = next(Individual.id_generator) #iterator that generates a new unique id everytime it is called
        self.age = age
        self.base_fertility_params = base_fertility_params
        self.isalive =True
        self.check_fertility()
        
    def check_alive(self, dt): #check if person is still alive after specified time period
        p_death = 1-np.e**(-self.mortality_rate()*dt)
        if np.random.random() < p_death:
            self.isalive = False
            self.death_age = self.age
            
    def check_fertility(self):
        if self.age >=15 and self.age <50: #assume males virility follow fertile rules in females
            self.isfertile = True
        else:
            self.isfertile = False

    def poly_4(self, a = -1.25493395e-06,  b = 1.79042907e-04, c=-9.63000029e-03,  d=2.09000305e-01, e =  -1.19246824e+00):
        if self.age < 15 or self.age > 50:
            return 0
        else:
            x = self.age
            return max(0, a *x**4 + b*x**3 + c*x**2 + d*x + e)
    
    def base_p_child(self, dt):
        '''baseline probabilities of having a child-- P(Child | no_children)'''
        #p_child =1-np.e**(-f_rate*dt)
        p_child = self.poly_4()
        #p_child = self.double_exp_fertility() #self.polynomial_fertility()
        return p_child 
    
    @classmethod
    def create_child(cls):
        if np.random.random() > 0.5:
            return Male(0)
        else: 
            return Female(0)

class Male(Individual):
    def __init__(self, *args, **kwargs):
        Individual.__init__(self, *args, **kwargs)
        self.sex = 1
    
    def mortality_rate(self, minimum_mortality = 0.001):
        #the age thresholds for child_fn and elder_fn are guestimated from the data -- will need to check this when applying to other countries
        #minimum mortality is also guestimated from the data -- this is mostly to do away with a 0 mortality rate
        if self.age <=12:
            m = max(minimum_mortality, np.e**(-4.038012826363843 + -0.4620757529135555*self.age))
        elif self.age > 30:
            m = max(minimum_mortality, np.e**(-9.57843000677779 + 0.08987439511446728*self.age))
        else:
            m = minimum_mortality
        return m
    
    def update(self, dt):
        self.check_alive(dt)
        if self.isalive:
            self.age += dt
            self.check_fertility()
    
class Female(Individual):
    def __init__(self, *args, **kwargs):
        Individual.__init__(self, *args, **kwargs)
        self.n_children = 0 #keep track of n_children only on female side
        self.sex =0
        self.t_since_birth = None #period of time where she does not give birth again -- rounding it to 1 year-- this could probably be turned off tbh
        self.t_secondary_births =[]
        self.birth_age = []
        self.p_birth_array = [] #temporary container
        self.ismarried=False
        
        self.just_gave_birth = False #hack-y tracker to keep track of number of births per time unit for calibrating to Births/individual/year

    def mortality_rate(self, minimum_mortality = 0.001):
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
        treated as independent of one another.. may want to change this'''
        if self.isfertile and self.ismarried:
            if self.t_since_birth != None:
                p_birth = np.max([0, self.base_p_child(dt) * self.birth_interval() * self.birth_desire()]) # prevent it from going negative
                
            else:
                p_birth = self.base_p_child(dt) #*self.birth_desire()
            self.p_birth_array.append(p_birth)
            
            if np.random.random() <= p_birth:
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
            self.p_birth_array.append(0)
            self.just_gave_birth = False
            return False
    
    def birth_interval(self):
        '''scaling factor that modifies p_child given t years after the previous birth'''
        weibull_exponentiated_params = { '15-19':  [7.008200010283388e-07, 0.18519982325420162, 4955315.919483488],
                                         '20-29': [0.617506311482934, 0.6131814950322227, 15.055481013015964],
                                         '30-39': [0.004432998576956875, 0.25617711748578326, 327.0658166547044],
                                         '40-49': [8.433674445639343e-10, 0.10913540296074727, 110510.2718509799]}
        if self.age >= 15 and self.age <20:
            p_birth_given_time = Weibull_exponentiated.cdf(weibull_exponentiated_params['15-19'], self.t_since_birth)
        elif self.age < 30:
            p_birth_given_time = Weibull_exponentiated.cdf(weibull_exponentiated_params['20-29'], self.t_since_birth)
        elif self.age < 40:
            p_birth_given_time = Weibull_exponentiated.cdf(weibull_exponentiated_params['30-39'], self.t_since_birth)
        elif self.age < 50:
            p_birth_given_time = Weibull_exponentiated.cdf(weibull_exponentiated_params['40-49'], self.t_since_birth)
        else:
            p_birth_given_time = 0
        
        return p_birth_given_time
            
    def birth_desire(self):#obtained from Bangladesh DHSS 2014 Table 6.1
        '''scaling factor that modifies p_child given they have n living children already
        Based on the percent distribution of women aged 15-49 who still desire children given n living children'''
        floor = 4.80536379e-01 
        params = [ 1.36874821,  floor,  1-floor, -4.101415]
        p_desire = sigmoid(params, self.n_children) #this data is on still LIVING children...
        return p_desire

    def update(self, dt):
        self.check_alive(dt)
        if self.isalive:            
            if self.t_since_birth != None:
                self.t_since_birth += dt
            self.birth(dt)
            self.age += dt #advance time
            self.check_fertility() #advance time            