# modular_ion_current_model.py

'''
motivation: 
there are many conductance based models (and their reduction, see ch.5 of Izhikevich)
They can all be simply constructed by combining various ion-currents, 
as well as respective gating variables

lets setup a modular system where specific models (E.g. I_Na,p+I_K model)
can be built from subcomponent ion-currents with parameters specified in some structure (.csv)


notes:
it is known that any reduced model with atleast 1 amplifying and 1 resonant gating variable
will yield a model capable of spiking behaviour. we can tag preset currents with amp/res type
and evaluate whether the modular model could have spiking.
this would be based on relative timeconstants, or if m(V)->m_inf(V) assumption is made for atelast 1 current


'''
# imports
import numpy as np
import pandas as pd

from neuron_functions import *

# global variables
temp = 20+273.15 #Kelvin

# example nernst potentials:
# these are not particular to any conductance model
E_Na = 54.9
E_K = -77.5

# pseudo-code

# where ion conductance/current is expressed as: I = [g*n*h*(V-E)]_i

# parameter-table headers
# current-name | Nernst Potential | g_max | exp (amp)^x | exp (res)^y 
# | other associated gate parameters... (null if no amp/res gate, or m_inf(V) assumption)
# | V_1/2 | k | Vmax | sigma | Camp | Cbase

# parameter data from chapter.2 of Izhikevich
ion_currents_data = pd.read_csv('ion_current_params.csv',index_col=0)

# example: (table data can be defaults, and then edited)
# Na_t |  ENa = nernst_full(T,1,50,440) | 120 mS/cm^2 | 4 | 1 | 
# | -40| 15 | -38 | 30 | 0.46 | 0.04 | -62|-7|-67|20|7.4|1.2

# leak current | ECl = nernst_full(T,-1,65,560) | 0.3 | 0 | 0 
# (0-exp means no gating variables, dismiss further calc by if-statement)

# a dedicated function under a generic class is created with these parameters
# the input to the class are the parameters (g_max, exponent, E_i for gn(V-Ei) format)
# the input to the function is simply membrane potential (V) & previous state variables
# or ion concentration (E.g. [Ca2+] for ion-gated currents (how to track ion conc.?))

# see figure 2.20 of Izhikevich
class volt_gate_model:
    # I=g(A^a)(B^b)(V-E)
    # where A-gate is the activation variable
    # and   B-gate is the inactivation variable
    # reduce A or B = True is to make A(V) = A_inf(V) assumption for reduced model, 
    # otherwise timestep change in gateA/B is included 
    def __init__(self,ion_currents_data,ion_current_name: str, nernst:float,reduceA: bool = False, reduceB: bool = False):
        self.nernst = nernst
        self.p = ion_currents_data.loc[ion_current_name].to_dict()
        
        self.reduceA = reduceA
        self.reduceB = reduceB
        self.expA = self.p["expA"]
        self.expB = self.p["expB"]
                

    def __call__(self, V,A,B):

        # produce current equation
            # pre-calc timestep gate equation x2 if reduceA/B=False
            # and use as A/B variable in current equation I=g(A^a)(B^b)(V-E)
            # where A and B have been A+dA/dt stepped-forward

        if self.expA != 0:
            if self.reduceA == False:
                self.tau_A = timeconstant(V,self.p["CbaseA"],self.p["CampA"],self.p["VmaxA"],self.p["sigmaA"])
                self.A_inf = gate_inf(V, self.p["VhalfA"],self.p["kA"])
                self.Astep = A + (self.A_inf-A)/self.tau_A
            else:
                self.Astep = self.A_inf
        elif self.expA == 0:
            self.Astep = 1 # to avoid calculating tau or inf with empty parameters from .csv params

        if self.expB != 0:
            if self.reduceB == False:
                self.tau_B = timeconstant(V,self.p["CbaseB"],self.p["CampB"],self.p["VmaxB"],self.p["sigmaB"])
                self.Bstep = B + (self.B_inf-B)/self.tau_B
                self.B_inf = gate_inf(V, self.p["VhalfB"],self.p["kB"])
            else:
                self.Bstep = self.B_inf
        elif self.expB == 0:
            self.Bstep = 1

        # current dV/dt fragment
        # I = -g(A^a)(B^b)(V-E)
        return -self.p["g_max"]*(self.Astep**self.expA)*(self.Bstep**self.expB)*(V-self.nernst)
         

    def gate_inf(self,V,Vhalf,k):
        return  1/(1+np.exp((Vhalf-V)/k))   
    
    def timeconstant(self,V,Cbase,Camp,Vmax,sigma):
        return Cbase + Camp*np.exp((-(Vmax-V)**2)/sigma**2)  
    

# check: sample calc timeconstants of each gate.
# later evaluate timeconstant ratio for spiking potential 

# building the differential equation:
# C*dV/dt =  I + fn1 + fn2 + ...
# dn_i/dt = (gate_inf(V,gate_id)-gate_0)/timeconstant(V,gate_id) 

# ... numerically run the model.

# toggle for reducing degrees of freedom, most commonly by removing gate for gate_inf(V)

# now that we have a script for taking model parameters and constructing them to be simulated
# we now build tools to analyze the simulations
# this has already been done very simply using graphs to present the simulation results
# but we also want to analyse phase portraits, bifurcations, equilibriums, eigenvalues, frequency, etc.