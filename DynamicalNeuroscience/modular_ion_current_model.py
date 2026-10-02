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
from neuron_functions import *
# pseudo-code

# where ion conductance/current is expressed as: I = [g*n*h*(V-E)]_i

# parameter-table headers
# current-name | Nernst Potential | g_max | exp (activator)^x | exp (in-activator)^y 
# | other associated gate parameters... (null if no (in)-activator gate, or m_inf(V) assumption)
# | V_1/2 | k | Vmax | sigma | Camp | Cbase

# example: (table data can be defaults, and then edited)
# Na_t |  ENa = nernst_full(T,1,50,440) | 120 mS/cm^2 | 4 | 1 | 
# | -40| 15 | -38 | 30 | 0.46 | 0.04 | -62|-7|-67|20|7.4|1.2

# leak current | ECl = nernst_full(T,-1,65,560) | 0.3 | 0 | 0 
# (0-exp means no gating variables, dismiss further calc by if-statement)

# a dedicated function under a generic class is created with these parameters
# the input to the class are the parameters (g_max, exponent, E_i for gn(V-Ei) format)
# the input to the function is simply membrane potential (V) & previous state variables
# or ion concentration (E.g. [Ca2+] for ion-gated currents (how to track ion conc.?))

# for every exp >0: make dgate/dt function (most the parameters are for this gating variable diff.)

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