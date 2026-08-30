# hodgkin huxley model.py

'''
Based on chapter 2 of Dynamical Systems in Neuroscience by Eugene M Izhikevich

the code bits of neuronModelNotebook.ipynb were adapted from this script
'''

# for squid axon ahh neuron 

import math
import matplotlib.pyplot as plt
import numpy as np
from nernstEquation import nernst_full
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import signal


# constants

#nernst potential data (exercise 1)
T=20+273.15 #Kelvin
Capa = 1 #mF/cm^2

def nernst_full(T,z,IonIn,IonOut): # K,,Mol,Mol
    R=8315 #mJ/(K*Mol)
    F=96480 #C/Mol
    E_ion=R*T/z/F*math.log(IonOut/IonIn) #mV
    return E_ion

# nernst potentials
# using ch.2 ex.1 data
EK = nernst_full(T,1,430,20)
ENa = nernst_full(T,1,50,440)
ECl = nernst_full(T,-1,65,560)


gKmax = 36 #mS/cm^2
gNamax = 120 #mS/cm^2
gL = 0.3 #mS/cm^2


# taken from 'parameters figure 2.20'
gate_m = {
    "Vhalf":-40,
    "k":15/2, 
    "Vmax":-38,
    "sigma":30,
    "Camp":0.46,
    "Cbase":0.04,
}
gate_h = {
    "Vhalf":-62,
    "k":-7,
    "Vmax":-67,
    "sigma":20,
    "Camp":7.4,
    "Cbase":1.2,
}
gate_n = {
    "Vhalf":-53,
    "k":15,
    "Vmax":-79,
    "sigma":50,
    "Camp":4.7,
    "Cbase":1.1,
}


#eq(2.11) Boltzmann equation
def gate_inf(V,gate):
    return  1/(1+np.exp((gate["Vhalf"]-V)/gate["k"])) 

#eq(2.12) gaussian
def timeconstant(V,gate):
    return gate["Cbase"] + gate["Camp"]*np.exp((-(gate["Vmax"]-V)**2)/gate["sigma"]**2) 

#eq(2.9) [gate] activation variable dynamics dn/dt
def gate_diff(V,gate_id,gate_0):
    return (gate_inf(V,gate_id)-gate_0)/timeconstant(V,gate_id) 

# current pulser (assumes same units as time series)
pulse_mag = 15 #mA
pulse_period = 10
pulse_length = 1
pulse_delay = 5 

# take a time-axis array and generate a PWM-style pulse array.
# pulses from zero to given magnitude
def current_pulser(magnitude,period,length,delay,timeseries):
    frequency = 1/period
    pwm_array = (signal.square(frequency*2*np.pi*t,duty=length/period)+1)/2*magnitude
    return np.where(timeseries >= delay, pwm_array, 0)


dt = 1/100 #ms
duration = 40 #ms

# initializing arrays
t = np.arange(0,duration+1,dt)
range = len(t)+1 #+1 accounts for i+1 operation within for loop
V = np.zeros(range)
m = np.zeros(range)
h = np.zeros(range)
n = np.zeros(range)
I_K = np.zeros(range)
I_Na = np.zeros(range)
I_L = np.zeros(range)


# Determining the resting membrane potential V_rest:

# Varray = np.linspace(-40,-90,40)
# Itest = np.zeros(len(Varray))
# for i,Vstep in enumerate(Varray):
#     Itest[i] = -(gKmax*(gate_inf(Vstep,gate_n)**4)*(Vstep-EK))-(gNamax*(gate_inf(Vstep,gate_m)**3)*(gate_inf(Vstep,gate_h))*(Vstep-ENa))-(gLmax*(Vstep-ECl))
    
# plt.figure()
# plt.plot(Varray,Itest, label="current flows")
# plt.axhline(y=0,color='r',linestyle='--')
# plt.axvline(x=-65,color='r',linestyle='--')
# plt.show()

V_0 = -65.5

# intial values of iterated variables
V[0] = V_0
m[0] = m_0 = gate_inf(V_0,gate_m)
h[0] = h_0 = gate_inf(V_0,gate_h)
n[0] = n_0 = gate_inf(V_0,gate_n)

I_K[0]=gKmax*(n_0**4)*(V[0]-EK)
I_Na[0]=gNamax*(m_0**3)*h_0*(V[0]-ENa)
I_L[0]=gL*(V[0]-ECl)

Iapplied = np.zeros(range)


# Simulating neuron membrane potential - iterating differential equation: 
for i,time in enumerate(t):

    # programming applied current pulses
    if time>=2   and time<=2.5:
        Iapplied[i]=15
    elif time>=10 and time<=10.5:
        Iapplied[i]=40

    # guiding equation:
    dV = (Iapplied[i]-I_K[i]-I_Na[i]-I_L[i])/Capa*dt
    V[i+1] = V[i]+dV

    dm = gate_diff(V[i],gate_m,m[i])*dt
    dh = gate_diff(V[i],gate_h,h[i])*dt
    dn = gate_diff(V[i],gate_n,n[i])*dt

    m[i+1]=m[i]+dm
    h[i+1]=h[i]+dh
    n[i+1]=n[i]+dn

    I_K[i+1]=gKmax*(n[i+1]**4)*(V[i+1]-EK)
    I_Na[i+1]=gNamax*(m[i+1]**3)*h[i+1]*(V[i+1]-ENa)
    I_L[i+1]=gL*(V[i+1]-ECl)

# plt.figure()
# plt.plot(t,V[0:len(V)-1], label="membrane potential")
# plt.xlabel("ms")
# plt.ylabel("mV")
# plt.legend()
# plt.axhline(y=V_0,color='r',linestyle='--')

fig = make_subplots(rows=2,cols=2)
fig.add_trace(go.Line(x=t,y=V[0:len(V)-1]),row=1,col=1)
fig.add_trace(go.Line(x=t,y=n[0:len(V)-1]),row=1,col=2)
fig.add_trace(go.Line(x=t,y=m[0:len(V)-1]),row=1,col=2)
fig.add_trace(go.Line(x=t,y=h[0:len(V)-1]),row=1,col=2)
fig.add_trace(go.Line(x=t,y=I_Na[0:len(V)-1]),row=2,col=2)
fig.add_trace(go.Line(x=t,y=I_K[0:len(V)-1]),row=2,col=2)
fig.add_trace(go.Line(x=t,y=Iapplied[0:len(V)-1]),row=2,col=1)
fig.show()

# plt.figure()
# plt.plot(t,m[0:len(V)-1])
# plt.figure()
# plt.plot(t,I_K[0:len(V)-1])

plt.show()
