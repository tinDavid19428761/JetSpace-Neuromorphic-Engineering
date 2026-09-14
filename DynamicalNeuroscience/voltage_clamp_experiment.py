# voltage_clamp_experiment.py
'''
variate applied current to hold voltage a constant 'clamp voltage'

modified hodgkin_huxley_model.py for constant dV=0, 
I_applied is now a calculated variable rather than a set input parameter
'''


from neuron_functions import *
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

T=20+273.15 #Kelvin
Capa = 1 #mF/cm^2

# parameters
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


dt = 1/50 #ms
duration = 40 #ms

# initializing data arrays
t = np.arange(0,duration+1,dt)
range = len(t)+1 #+1 accounts for i+1 operation within for-loop
V = np.zeros(range)
m = np.zeros(range)
h = np.zeros(range)
n = np.zeros(range)
I_K = np.zeros(range)
I_Na = np.zeros(range)
I_L = np.zeros(range)


Varray = np.linspace(-40,-90,40)
Itest = np.zeros(len(Varray))
for i,Vstep in enumerate(Varray):
    Itest[i] = -(gKmax*(gate_inf(Vstep,gate_n)**4)*(Vstep-EK))-(gNamax*(gate_inf(Vstep,gate_m)**3)*(gate_inf(Vstep,gate_h))*(Vstep-ENa))-(gL*(Vstep-ECl))

# pre-step potential: (held constant by dynamic I_applied)
V_0 = np.interp(0,Itest,Varray)
V_s=20
print(F"resting membrane potential: {V_0:.4f}")


# Programming voltage step
V.fill(V_0) #pre step voltage
V[round(range*0.2):]=V_s

# intial values of iterated variables
m[0] = m_0 = gate_inf(V_0,gate_m)
h[0] = h_0 = gate_inf(V_0,gate_h)
n[0] = n_0 = gate_inf(V_0,gate_n)

I_K[0]=gKmax*(n_0**4)*(V[0]-EK)
I_Na[0]=gNamax*(m_0**3)*h_0*(V[0]-ENa)
I_L[0]=gL*(V[0]-ECl)

Iapplied = np.zeros(range)


# Simulating neuron membrane potential - iterating differential equation: 
for i,time in enumerate(t):

    dm = gate_diff(V[i],gate_m,m[i])*dt
    dh = gate_diff(V[i],gate_h,h[i])*dt
    dn = gate_diff(V[i],gate_n,n[i])*dt

    m[i+1]=m[i]+dm
    h[i+1]=h[i]+dh
    n[i+1]=n[i]+dn

    I_K[i+1]=gKmax*(n[i+1]**4)*(V[i+1]-EK)
    I_Na[i+1]=gNamax*(m[i+1]**3)*h[i+1]*(V[i+1]-ENa)
    I_L[i+1]=gL*(V[i+1]-ECl)

    # dV, instead of modifying V, will modify I_applied to keep dV=0
    # by counteracting any ion currents
    Iapplied[i+1] = sum([I_K[i+1],I_Na[i+1],I_L[i+1]])


#post processing
sumIonCurrents = I_K+I_Na+I_L
sumCurrent = sumIonCurrents+Iapplied # this will always be zero

fig = make_subplots(rows=2,cols=2)

fig.add_trace(go.Scatter(x=t,y=V[0:len(V)-1], name="clamped membrane potential"),row=1,col=1)
fig.add_trace(go.Scatter(x=t,y=Iapplied[0:len(V)-1], name="Applied Current"),row=2,col=1)
fig.add_trace(go.Scatter(x=t,y=sumIonCurrents[0:len(V)-1], name="Ion currents"),row=2,col=1)

fig.add_trace(go.Scatter(x=t,y=h[0:len(V)-1], name="h - gate parameter"),row=1,col=2)
fig.add_trace(go.Scatter(x=t,y=n[0:len(V)-1], name="n - gate parameter"),row=1,col=2)
fig.add_trace(go.Scatter(x=t,y=m[0:len(V)-1], name="m - gate parameter"),row=1,col=2)


fig.update_traces(textposition="top center")
fig.show()

if V_s > V_0:
    # instantaneous current-voltage I_0
    # steady-state current-voltage I_infinity
    I_0 = min(Iapplied)
    I_inf = max(Iapplied)
    print("V_s > V_0")
elif V_s < V_0:
    I_0 = max(Iapplied)
    I_inf = min(Iapplied)
    print("V_s < V_0")

print(f"I0= {I_0:.2f}")
print(f"Iinfinity= {I_inf:.2f}")