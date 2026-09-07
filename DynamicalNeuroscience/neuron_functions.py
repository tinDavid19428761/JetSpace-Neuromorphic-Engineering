# neuron_functions.py
import numpy as np
from scipy import signal
def nernst_full(T,z,IonIn,IonOut): # K,,Mol,Mol
    R=8315 #mJ/(K*Mol)
    F=96480 #C/Mol
    E_ion=R*T/z/F*np.log(IonOut/IonIn) #mV
    return E_ion

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
# pulse_mag = 15 #mA
# pulse_period = 10
# pulse_length = 1
# pulse_delay = 5 

# take a time-axis array and generate a PWM-style pulse array.
# pulses from zero to given magnitude
def current_pulser(magnitude,period,length,delay,timeseries):
    frequency = 1/period
    pwm_array = (signal.square(frequency*2*np.pi*timeseries,duty=length/period)+1)/2*magnitude
    return np.where(timeseries >= delay, pwm_array, 0)

