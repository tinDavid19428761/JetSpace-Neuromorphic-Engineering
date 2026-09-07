# hodgkin huxley validation.py

'''validate model against figure 2.20 in 'dynamical systems in neuroscience book '''

from neuron_functions import gate_inf,timeconstant
import matplotlib.pyplot as plt
import numpy as np

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


# gate plots
V = np.linspace(-40-65,100-65,15)

plt.figure()

plt.plot(V,gate_inf(V,gate_m), label="m_inf", color="blue")
plt.plot(V,gate_inf(V,gate_h), label="h_inf", color="red")
plt.plot(V,gate_inf(V,gate_n), label="n_inf", color="green")
plt.legend()

V2 = np.linspace(-100,100,14*3)

plt.figure()
plt.plot(V2,timeconstant(V2,gate_m), label="tau_m", color="blue")
plt.plot(V2,timeconstant(V2,gate_h), label="tau_h", color="red")
plt.plot(V2,timeconstant(V2,gate_n), label="tau_n", color="green")
plt.legend()

plt.show()
