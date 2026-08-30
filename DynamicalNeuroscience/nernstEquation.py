# nernst equation.py

import math

def nernst_full(T,z,IonIn,IonOut): # K,,Mol,Mol
    R=8315 #mJ/(K*Mol)
    F=96480 #C/Mol
    E_ion=R*T/z/F*math.log(IonOut/IonIn) #mV
    return E_ion

def nernst_abr(z,IonOut,IonIn): #@T=310K or 36.85degC
    E_ion = 62/z*math.log10(IonOut/IonIn)
    return E_ion

#test of E_K+
TdegC = 20 #degC
TKel=TdegC+273.15 #Kelvin
charge = 1
Xin = 430
Xout = 20

# print(f"{nernst_abr(1,20,430):.4f} (mV)")
def nernst_example_print():
    print(f"K+ {nernst_full(TKel,charge,Xin,Xout):.4f} (mV)")
    print(f"Na+ {nernst_full(TKel,charge,50,440):.4f} (mV)")
    print(f"Cl- {nernst_full(TKel,-charge,65,560):.4f} (mV)")

nernst_example_print()