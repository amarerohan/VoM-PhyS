import pandas as pd
import numpy as np


sensitivity_para = 'new'


# Parameters
dx = dy =6.4e-5 #1e-5   # metres
dz = 0.333e-3     # metres


e = 0.0032 # metres   # finite radius of influence
print('e = ',e)
E = e

a = 1e-6       # Perfusion proportionaltiy factor

Ka = 1e-12     # Permeability of arterial compartment tissue
Kv = 5e-10 #1.01*

myu = 3e-3        # Viscosity of blood

gamma_a = 1e-14  # This gamma is used for the terminal boundary condition equation
gamma_v = 1e-14

P_ba = 10600 #Pa
P_bv = 1600 #Pa

Lambda_a = Ka/myu
Lambda_v = Kv/myu


a_inlet1 = 0
a_inlet2 = 24

v_inlet1 = 0
v_inlet2 = 20

a_element = pd.read_csv('arteries_element_database.csv')
v_element = pd.read_csv('veins_element_database.csv')
    
a_outlets = pd.read_csv('arteries_outlet_coordinates_3D_shifted.csv')
v_outlets = pd.read_csv('veins_outlet_coordinates_3D_shifted.csv')

#v_element.iloc[94,3] = 0.001
#v_element.iloc[99,3] = 0.001

    
ca_t = 'constants/Ca_3D_dx_dy_6.4e-05_e_'  + str(e) + '.npy'
cv_t = 'constants/Cv_3D_dx_dy_6.4e-05_e_'  + str(e) + '.npy'

title_nbr_a = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_a_dx_dy_6.4e-05_e_' + str(e) + '_new.npy'
title_nbr_v = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_v_dx_dy_6.4e-05_e_' + str(e) + '_new.npy'
    
Ca = np.load(ca_t)
Cv = np.load(cv_t)

nbrhd_a = np.load(title_nbr_a, allow_pickle=True).tolist()
nbrhd_v = np.load(title_nbr_v, allow_pickle=True).tolist()

dom = np.load('tongue_3D.npy')

ny,nx,nz = np.shape(dom)
c_dom = np.load('c_dom.npy')

Ntvx = np.max(c_dom[:,:,:]) + 1

# Number of elements in the arterial blood vessel network
nPa = len(a_element)
nNa = int(max((max(a_element.iloc[:,1]),max(a_element.iloc[:,2])))) + 1

# Number of elements in the venal blood vessel network
nPv = len(v_element)
nNv = int(max((max(v_element.iloc[:,1]),max(v_element.iloc[:,2])))) + 1

N_unk = 2*Ntvx + nNa + nNv