import numpy as np
import pandas as pd

sensitivity_para = 'new'

fluid_para = ['myu','a','Ka','Kv','Ga','Gv', 'new']
heat_para = ['q', 'K', 'rho', 'rho_tissue', 'Cp', 'T_in', 'h_amb']
# # # # Parameters # # #
dx = 6.4e-05
dy = 6.4e-05
dz = 0.333e-3
a_ele_in = [0,23]
v_ele_out = [0,19]
e_ref = 0.0032
e = e_ref


# Thermal Parameters
q = 1000.0
K = 0.49

h_amb = 20.0
h_tv = 10.0
hx = 1E-10

Cp = 3421 

T_in = 10.00
T_amb = 10.0

time = 120
dt = 1

# Fluid Parameters
Ka = 1e-12
Kv = 5e-10
myu = 3e-3
rho = 1050 
rho_tissue = 1050 #1090 
rho_avg = (rho + rho_tissue)/2
gamma_a = 1e-14
gamma_v = 1e-14
Lambda_a = Ka/myu
Lambda_v = Kv/myu

# # # # Load Documents # # #
na_title = 'nbrhd_matrices/'+str(e_ref)+'/nbrhd_3D_a_dx_dy_6.4e-05_e_' + str(e) + '_new.npy'
nv_title = 'nbrhd_matrices/'+str(e_ref)+'/nbrhd_3D_v_dx_dy_6.4e-05_e_' + str(e) + '_new.npy'
nbrhd_a = np.load(na_title, allow_pickle=True).tolist()
nbrhd_v = np.load(nv_title, allow_pickle=True).tolist()

ca_title = 'constants/Ca_3D_dx_dy_6.4e-05_e_' + str(e) + '.npy'
cv_title = 'constants/Cv_3D_dx_dy_6.4e-05_e_' + str(e) + '.npy'
Ca = np.load(ca_title)
Cv = np.load(cv_title)

a_element = pd.read_csv('arteries_element_database.csv')
a_outlets = pd.read_csv('arteries_outlet_coordinates_3D_shifted.csv')
v_element = pd.read_csv('veins_element_database.csv')
v_outlets = pd.read_csv('veins_outlet_coordinates_3D_shifted.csv')



# # creating domain

# arteries = np.load('Frog_Tongue_Arteries_reference.npy')
# veins = np.load('Frog_Tongue_Veins_reference.npy')
# tongue = np.load('Frog_Tongue_Tissue.npy')
dom = np.load('tongue_3D.npy')
ny,nx,nz = np.shape(dom)


av_ele_ref = np.zeros((ny,nx,2),dtype = int)

av_ele_ref[:,:,0] = np.load('Frog_Arteries.npy')
av_ele_ref[:,:,1] = np.load('Frog_Veins.npy')
# np.save('av_ele_ref.npy',av_ele_ref)

# e = e*dx

if(sensitivity_para in fluid_para):

    arterial_compartment_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Compartment_Pressure_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'   
    venal_compartment_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Compartment_Pressure_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'        
    arterial_nodal_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Nodal_Pressure_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'   
    venal_nodal_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Nodal_Pressure_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'
    arterial_flow_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Elemental_Flow_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'   
    venal_flow_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Elemental_Flow_'+str(e_ref)+'_sensitivity_'+sensitivity_para+'_.npy'   

elif(sensitivity_para in heat_para):
    arterial_compartment_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Compartment_Pressure_'+str(e_ref)+'_.npy'   
    venal_compartment_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Compartment_Pressure_'+str(e_ref)+'_.npy'        
    arterial_nodal_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Nodal_Pressure_'+str(e_ref)+'_.npy'   
    venal_nodal_pressure_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Nodal_Pressure_'+str(e_ref)+'_.npy'
    arterial_flow_title = 'new_method_flow_rcd/'+str(e_ref)+'/Arterial_Elemental_Flow_'+str(e_ref)+'_.npy'   
    venal_flow_title = 'new_method_flow_rcd/'+str(e_ref)+'/Venal_Elemental_Flow_'+str(e_ref)+'_.npy'   
    
    
QA = np.load(arterial_flow_title)
QV = np.load(venal_flow_title)

    
P_dom_a = np.zeros((ny,nx,nz),dtype = float)
P_dom_v = np.zeros((ny,nx,nz),dtype = float)
P_dom_a = np.load(arterial_compartment_pressure_title)
P_dom_v = np.load(venal_compartment_pressure_title)

Pa = np.load(arterial_nodal_pressure_title)
Pv = np.load(venal_nodal_pressure_title)



# # Number of unknowns and their order
n_ae = np.max(av_ele_ref[:,:,0]) + 1    # Number of arterial elements
n_ve = np.max(av_ele_ref[:,:,1]) + 1    # Numbere of venal elements
n_a_i = len(a_ele_in)
n_v_o = len(v_ele_out)


ny,nx,nz = np.shape(dom)
#c_dom = np.zeros((ny,nx,nz),dtype = int)
c_dom = np.load('c_dom.npy')
'''
voxel_counter = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            if(dom[j,i,k] == 1):
                c_dom[j,i,k] = voxel_counter
                voxel_counter = voxel_counter + 1
'''
# del(voxel_counter)
n_vxl = np.max(c_dom[:,:,:]) + 1

N_unk = n_vxl + n_ae + n_ve # Total number of unknowns

a_cross = np.load('Arterial_Cross_flow_neighbour.npy',allow_pickle=True).tolist()
v_cross = np.load('Venal_Cross_flow_neighbour.npy',allow_pickle=True).tolist()
ae_nbr = np.load('Arterial_Element_Tissue_neighbour.npy',allow_pickle=True).tolist()
ve_nbr = np.load('Venal_Element_Tissue_neighbour.npy',allow_pickle=True).tolist()


