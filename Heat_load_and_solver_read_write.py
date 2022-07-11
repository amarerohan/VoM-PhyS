import numpy as np
import Heat_solver_parameter_file_loader as para
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import time as t
import scipy.sparse.linalg as spla
import os
import pandas as pd

font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_size(15)
# font.set_style('italic')
f = 20


# # # Parameters # # #
t1 = t.time()

sensitivity_para = para.sensitivity_para


dx = para.dx
dy = para.dy
dz = para.dz
a_ele_in = para.a_ele_in
v_ele_out = para.v_ele_out
a_node_in = [0,24]
v_node_out = [0,20]
nx , ny , nz = para.nx , para.ny , para.nz

# Thermal Parameters
q = para.q
K = para.K

h_amb = para.h_amb
h_tv = para.h_tv

hx = para.hx

Cp = para.Cp

T_in = para.T_in
T_amb = para.T_amb
time = para.time

dt = para.dt

# Fluid Parameters
E = para.e_ref
e = E
Ka = para.Ka
Kv = para.Kv
myu = para.myu
rho = para.rho
rho_tissue = para.rho_tissue
gamma_a = para.gamma_a
gamma_v = para.gamma_v
Lambda_a = Ka/myu
Lambda_v = Kv/myu

# Load some important functions  and variables
dom = para.dom
c_dom = para.c_dom
Pa_com = para.P_dom_a
Pv_com = para.P_dom_v
a_element = para.a_element
v_element = para.v_element
Ca = para.Ca
Cv = para.Cv
nbrhd_a = para.nbrhd_a
nbrhd_v = para.nbrhd_v
a_outlets = para.a_outlets
v_outlets = para.v_outlets
Pa = para.Pa
Pv = para.Pv
QA = para.QA
QV = para.QV

a_out_np = a_outlets.to_numpy()
v_out_np = v_outlets.to_numpy()


rcd_voxel = pd.read_csv('new_method_heat_rcd/voxel_heat_rcd_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_test_.csv')
b_voxel = pd.read_csv('new_method_heat_rcd/voxel_heat_B_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_test_.csv')

rcd_other = pd.read_csv('new_method_heat_rcd/other_heat_rcd_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_test_.csv')
b_other = pd.read_csv('new_method_heat_rcd/other_heat_B_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_test_.csv')

r = rcd_voxel.iloc[:,0].tolist() + rcd_other.iloc[:,0].tolist()
c = rcd_voxel.iloc[:,1].tolist() + rcd_other.iloc[:,1].tolist()
d = rcd_voxel.iloc[:,2].tolist() + rcd_other.iloc[:,2].tolist()

B_index = b_voxel.iloc[:,0].tolist() + b_other.iloc[:,0].tolist()
B_data  = b_voxel.iloc[:,1].tolist() + b_other.iloc[:,1].tolist()


'''
row1 = np.load('new_method_heat_rcd/' + 'row_voxel_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
col1 = np.load('new_method_heat_rcd/' + 'col_voxel_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
data1 = np.load('new_method_heat_rcd/' + 'data_voxel_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()

row2 = np.load('new_method_heat_rcd/' + 'row_other_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
col2 = np.load('new_method_heat_rcd/' + 'col_other_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
data2 = np.load('new_method_heat_rcd/' + 'data_other_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()

B_index_1 = np.load('new_method_heat_rcd/' + 'Bi_voxel_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
B_index_2 = np.load('new_method_heat_rcd/' + 'Bi_other_heat_e_' + str(para.e) + '_dx_6.4e-05_.npy').tolist()

B_data_1 = np.load('new_method_heat_rcd/' + 'B_voxel_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()
B_data_2 = np.load('new_method_heat_rcd/' + 'B_other_heat_e_' + str(para.e) + '_dx_6.4e-05_' + sensitivity_para + '_.npy').tolist()


r = row1 + row2
c = col1 + col2
d = data1  + data2
B_index = B_index_1 + B_index_2
B_data = B_data_1 + B_data_2

'''

B_col = np.zeros((len(B_index)),dtype = int).tolist()

A = sp.csc_matrix((d,(r,c)),shape=(para.N_unk,para.N_unk))

B = sp.csr_matrix((B_data,(B_index,B_col)),shape=(para.N_unk,1))
B = B.todense()

B_a = B[0:para.n_vxl]
B_b = B[para.n_vxl:-1]

path = 'new_method_heat_rcd/' + str(para.e_ref)
try:
    os.mkdir(path)
except OSError as error:
    print(error)
    
r_title = path + '/row_k_'+str(K)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+'_Tamb_'+str(T_amb)+'_Tin_'+str(T_in)+ '_'+ sensitivity_para + '_.npy'
c_title = path + '/col_k_'+str(K)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+'_Tamb_'+str(T_amb)+'_Tin_'+str(T_in)+ '_'+ sensitivity_para + '_.npy'
d_title = path + '/data_k_'+str(K)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+'_Tamb_'+str(T_amb)+'_Tin_'+str(T_in)+ '_'+ sensitivity_para + '_.npy'
B_title = path + '/B_k_'+str(K)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+'_Tamb_'+str(T_amb)+'_Tin_'+str(T_in)+ '_'+ sensitivity_para + '_.npy'

np.save(r_title,r)
np.save(c_title,c)
np.save(d_title,d)
np.save(B_title,B)

t2 = t.time()

if(e > 0.001):
    r0 = np.load('new_method_heat_rcd/0.001/row_k_0.49_ht_10.0_h_amb_20.0.npy')
    c0 = np.load('new_method_heat_rcd/0.001/col_k_0.49_ht_10.0_h_amb_20.0.npy')
    d0 = np.load('new_method_heat_rcd/0.001/data_k_0.49_ht_10.0_h_amb_20.0.npy')
    A0 = sp.csc_matrix((d0,(r0,c0)),shape=(para.N_unk,para.N_unk))
    LU = spla.splu(A0)
elif(e<=0.001):
    LU = spla.splu(A)
M = spla.LinearOperator(np.shape(LU),LU.solve)

X0 = np.load('new_method_heat_rcd/0.001/heat_solution_k_0.49_ht_10.0_h_amb_20.0.npy')
X = spla.gmres(A,B,M=M,x0=X0,tol=1e-10)


t3 = t.time()
                 
print(e)
print('error code = ',X[1])
print('time taken to solve the matrix = ',round((t3-t2)/60,3),' mins')
X = X[0]
print('maximum of X = ',round(np.max(X),3),' Minimum of X = ',round(np.min(X),3))
print('Using A0 of 0.001 for LU')

    
np.save(path + '/heat_solution_k_'+str(K)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+'_Tamb_'+str(T_amb)+'_Tin_'+str(T_in)+ '_'+ sensitivity_para + '_.npy',X)

T = np.zeros((ny,nx,nz),dtype = float)
T[:,:,:] = 0
Ta = np.copy(T)
Tv = np.copy(T)
for k in range(para.nz):
    for j in range(para.ny):
        for i in range(para.nx):
            if(dom[j,i,k] == 1):
                T[j,i,k] = X[c_dom[j,i,k]]
            if(dom[j,i,k] == 2):
                T[j,i,k] = X[para.n_vxl + para.av_ele_ref[j,i,0]]
                Ta[j,i,k] = X[para.n_vxl + para.av_ele_ref[j,i,0]]
            if(dom[j,i,k] == 3):
                T[j,i,k] = X[para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1]]
                Tv[j,i,k] = X[para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1]]




# Energy Conservation
surf_T = []
energy = []
tissue_volume = 0
for k in range(para.nz):
    for j in range(para.ny):
        for i in range(para.nx):
            if(dom[j,i,k] == 1):
                
                tissue_volume = tissue_volume + dy*dx*dz
                
                if(dom[j+1,i,k] == 0):
                    UA = (1/(dy/(2*K) + 1/h_amb))*dx*dz
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell) 
                    surf_T.append(T[j,i,k])
                
                if(dom[j-1,i,k] == 0):
                    UA = (1/(dy/(2*K) + 1/h_amb))*dx*dz
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell)
                    surf_T.append(T[j,i,k])
                    
                if(dom[j,i+1,k] == 0):
                    UA = (1/(dx/(2*K) + 1/h_amb))*dy*dz
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell)
                    surf_T.append(T[j,i,k])
                    
                if(dom[j,i-1,k] == 0):
                    UA = (1/(dx/(2*K) + 1/h_amb))*dy*dz
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell)
                    surf_T.append(T[j,i,k])
                
                if(dom[j,i,k+1] == 0):
                    UA = (1/(dz/(2*K) + 1/h_amb))*dx*dy
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell)
                    surf_T.append(T[j,i,k])
                    
                if(dom[j,i,k-1] == 0):
                    UA = (1/(dz/(2*K) + 1/h_amb))*dx*dy
                    q_cell = UA*(T[j,i,k] - T_amb)
                    energy.append(q_cell)
                    surf_T.append(T[j,i,k])
                

q_in = abs(para.QA[a_ele_in[0]])*rho*Cp*(T_in) + abs(para.QA[a_ele_in[1]])*rho*Cp*(T_in)
q_out = abs(para.QV[v_ele_out[0]])*rho*Cp*X[para.n_vxl + para.n_ae + v_ele_out[0]] + abs(para.QV[v_ele_out[1]])*rho*Cp*X[para.n_vxl + para.n_ae + para.v_ele_out[1]]
energy_in = q_in + tissue_volume*q
energy_out = abs(sum(energy)) + abs(q_out)
energy_convection = sum(energy)
energy_mass_out = q_out
energy_mass_in = q_in
energy_gen = tissue_volume*q

mass_in = abs(para.QA[a_ele_in[0]]*rho) + abs(para.QA[a_ele_in[1]]*rho)
mass_out = abs(para.QV[v_ele_out[0]]*rho) + abs(para.QV[v_ele_out[1]]*rho) 

energy_error = abs(energy_in) - abs(energy_out)

print('\nenergy conservation error = ',round(energy_error,5))

print('\n')
print('mass in  = ', round(mass_in,5))
print('mass out = ', round(mass_out,5))
print('\n')
print('energy generated = ',round(energy_gen,5))
print('energy mass in  = ',round(energy_mass_in,5))
print('energy mass out = ',round(energy_mass_out,5))
print('energy convection = ',round(energy_convection,5))
print('\n')
print('energy in   = ',round(energy_in,5))
print('energy out  = ',round(energy_out,5))
print('energy difference = ',round(energy_error,5))
print('\n')

print('Inlet Temperature = ',round(X[para.n_vxl + a_ele_in[0]],3),round(X[para.n_vxl + a_ele_in[1]],3))
print('Outlet Temperature = ',round(X[para.n_vxl + para.n_ae + v_ele_out[0]],3),round(X[para.n_vxl + para.n_ae + v_ele_out[1]],3))


TT = X[0:para.n_vxl]
TA = X[para.n_vxl:para.n_vxl + para.n_ae]
TV = X[para.n_vxl + para.n_ae:para.N_unk]

print('Maximum Arterial Temperature = ',round(np.max(TA[:]),5))
print('Maximum Venal Temperature = ',round(np.max(TV[:]),5))
print('Maximum Tissue Temperature = ',round(np.max(TT[:]),5),'\n')


'''
'''
# # # make directory to save data
path = 'new_method_heat_rcd/' + str(para.e_ref)
try:
    os.mkdir(path)
except OSError as error:
    print(error)
    
TAtitle = path + '/Temperature_Arterial_Element_e_'+str(e)+'_q_'+str(q)+'_Tin_'+str(T_in)+'_T_amb_'+str(T_amb)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+ '_'+ sensitivity_para + '.npy'
TVtitle = path + '/Temperature_Venal_Element_e_'+str(e)+'_q_'+str(q)+'_Tin_'+str(T_in)+'_T_amb_'+str(T_amb)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+ '_'+ sensitivity_para +  '.npy'
TTtitle = path + '/Temperature_Tissue_Element_e_'+str(e)+'_q_'+str(q)+'_Tin_'+str(T_in)+'_T_amb_'+str(T_amb)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+ '_'+ sensitivity_para +  '.npy'
T_title = path +'/Temperature_Domain_e_'+str(e)+'_q_'+str(q)+'_Tin_'+str(T_in)+'_T_amb_'+str(T_amb)+'_ht_'+str(h_tv)+'_h_amb_'+str(h_amb)+ '_'+ sensitivity_para + '.npy'

np.save(TAtitle,TA)
np.save(TVtitle,TV)
np.save(TTtitle,TT)
np.save(T_title,T)

#print('Maximum Temperature at: ', np.where(T[:,:,:] == np.amax(T)))
