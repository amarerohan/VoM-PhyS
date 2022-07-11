import numpy as np
from multiprocessing import Pool
import Heat_solver_parameter_file_loader as para
import scipy.sparse as sp
import matplotlib.pyplot as plt
import time as t
import csv

sensitivity_para = para.sensitivity_para

# # # Parameters # # #
t1 = t.time()

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
rho_avg = para.rho_avg
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


row1 = [] ; col1 = [] ; data1 = []
row2 = [] ; col2 = [] ; data2 = []

B_index = [] ; B_data = []



def conduction(ds1 , ds2, ds3):
    UA = 1/(ds1/(2*K) + ds1/(2*K))*ds2*ds3
    UAT = 0
    return(UA , UAT)

def convection(ds1 , ds2 , ds3 , T , h):
    UA = 1/(1/h + ds1/(2*K))*ds2*ds3
    UAT = UA*T
    return(UA , UAT)

def eta(x,C):
    global e
    if(x/e < 1):
        eta_x = C*np.exp(1/((abs(x/e))**2 - 1))#/(e**3)
    else:
        eta_x = 0
    return eta_x


def venal_mass_flow_to_voxel(j,i,k):
     m_dot = 0
     for b in range(len(v_outlets)):
         if([j,i,k] in nbrhd_v[b]):
             ele = v_element.loc[v_element['2'] == v_outlets.iloc[b,0]]
             s = np.sqrt(((i-v_outlets.iloc[b,1])*dx)**2 + ((j-v_outlets.iloc[b,2])*dy)**2 +((k-v_outlets.iloc[b,3])*dz)**2)
             n_e_x = eta(s,Cv[b,0])
             Q = abs(para.QV[ele.iloc[0,0]])  
             m_dot = m_dot + Q*rho*n_e_x 
            
             # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + ele.iloc[0,0]) ; data.append(-m_dot*Cp/2)
        
     return(-m_dot)


def arterial_heat_func():
    global row1
    global col1
    global data1
    global B_index
    global B_data
    
    row = [] ; col = [] ; data = []

    # Arterial Mass Balance and heat transfer
    for i in range(para.n_ae):
        # Convection with surrounding tissue voxel
        UA_sum = 0
        for j in range(len(para.ae_nbr[i])):
            ds = para.ae_nbr[i][j][3]
            dA = para.ae_nbr[i][j][4]
            y,x,z = para.ae_nbr[i][j][0], para.ae_nbr[i][j][1], para.ae_nbr[i][j][2]
            UA = 1/(1/h_tv + ds/(2*K))*dA 
            UA_sum = UA_sum + UA
            
            # row.append(para.n_vxl + i) ; col.append(c_dom[y,x,z]) ; data.append(-UA)
            wr_rcd.writerow([para.n_vxl + i, c_dom[y,x,z], -UA])
        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + i, UA_sum])
        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + i) ; data.append(UA_sum) 
        
        # Cross flow heat excahnge with surrounding vein if any
        UAx_sum = 0
        for j in range(len(para.a_cross[i])):
            ds = para.a_cross[i][j][1]
            dA = para.a_cross[i][j][2]
            UAx = 1/(1/hx + 1/hx)*dA
            UAx_sum = UAx_sum + UAx
            
            # row.append(para.n_vxl + i) ; col.append(para.n_vxl + para.n_ae + para.a_cross[i][j][0]) ; data.append(-UAx)
            wr_rcd.writerow([para.n_vxl + i, para.n_vxl + para.n_ae + para.a_cross[i][j][0], -UAx])
        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + i, UAx_sum])
        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + i) ; data.append(UAx_sum) 
       
        # Mass Balance across elements and nodes
        ele = a_element.loc[a_element['0'] == i]
        sum_m_dotCp = 0
        ele_volume = ((np.pi/4)*(ele.iloc[0,3]*dx)**2)*ele.iloc[0,4]*dx
        # phi = rho*Cp*ele_volume/dt
        # phi_B[para.n_vxl + i] = phi
        if(ele.iloc[0,1] in a_node_in or ele.iloc[0,2] in a_node_in):
            k1 = np.pi*((ele.iloc[0,3]*dx)**4)/(128*myu*(ele.iloc[0,4]*dx))
            m_ele = abs(k1*(Pa[ele.iloc[0,1]] - Pa[ele.iloc[0,2]])*rho)
            
            # B_index.append(para.n_vxl + i)
            # B_data.append(m_ele*Cp*T_in)
            wr_b.writerow([para.n_vxl + i, m_ele*Cp*T_in])
            # B[para.n_vxl + i] = m_ele*Cp*T_in
            # row.append(para.n_vxl + i) ; col.append(para.n_vxl + i) ; data.append(m_ele*Cp) 
            wr_rcd.writerow([para.n_vxl + i, para.n_vxl + i, m_ele*Cp])
            
        
        if(ele.iloc[0,1] in a_node_in or ele.iloc[0,2] in a_node_in):
            k1 = np.pi*((ele.iloc[0,3]*dx)**4)/(128*myu*(ele.iloc[0,4]*dx))
            Q = k1*(Pa[ele.iloc[0,1]] - Pa[ele.iloc[0,2]])
            m_dot = Q*rho
            sum_m_dotCp = sum_m_dotCp + m_dot*Cp
            
        
        if((ele.iloc[0,1] not in a_node_in) and (ele.iloc[0,2] not in a_node_in)):
            k1 = np.pi*((ele.iloc[0,3]*dx)**4)/(128*myu*(ele.iloc[0,4]*dx))
            m_ele = abs(k1*(Pa[ele.iloc[0,1]] - Pa[ele.iloc[0,2]])*rho)
            sum_m_dotCp = sum_m_dotCp + m_ele*Cp 
            
            if(Pa[ele.iloc[0,1]] > Pa[ele.iloc[0,2]]): 
                
                # Node 1 is inlet for the element
                ele_1 = a_element.loc[((a_element['1'] == ele.iloc[0,1]) | (a_element['2'] == ele.iloc[0,1])) & (a_element['0'] != ele.iloc[0,0])]
                sum_m_dot_in = 0
               
                for j in range(len(ele_1)):
                    if((ele.iloc[0,1] == ele_1.iloc[j,1]) and (Pa[ele_1.iloc[j,2]] > Pa[ele.iloc[0,1]])):
                       
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_1.iloc[j,2]] - Pa[ele_1.iloc[j,1]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot 
                            
                    elif((ele.iloc[0,1] == ele_1.iloc[j,2]) and (Pa[ele_1.iloc[j,1]] > Pa[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_1.iloc[j,1]] - Pa[ele_1.iloc[j,2]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot 
                
                for j in range(len(ele_1)):
                    
                    if((ele.iloc[0,1] == ele_1.iloc[j,1]) and (Pa[ele_1.iloc[j,2]] > Pa[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_1.iloc[j,2]] - Pa[ele_1.iloc[j,1]])*rho
                        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + ele_1.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + ele_1.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                    
                    elif((ele.iloc[0,1] == ele_1.iloc[j,2]) and (Pa[ele_1.iloc[j,1]] > Pa[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_1.iloc[j,1]] - Pa[ele_1.iloc[j,2]])*rho
                        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + ele_1.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + ele_1.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
               
                # row.append(para.n_vxl + i) ; col.append(para.n_vxl + i) ; data.append(m_ele*Cp) 
                wr_rcd.writerow([para.n_vxl + i, para.n_vxl + i, m_ele*Cp])
               
            if(Pa[ele.iloc[0,2]] > Pa[ele.iloc[0,1]]):
                # print(ele.iloc[0,0])
                ele_2 = a_element.loc[((a_element['1'] == ele.iloc[0,2]) | (a_element['2'] == ele.iloc[0,2])) & (a_element['0'] != ele.iloc[0,0])]
                sum_m_dot_in = 0
               
                for j in range(len(ele_2)):
                    if( (ele.iloc[0,2] == ele_2.iloc[j,1]) and (Pa[ele_2.iloc[j,2]] > Pa[ele_2.iloc[j,1]]) ):
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_2.iloc[j,2]] - Pa[ele_2.iloc[j,1]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot
                    
                    if((ele.iloc[0,2] == ele_2.iloc[j,2]) and (Pa[ele_2.iloc[j,1]] > Pa[ele_2.iloc[j,2]])):
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_2.iloc[j,1]] - Pa[ele_2.iloc[j,2]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot
                
                for j in range(len(ele_2)):
                    if((ele.iloc[0,2] == ele_2.iloc[j,1]) and (Pa[ele_2.iloc[j,2]] > Pa[ele_2.iloc[j,1]])):
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_2.iloc[j,2]] - Pa[ele_2.iloc[j,1]])*rho 
                        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + ele_2.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + ele_2.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                    if((ele.iloc[0,2] == ele_2.iloc[j,2]) and (Pa[ele_2.iloc[j,1]] > Pa[ele_2.iloc[j,2]])):
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pa[ele_2.iloc[j,1]] - Pa[ele_2.iloc[j,2]])*rho 
                        # row.append(para.n_vxl + i) ; col.append(para.n_vxl + ele_2.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + i, para.n_vxl + ele_2.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                # row.append(para.n_vxl + i) ; col.append(para.n_vxl + i) ; data.append(m_ele*Cp) 
                wr_rcd.writerow([para.n_vxl + i, para.n_vxl + i, m_ele*Cp])
                
    row1 = row
    col1 = col
    data1 = data

def venal_heat_func():
    global row2 
    global col2
    global data2
    
    row = [] ; col = [] ; data = [] 
    # Venal Mass Balance and Heat Tranfer
    for i in range(para.n_ve):
        
        # Convection with surrounding tissue voxel
        UA_sum = 0
        for j in range(len(para.ve_nbr[i])):
            ds = para.ve_nbr[i][j][3]
            dA = para.ve_nbr[i][j][4]
            y,x,z = para.ve_nbr[i][j][0], para.ve_nbr[i][j][1], para.ve_nbr[i][j][2]
            UA = 1/(1/h_tv + ds/(2*K))*dA 
            UA_sum = UA_sum + UA
            
            # row.append(para.n_vxl + para.n_ae + i) ; col.append(c_dom[y,x,z]) ; data.append(-UA)
            wr_rcd.writerow([para.n_vxl + para.n_ae + i, c_dom[y,x,z], -UA])
        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + i, UA_sum])
        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + i) ; data.append(UA_sum) 
        
        # Cross flow heat excahnge with surrounding vein if any
        UAx_sum = 0
        for j in range(len(para.v_cross[i])):
            ds = para.v_cross[i][j][1]
            dA = para.v_cross[i][j][2]
            UAx = 1/(1/hx + 1/hx)*dA
            UAx_sum = UAx_sum + UAx
            
            # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.v_cross[i][j][0]) ; data.append(-UAx)
            wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.v_cross[i][j][0], -UAx])
        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + i, UAx_sum])
        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + i) ; data.append(UAx_sum)
        
        # Mass Balance across elements and nodes
        ele = v_element.loc[v_element['0'] == i]
        sum_m_dotCp = 0
        ele_volume = ((np.pi/4)*(ele.iloc[0,3]*dx)**2)*ele.iloc[0,4]*dx
        # phi = rho*Cp*ele_volume/dt
        # phi_B[para.n_vxl + para.n_ae + i] = phi
        if((ele.iloc[0,1] not in v_out_np[:,0]) and (ele.iloc[0,2] not in v_out_np[:,0])):
            k1 = np.pi*((ele.iloc[0,3]*dx)**4)/(128*myu*(ele.iloc[0,4]*dx)) 
            m_ele = abs(k1*(Pv[ele.iloc[0,1]] - Pv[ele.iloc[0,2]])*rho) 
            sum_m_dotCp = sum_m_dotCp + m_ele*Cp 
            
            if(Pv[ele.iloc[0,1]] > Pv[ele.iloc[0,2]]): 
                
                # Node 1 is inlet for the element
                ele_1 = v_element.loc[((v_element['1'] == ele.iloc[0,1]) | (v_element['2'] == ele.iloc[0,1])) & (v_element['0'] != ele.iloc[0,0])]
                sum_m_dot_in = 0
                
                for j in range(len(ele_1)):
                    if((ele.iloc[0,1] == ele_1.iloc[j,1]) and (Pv[ele_1.iloc[j,2]] > Pv[ele.iloc[0,1]])):
                       
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_1.iloc[j,2]] - Pv[ele_1.iloc[j,1]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot 
                            
                    elif((ele.iloc[0,1] == ele_1.iloc[j,2]) and (Pv[ele_1.iloc[j,1]] > Pv[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_1.iloc[j,1]] - Pv[ele_1.iloc[j,2]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot 
                
                for j in range(len(ele_1)):
                    
                    if((ele.iloc[0,1] == ele_1.iloc[j,1]) and (Pv[ele_1.iloc[j,2]] > Pv[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_1.iloc[j,2]] - Pv[ele_1.iloc[j,1]])*rho
                        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + ele_1.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + ele_1.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                    
                    elif((ele.iloc[0,1] == ele_1.iloc[j,2]) and (Pv[ele_1.iloc[j,1]] > Pv[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_1.iloc[j,3]*dx)**4)/(128*myu*(ele_1.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_1.iloc[j,1]] - Pv[ele_1.iloc[j,2]])*rho
                        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + ele_1.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp)
                        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + ele_1.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                        
               
                # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + i) ; data.append(m_ele*Cp) 
                wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + i, m_ele*Cp])
                
            if(Pv[ele.iloc[0,1]] < Pv[ele.iloc[0,2]]): 
                
                # Node 1 is inlet for the element
                ele_2 = v_element.loc[((v_element['1'] == ele.iloc[0,2]) | (v_element['2'] == ele.iloc[0,2])) & (v_element['0'] != ele.iloc[0,0])]
                sum_m_dot_in = 0
                
                for j in range(len(ele_2)):
                    if((ele.iloc[0,2] == ele_2.iloc[j,1]) and (Pv[ele_2.iloc[j,2]] > Pv[ele.iloc[0,1]])):
                       
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_2.iloc[j,2]] - Pv[ele_2.iloc[j,1]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot 
                            
                    elif((ele.iloc[0,2] == ele_2.iloc[j,2]) and (Pv[ele_2.iloc[j,1]] > Pv[ele.iloc[0,1]])):
                        
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_2.iloc[j,1]] - Pv[ele_2.iloc[j,2]])*rho
                        sum_m_dot_in = sum_m_dot_in + m_dot
                
                for j in range(len(ele_2)):
                    
                    if((ele.iloc[0,2] == ele_2.iloc[j,1]) and (Pv[ele_2.iloc[j,2]] > Pv[ele.iloc[0,2]])):
                        
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_2.iloc[j,2]] - Pv[ele_2.iloc[j,1]])*rho
                        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + ele_2.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + ele_2.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                    
                    elif((ele.iloc[0,2] == ele_2.iloc[j,2]) and (Pv[ele_2.iloc[j,1]] > Pv[ele.iloc[0,2]])):
                        
                        k1 = np.pi*((ele_2.iloc[j,3]*dx)**4)/(128*myu*(ele_2.iloc[j,4]*dx))
                        m_dot = k1*(Pv[ele_2.iloc[j,1]] - Pv[ele_2.iloc[j,2]])*rho
                        # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + ele_2.iloc[j,0]) ; data.append(-m_ele*m_dot/sum_m_dot_in*Cp) 
                        wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + ele_2.iloc[j,0], -m_ele*m_dot/sum_m_dot_in*Cp])
                        
                        
                # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + i) ; data.append(m_ele*Cp) 
                wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + i, m_ele*Cp])
                
        if(ele.iloc[0,2] in v_out_np[:,0]):
            sum_mdot_in = abs(para.QV[ele.iloc[0,0]])*rho
            out_coord = v_outlets.loc[v_outlets['0'] == ele.iloc[0,2]]
            out_n = out_coord.index.values[0]
            
            k1 = np.pi*((ele.iloc[0,3]*dx)**4)/(128*myu*(ele.iloc[0,4]*dx)) 
            m_ele = abs(k1*(Pv[ele.iloc[0,1]] - Pv[ele.iloc[0,2]])*rho)
            
            # for b in range(len(nbrhd_v[out_n])):
            #    s = np.sqrt( ((out_coord.iloc[0,2] - nbrhd_v[out_n][b][0])*dx)**2 + ((out_coord.iloc[0,1] - nbrhd_v[out_n][b][1])*dy)**2 + ((out_coord.iloc[0,3] - nbrhd_v[out_n][b][2])*dz)**2) 
            #   Q = abs(para.QV[ele.iloc[0,0]])#gamma_v/myu*(Pv_com[nbrhd_v[out_n][b][0],nbrhd_v[out_n][b][1],nbrhd_v[out_n][b][2]] - Pv[ele.iloc[0,2]]) 
            #    n_e_x = eta(s,Cv[out_n,0]) 
            #    m_dot = n_e_x*Q*rho 
            #    sum_mdot_in = sum_mdot_in + m_dot 
              
            for b in range(len(nbrhd_v[out_n])):
                s = np.sqrt( ((out_coord.iloc[0,2] - nbrhd_v[out_n][b][0])*dx)**2 + ((out_coord.iloc[0,1] - nbrhd_v[out_n][b][1])*dy)**2 + ((out_coord.iloc[0,3] - nbrhd_v[out_n][b][2])*dz)**2) 
                Q = abs(para.QV[ele.iloc[0,0]])#gamma_v/myu*(Pv_com[nbrhd_v[out_n][b][0],nbrhd_v[out_n][b][1],nbrhd_v[out_n][b][2]] - Pv[ele.iloc[0,2]]) 
                n_e_x = eta(s,Cv[out_n,0]) 
                m_dot = n_e_x*Q*rho 
                
                # row.append(para.n_vxl + para.n_ae + i) ; col.append(c_dom[nbrhd_v[out_n][b][0],nbrhd_v[out_n][b][1],nbrhd_v[out_n][b][2]]) ; data.append(-m_ele*m_dot/sum_mdot_in*Cp) ; 
                wr_rcd.writerow([para.n_vxl + para.n_ae + i, c_dom[nbrhd_v[out_n][b][0],nbrhd_v[out_n][b][1],nbrhd_v[out_n][b][2]], -m_ele*m_dot/sum_mdot_in*Cp])
            wr_rcd.writerow([para.n_vxl + para.n_ae + i, para.n_vxl + para.n_ae + i, m_ele*Cp])
                
            # row.append(para.n_vxl + para.n_ae + i) ; col.append(para.n_vxl + para.n_ae + i) ; data.append(m_ele*Cp) 
            
    row2 = row
    col2 = col
    data2 = data
    
if __name__ == "__main__":
    
    
    xl_title = 'new_method_heat_rcd/other_heat_rcd_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'test_.csv'
    wb_rcd = open(xl_title, 'w')
    wr_rcd = csv.writer(wb_rcd)
    wr_rcd.writerow(['row','col','data'])
    
    xl_title_b = 'new_method_heat_rcd/other_heat_B_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'test_.csv'
    wb_b = open(xl_title_b, 'w')
    wr_b = csv.writer(wb_b)
    wr_b.writerow(['Bi', 'B'])
    
    arterial_heat_func()
    venal_heat_func()
    
    wb_rcd.close()
    wb_b.close()
    
    t2 = t.time()
    
    print(' time = ',round((t2-t1)/60,3),' min')
    
    '''
    arterial_heat_func()
    venal_heat_func()
    
    row = row1 + row2
    col = col1 + col2
    data = data1 + data2
    
    A = sp.csc_matrix((data,(row,col)), shape=(para.N_unk,para.N_unk))
    
    plt.figure(figsize=(10,10))
    plt.spy(A)
    plt.savefig('other_heat.png')
    
    ra_title = 'new_method_heat_rcd/row_other_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    ca_title = 'new_method_heat_rcd/col_other_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    da_title = 'new_method_heat_rcd/data_other_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    
    Bi_title = 'new_method_heat_rcd/Bi_other_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    B_title = 'new_method_heat_rcd/B_other_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    # phi_B_title = 'new_method_heat_rcd/phi_B_other_heat_e_'+str(e)+'_dx_'+str(dx)+'_.npy'
    
    np.save(ra_title,row)
    np.save(ca_title,col)
    np.save(da_title,data)
    np.save(Bi_title,B_index)
    np.save(B_title,B_data)
    # np.save(phi_B_title,phi_B)

    t2 = t.time()
    print('Other \ntime = ',round((t2-t1)/60,3),' min')
    '''
    
    
    
    
 