import numpy as np
import scipy as sp
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import Heat_solver_parameter_file_loader as para
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import os
from multiprocessing import Pool, Manager
import time as t
import csv

t1 = t.time()

sensitivity_para = para.sensitivity_para

# # # Parameters # # #
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

N_unk = para.N_unk

# # # Functions # # # 
def mass_flow_across_voxels(nbr , cell , ds , ds2 , ds3):
    global Pa_com
    global Pv_com
    
    t_ij_a = 2*ds2*ds3*(1/(ds/Lambda_a + ds/Lambda_a))  
    Qa = t_ij_a*(Pa_com[nbr[0],nbr[1],nbr[2]] - Pa_com[cell[0],cell[1],cell[2]]) 
    ma = rho_avg*Qa
    
    t_ij_v = 2*ds2*ds3*(1/(ds/Lambda_v + ds/Lambda_v))
    Qv = t_ij_v*(Pv_com[nbr[0],nbr[1],nbr[2]] - Pv_com[cell[0],cell[1],cell[2]])
    mv = rho_avg*Qv
    
    m = ma + mv 
    return (m)


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

def arterial_mass_flow_to_voxel(j,i,k,row,col,data, wr_rcd):
    # global row
    # global col
    # global data
    m_dot = 0
    for b in range(len(a_outlets)):
        if([j,i,k] in nbrhd_a[b]):
            
            ele = a_element.loc[a_element['2']==a_outlets.iloc[b,0]]
            s = np.sqrt(((i-a_outlets.iloc[b,1])*dx)**2 + ((j-a_outlets.iloc[b,2])*dy)**2 +((k-a_outlets.iloc[b,3])*dz)**2)
            n_e_x = eta(s,Ca[b,0]) 
            Q = abs(para.QA[ele.iloc[0,0]]) 
            m_dot = m_dot + Q*rho*n_e_x 
            
            # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + ele.iloc[0,0]) ; data.append(Q*rho*n_e_x*Cp) 
            wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + ele.iloc[0,0], Q*rho*n_e_x*Cp])
    
    return(m_dot)

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


def func(input_values):
    nx,ny_array,nz,wr_rcd, wr_b = input_values
    ny1 = ny_array[0]
    ny2 = ny_array[1]
    
    row = [] ; col = [] ; data = []
    B_index = [] ; B_data = []
    
    
    for k in range(nz):
        for j in range(ny1,ny2):
            for i in range(nx):
                if(dom[j,i,k] == 1):   # Voxel is tissue not blood vessel
                    # vxl_volume = dx*dy*dz
                    # phi = rho_tissue*Cp*vxl_volume/dt
                
                    mav_n = mav_s = mav_e = mav_w = mav_f = mav_b = 0.0
                    UA_n = UA_s = UA_e = UA_w = UA_f = UA_b = 0.0
                    UAT_n = UAT_s = UAT_e = UAT_w = UAT_f = UAT_b = 0.0
                    UA_na = UA_sa = UA_ea = UA_wa = UA_fa = UA_ba = 0.0
                    UA_nv = UA_sv = UA_ev = UA_wv = UA_fv = UA_bv = 0.0
                    UAT_na = UAT_sa = UAT_ea = UAT_wa = UAT_fa = UAT_ba = 0.0
                    UAT_nv = UAT_sv = UAT_ev = UAT_wv = UAT_fv = UAT_bv = 0.0
                    mass_cons = 0
                    
                    
                    #if(j < ny-1):
                    # North is a voxel
                    if(dom[j+1,i,k] == 1):
                        mav_n = mass_flow_across_voxels([j+1,i,k] , [j,i,k] , dy , dx, dz)
                        UA_n , UAT_n= conduction(dy,dx,dz)
                        mass_cons = mass_cons + mav_n
                        if(mav_n > 0):
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j+1,i,k], mav_n*Cp])
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j+1,i,k]) ; data.append(mav_n*Cp)
                        if(mav_n < 0):
                            mav_n = 0
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j+1,i,k], UA_n])
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j+1,i,k]) ; data.append(UA_n)
                    
                    if(dom[j+1,i,k] == 0):
                        mav_n = 0
                        UA_n , UAT_n= convection(dy,dx,dz,T_amb,h_amb) 
                        
                    # North is an artery
                    if(dom[j+1,i,k] == 2):
                        UA_na , UAT_na = convection(dy,dx,dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j+1,i,0]) ; data.append(UA_na)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j+1,i,0], UA_na])
                        
                    # North is a vein
                    if(dom[j+1,i,k] == 3): 
                        UA_nv , UAT_nv = convection(dy,dx,dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j+1,i,1]) ; data.append(UA_nv)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j+1,i,1], UA_nv])
                       
                        
                       
                        
                    # if(j > 0):
                    # South is a voxel
                    if(dom[j-1,i,k] == 1):
                        mav_s = mass_flow_across_voxels([j-1,i,k] , [j,i,k] , dy , dx, dz)
                        UA_s , UAT_s = conduction(dy , dx , dz)
                        mass_cons = mass_cons + mav_s
                        if(mav_s > 0):
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j-1,i,k]) ; data.append(mav_s*Cp)
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j-1,i,k], mav_s*Cp])
                        if(mav_s < 0):
                            mav_s = 0
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j-1,i,k]) ; data.append(UA_s)
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j-1,i,k], UA_s])
                    
                    if(dom[j-1,i,k] == 0):
                        mav_s = 0
                        UA_s , UAT_s= convection(dy,dx,dz,T_amb,h_amb)
                        
                    # South is an artery
                    if(dom[j-1,i,k] == 2):
                        UA_sa , UAT_sa = convection(dy, dx , dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j-1,i,0]) ; data.append(UA_sa)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j-1,i,0], UA_sa])
                    
                    # South is a vein
                    if(dom[j-1,i,k] == 3):
                        UA_sv , UAT_sv = convection(dy,dx,dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j-1,i,1]) ; data.append(UA_sv)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j-1,i,1], UA_sv])
                            
                        
                        
                        
                    # if(i < nx-1):
                    # East is a voxel
                    if(dom[j,i+1,k] == 1):
                        mav_e = mass_flow_across_voxels([j,i+1,k] , [j,i,k] , dx ,dy, dz)
                        UA_e , UAT_e = conduction(dx , dy , dz)
                        mass_cons = mass_cons + mav_e
                        if(mav_e > 0):
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i+1,k]) ; data.append(mav_e*Cp)
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i+1,k], mav_e*Cp])
                        if(mav_e < 0):
                            mav_e = 0
                            
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i+1,k]) ; data.append(UA_e)
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i+1,k], UA_e])
                    
                    if(dom[j,i+1,k] == 0):
                        mav_e = 0
                        UA_e , UAT_e= convection(dx,dy,dz,T_amb,h_amb)
                        
                    # East is an artery
                    if(dom[j,i+1,k] == 2):
                        UA_ea , UAT_ea = convection(dx, dy , dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j,i+1,0]) ; data.append(UA_ea)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j,i+1,0], UA_ea])
                    
                    # East is a vein
                    if(dom[j,i+1,k] == 3):
                        UA_ev , UAT_ev = convection(dx,dy,dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j,i+1,1]) ; data.append(UA_ev) 
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j,i+1,1], UA_ev])
                            
                        
                        
                        
                    # if(i > 0):
                    # West is a voxel
                    if(dom[j,i-1,k] == 1):
                        mav_w = mass_flow_across_voxels([j,i-1,k] , [j,i,k] , dx ,dy, dz)
                        UA_w , UAT_w = conduction(dx , dy , dz)
                        mass_cons = mass_cons + mav_w
                        if(mav_w > 0):
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i-1,k]) ; data.append(mav_w*Cp)
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i-1,k], mav_w*Cp])
                        if(mav_w < 0):
                            mav_w = 0
                            
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i-1,k]) ; data.append(UA_w)
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i-1,k], UA_w])
                    
                    if(dom[j,i-1,k] == 0):
                        mav_w = 0
                        UA_w , UAT_w= convection(dx,dy,dz,T_amb,h_amb)
                        
                    # West is an artery
                    if(dom[j,i-1,k] == 2):
                        UA_wa , UAT_wa = convection(dx, dy , dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j,i-1,0]) ; data.append(UA_wa)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j,i-1,0], UA_wa])
                    
                    # West is a vein
                    if(dom[j,i-1,k] == 3):
                        UA_wv , UAT_wv = convection(dx,dy,dz,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j,i-1,1]) ; data.append(UA_wv) 
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j,i-1,1], UA_wv])
                    
                    
                    
                    
                    # if(k<nz-1):
                    # front is voxel
                    if(dom[j,i,k+1] == 1):
                        mav_f = mass_flow_across_voxels([j,i,k+1] , [j,i,k] , dz ,dy, dx)
                        UA_f , UAT_f = conduction(dz , dy , dx)
                        mass_cons = mass_cons + mav_f
                        if(mav_f > 0):
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i,k+1]) ; data.append(mav_f*Cp)
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i,k+1], mav_f*Cp])
                        if(mav_f < 0):
                            mav_f = 0
                            
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i,k+1]) ; data.append(UA_f)
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i,k+1], UA_f])
                        
                    if(dom[j,i,k+1] == 0):
                        mav_f = 0
                        UA_f , UAT_f = convection(dz , dx, dy, T_amb,h_amb)
                        
                    if(dom[j,i,k+1] == 2):
                        UA_fa , UAT_fa = convection(dz,dx,dy,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j,i,0]) ; data.append(UA_fa)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j,i,0], UA_fa])
                        
                    if(dom[j,i,k+1] == 3):
                        UA_fv , UAT_fv = convection(dz,dy,dx,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1]) ; data.append(UA_fv) 
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1], UA_fv])
                    
                    
                    
                    
                    # if(k>0):
                    # back is voxel
                    if(dom[j,i,k-1] == 1):
                        mav_b = mass_flow_across_voxels([j,i,k-1] , [j,i,k] , dz ,dy, dx)
                        UA_b , UAT_b = conduction(dz , dy , dx)
                        mass_cons = mass_cons + mav_b
                        if(mav_b > 0):
                            # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i,k-1]) ; data.append(mav_b*Cp)
                            wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i,k-1], mav_b*Cp])
                        if(mav_b < 0):
                            mav_b = 0
                            
                        # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i,k-1]) ; data.append(UA_b)
                        wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i,k-1], UA_b])
                    if(dom[j,i,k-1] == 0):
                        mav_b = 0
                        UA_b , UAT_b = convection(dz , dx, dy, T_amb,h_amb)
                        
                    if(dom[j,i,k-1] == 2):
                        UA_ba , UAT_ba = convection(dz,dx,dy,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.av_ele_ref[j,i,0]) ; data.append(UA_ba)
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.av_ele_ref[j,i,0], UA_ba])
                        
                    if(dom[j,i,k-1] == 3):
                        UA_bv , UAT_bv = convection(dz,dy,dx,0,h_tv)
                        
                        # row.append(c_dom[j,i,k]) ; col.append(para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1]) ; data.append(UA_bv) 
                        wr_rcd.writerow([c_dom[j,i,k], para.n_vxl + para.n_ae + para.av_ele_ref[j,i,1], UA_bv])
                
                    
                    
                    # Mass flow directly between voxel and artery or voxel and veiin
                    mA = arterial_mass_flow_to_voxel(j,i,k,row,col,data, wr_rcd)
                    # mV = venal_mass_flow_to_voxel(j,i,k)
                    sum_mav = mav_n + mav_s  + mav_e + mav_w + mav_f + mav_b + mA 
                    sum_UA = UA_n + UA_s + UA_e + UA_w + UA_f + UA_b + UA_na + UA_sa + UA_ea + UA_wa + UA_fa + UA_ba + UA_nv + UA_sv + UA_ev + UA_wv + UA_bv + UA_fv
                    sum_UAT = UAT_n + UAT_s + UAT_e + UAT_w + UAT_f + UAT_b + UAT_na + UAT_sa + UAT_ea + UAT_wa + UAT_ba + UAT_fa + UAT_nv + UAT_sv + UAT_ev + UAT_wv + UAT_bv + UAT_fv
                    
                    #mass_cons = mass_cons + mA + mV
                    #if(((mass_cons) > 1e-19)):
                    #    print([j,i,k], round(mass_cons,20))
                    # row.append(c_dom[j,i,k]) ; col.append(c_dom[j,i,k]) ; data.append(-(sum_UA) - sum_mav*Cp) 
                    wr_rcd.writerow([c_dom[j,i,k], c_dom[j,i,k], -(sum_UA) - sum_mav*Cp])
                    
                    # B_index.append(c_dom[j,i,k]) 
                    # B_data.append(-q*dx*dy*dz - sum_UAT)
                    wr_b.writerow([c_dom[j,i,k], -q*dx*dy*dz - sum_UAT])
                    
                    # B[c_dom[j,i,k],0] = -q*dx*dy*dz - sum_UAT
                    # phi_B[c_dom[j,i,k],0] = -phi
    
    
    return(row,col,data,B_index,B_data)


if __name__ == "__main__":
    xl_title = 'new_method_heat_rcd/voxel_heat_rcd_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'test_.csv'
    wb_rcd = open(xl_title, 'w')
    wr_rcd = csv.writer(wb_rcd)
    wr_rcd.writerow(['row','col','data'])
    
    xl_title_b = 'new_method_heat_rcd/voxel_heat_B_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'test_.csv'
    wb_b = open(xl_title_b, 'w')
    wr_b = csv.writer(wb_b)
    wr_b.writerow(['Bi', 'B'])
    
    test = func([nx,[0,ny],nz,wr_rcd, wr_b])
    
    wb_rcd.close()
    wb_b.close()
    
    t2 = t.time()
    
    print(' time = ',round((t2-t1)/60,3),' min')
    
    '''
    
    rcdb = []
    dny = 32
    pool_n = 32
    ny_list = []
    
    for i in range(dny):
        ny_list.append([int(i*ny/dny),int((i+1)*ny/dny)])
    
    p = Pool(pool_n)
    input_value = [(nx,ny_list[x],nz) for x in range(dny)]
    rcdb.append(p.map(func,input_value))
    p.close()
    p.join()
    
    row = []
    col = []
    data = []
    B_index = []
    B_data = []
    N_unk = para.N_unk
    
    for i in range(len(rcdb[0])):
        row = row + rcdb[0][i][0]
        col = col + rcdb[0][i][1]
        data = data + rcdb[0][i][2]
        B_index = B_index + rcdb[0][i][3]
        B_data = B_data + rcdb[0][i][4]
        
    r_title = 'new_method_heat_rcd/row_voxel_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    c_title = 'new_method_heat_rcd/col_voxel_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    d_title = 'new_method_heat_rcd/data_voxel_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    bi_title = 'new_method_heat_rcd/Bi_voxel_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    b_title = 'new_method_heat_rcd/B_voxel_heat_e_'+str(e)+'_dx_'+str(dx)+ '_'+ sensitivity_para + '_.npy'
    
    np.save(r_title,row)
    np.save(c_title,col)
    np.save(d_title,data)
    np.save(bi_title,B_index)
    np.save(b_title,B_data)
    
    t2 = t.time()
    
    print('Voxel Heat Equations\n pool number = '+ str(pool_n)+' dny = '+str(dny)+' time = ',round((t2-t1)/60,3),' min')
    
    print('max of row = ',max(row))
    print('max of col = ',max(col))
    
    # A = sp.csc_matrix((data,(row,col)), shape=(N_unk,N_unk))
    
    # plt.figure(figsize=(10,10))
    # plt.spy(A)
    # plt.savefig('heat_voxels.png')
    
    '''
     
