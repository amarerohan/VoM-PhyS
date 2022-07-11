import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import Pool, Manager
import time
import Flow_solver_parameter_file_loader as para
import csv

sensitivity_para = para.sensitivity_para
t1 = time.time()

# Parameters
dx = para.dx
dy = para.dy
dz = para.dz


e = para.e
print('e = ',e)
E = e

a = para.a

Ka = para.Ka
Kv = para.Kv

myu = para.myu

gamma_a = para.gamma_a
gamma_v = para.gamma_v

P_ba = para.P_ba
P_bv = para.P_bv

Lambda_a = para.Lambda_a
Lambda_v = para.Lambda_v


a_inlet1 = para.a_inlet1
a_inlet2 = para.a_inlet2

v_inlet1 = para.v_inlet1
v_inlet2 = para.v_inlet2

a_element = para.a_element
v_element = para.v_element
    
a_outlets = para.a_outlets
v_outlets = para.v_outlets

# v_element.iloc[94,3] = 0.001
# v_element.iloc[99,3] = 0.001

    
Ca = para.Ca # np.load(ca_t)
Cv = para.Cv #np.load(cv_t)



nbrhd_a = para.nbrhd_a # np.load(title_nbr_a, allow_pickle=True).tolist()
nbrhd_v = para.nbrhd_v # np.load(title_nbr_v, allow_pickle=True).tolist()


dom = para.dom

ny,nx,nz = np.shape(dom)
c_dom = para.c_dom



Ntvx = para.Ntvx


# Number of elements in the arterial blood vessel network
nPa = len(a_element)
nNa = int(max((max(a_element.iloc[:,1]),max(a_element.iloc[:,2])))) + 1

# Number of elements in the venal blood vessel network
nPv = len(v_element)
nNv = int(max((max(v_element.iloc[:,1]),max(v_element.iloc[:,2])))) + 1


# Matrix generation

def eta(x,nbrhd,C):
    global e
    if(x/e < 1):
        eta_x = C*np.exp(1/((abs(x/e))**2 -1))#/(e**3)#*(dx*dy*dz) 
    else:
        eta_x = 0
    
    return eta_x 



def arterial_nodal_mass_balance():
    # row = [] ; col = [] ; data = []
    eqn_n = 2*Ntvx 
    a_out = a_outlets.to_numpy()
    for i in range(0,nNa):
        sum_ij = 0
        if(((i != a_inlet1 and i!= a_inlet2)) and (i not in a_out[:,0])):
            for b in range(len(a_element)):
                if(i == a_element.iloc[b,1]): 
                    phi = np.pi*(a_element.iloc[b,3]*dx)**4/(128*myu*(a_element.iloc[b,4])*dx)
                    sum_ij = sum_ij + phi
                    ws.writerow([eqn_n, 2*Ntvx + a_element.iloc[b,2], -phi])
                    # data.append(-phi) ; row.append(eqn_n) ; col.append(2*Ntvx + a_element.iloc[b,2]) 
                if(i == a_element.iloc[b,2]):
                    phi = np.pi*(a_element.iloc[b,3]*dx)**4/(128*myu*(a_element.iloc[b,4])*dx)
                    sum_ij = sum_ij + phi
                    ws.writerow([eqn_n, 2*Ntvx + a_element.iloc[b,1], -phi])
                    # data.append(-phi) ; row.append(eqn_n) ; col.append(2*Ntvx + a_element.iloc[b,1]) 
            
            ws.writerow([eqn_n, 2*Ntvx + i, sum_ij])
            # data.append(sum_ij) ; row.append(eqn_n) ; col.append(2*Ntvx + i) 
            eqn_n = eqn_n + 1
    
    return()

# Veins
def venal_nodal_mass_balance():
    # row = [] ; col = [] ; data = []
    
    eqn_n = 2*Ntvx + (nNa - 2 - len(a_outlets))
    v_out = v_outlets.to_numpy()
    for i in range(0,nNv):
        sum_ij = 0
        if((i != v_inlet1 and i!= v_inlet2) and i not in v_out[:,0]):
            for b in range(len(v_element)):
                if(i == v_element.iloc[b,1]): 
                    phi = np.pi*(v_element.iloc[b,3]*dx)**4/(128*myu*(v_element.iloc[b,4]*dx))
                    sum_ij = sum_ij + phi
                    ws.writerow([eqn_n, 2*Ntvx + nNa + v_element.iloc[b,2], -phi])
                    # data.append(-phi) ; row.append(eqn_n) ; col.append(2*Ntvx + nNa + v_element.iloc[b,2])
                if(i == v_element.iloc[b,2]):
                    phi = np.pi*(v_element.iloc[b,3]*dx)**4/(128*myu*(v_element.iloc[b,4]*dx))
                    sum_ij = sum_ij + phi
                    ws.writerow([eqn_n, 2*Ntvx + nNa + v_element.iloc[b,1], -phi])
                    # data.append(-phi) ; row.append(eqn_n) ; col.append(2*Ntvx + nNa + v_element.iloc[b,1])
            ws.writerow([eqn_n, 2*Ntvx + nNa + i, sum_ij])
            # data.append(sum_ij) ; row.append(eqn_n) ; col.append(2*Ntvx + nNa + i)
            eqn_n = eqn_n + 1
    
    return()

# Terminal nodes boundary conditions
# Aretries
def terminal_arteries():
    # row = [] ; col = [] ; data = []
    
    eqn_n = 2*Ntvx + (nNa - 2 - len(a_outlets)) + (nNv - 2 - len(v_outlets))
    for o in range(len(a_outlets)):
        i = a_outlets.iloc[o,0]
        ele = a_element.loc[a_element['2'] == i]
        k = np.pi*(ele.iloc[0,3]*dx)**4/(128*myu*(ele.iloc[0,4]*dx))
        
        ws.writerow([eqn_n, 2*Ntvx + ele.iloc[0,1], k])
        ws.writerow([eqn_n, 2*Ntvx + i, -(k + gamma_a/myu)])
        # data.append(k) ; row.append(eqn_n) ; col.append(2*Ntvx + ele.iloc[0,1])
        # data.append(-(k + gamma_a/myu)) ; row.append(eqn_n) ; col.append(2*Ntvx + i)
        
        for j in range(len(nbrhd_a[o])):
            y,x,z = nbrhd_a[o][j]
            s = np.sqrt(((a_outlets.iloc[o,1] -  x)*dx)**2 + ((a_outlets.iloc[o,2] -  y)*dy)**2 + ((a_outlets.iloc[o,3] - z)*dz)**2)
            n_e_x = eta(s,nbrhd_a[o],Ca[o,0]) 
            chi = n_e_x*gamma_a/myu#*(1/e**3)*dx*dy*dz  
            ws.writerow([eqn_n, c_dom[y,x,z], chi])
            # data.append(chi) ; row.append(eqn_n) ; col.append(c_dom[y,x,z]) 
        
        eqn_n = eqn_n + 1
    
    return()
    
# veins
def terminal_veins():
    # row = [] ; col = [] ; data = []
    
    eqn_n = 2*Ntvx + (nNa - 2) + (nNv - 2 - len(v_outlets))
    for o in range(len(v_outlets)):
        i = v_outlets.iloc[o,0] 
        ele = v_element.loc[v_element['2'] == i]
        k = np.pi*(ele.iloc[0,3]*dx)**4/(128*myu*(ele.iloc[0,4]*dx))
        
        ws.writerow([eqn_n, 2*Ntvx + nNa +  ele.iloc[0,1], k])
        ws.writerow([eqn_n, 2*Ntvx + nNa + i, -(k + gamma_v/myu)])
        # data.append(k) ; row.append(eqn_n) ; col.append(2*Ntvx + nNa +  ele.iloc[0,1])
        # data.append(-(k + gamma_v/myu)) ; row.append(eqn_n) ; col.append(2*Ntvx + nNa + i)
        
        for j in range(len(nbrhd_v[o])):
            y,x,z = nbrhd_v[o][j]
            s = np.sqrt(((v_outlets.iloc[o,1] -  x)*dx)**2 + ((v_outlets.iloc[o,2] -  y)*dy)**2 + ((v_outlets.iloc[o,3] - z)*dz)**2)
            n_e_x = eta(s,nbrhd_v[o],Cv[o,0])
            chi = n_e_x*gamma_v/myu#*(1/e**3)*dx*dy*dz 
            ws.writerow([eqn_n, c_dom[y,x,z] + Ntvx, chi])
            # data.append(chi) ; row.append(eqn_n) ; col.append(c_dom[y,x,z] + Ntvx)
        
        eqn_n = eqn_n + 1
    
    return()

# Boundary conditions

def inlet_condition():
    # row = [] ; col = [] ; data = []
    
    eqn_n = 2*Ntvx + (nNa - 2) + (nNv - 2)
    
    ws.writerow([eqn_n, 2*Ntvx + a_inlet1, 1]) ; eqn_n = eqn_n + 1
    ws.writerow([eqn_n, 2*Ntvx + a_inlet2, 1]) ; eqn_n = eqn_n + 1
    ws.writerow([eqn_n, 2*Ntvx + nNa + v_inlet1, 1]) ; eqn_n = eqn_n + 1
    ws.writerow([eqn_n, 2*Ntvx + nNa + v_inlet2, 1])
    
    # row.append(eqn_n) ; col.append(2*Ntvx + a_inlet1) ; data.append(1) ; eqn_n = eqn_n + 1
    # row.append(eqn_n) ; col.append(2*Ntvx + a_inlet2) ; data.append(1) ; eqn_n = eqn_n + 1
    # row.append(eqn_n) ; col.append(2*Ntvx + nNa + v_inlet1) ; data.append(1) ; eqn_n = eqn_n + 1
    # row.append(eqn_n) ; col.append(2*Ntvx + nNa + v_inlet2) ; data.append(1) 
    
    return()
    
if __name__ == "__main__":
    
    xl_title = 'new_method_flow_rcd/other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.csv'
    wb = open(xl_title, 'w')
    ws = csv.writer(wb)
    ws.writerow(['row', 'col', 'data'])
    
    
    arterial_nodal_mass_balance()
    venal_nodal_mass_balance()
    terminal_arteries()
    terminal_veins()
    inlet_condition()
    
    
    wb.close()
    
    t2 = time.time()
    
    print('time taken for the entire code = ', round((t2-t1)/60,2), ' minutes')
    
    # r1 , c1, d1 = arterial_nodal_mass_balance()
    # r2, c2, d2 = venal_nodal_mass_balance()
    # r3, c3, d3 = terminal_arteries()
    # r4, c4, d4 = terminal_veins()
    # r5, c5, d5 = inlet_condition()
    
    # row = r1 + r2 + r3 + r4 + r5
    # col = c1 + c2 + c3 + c4 + c5
    # data = d1 + d2 + d3 + d4 + d5
    
    # r_title = 'new_method_flow_rcd/row_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy'
    # c_title = 'new_method_flow_rcd/col_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy'
    # d_title = 'new_method_flow_rcd/data_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy'
    
    # np.save(r_title,row)
    # np.save(c_title,col)
    # np.save(d_title,data)
    
    # t2 = time.time()
    print('Other Equations'+' time = ',round((t2-t1)/60,3),' min')
    
    # print('max of row = ',max(row))
    # print('max of col = ',max(col))

 
