import numpy as np
import pandas as pd
from multiprocessing import Pool, Manager
import time
import Flow_solver_parameter_file_loader as para
import scipy.sparse as sp
import matplotlib.pyplot as plt
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

t2 = time.time()
print('time taken to load all the files = ', t2-t1)
# Number of unknowns
# N_unk = 2*Ntvx + nNa + nNv

# Matrix generation

def eta(x,nbrhd,C):
    global e
    if(x/e < 1):
        eta_x = C*np.exp(1/((abs(x/e))**2 -1))#/(e**3)#*(dx*dy*dz) 
    else:
        eta_x = 0
    
    return eta_x #1/len(nbrhd)



def func(input_value):
    nx,ny_array,nz,ws = input_value
    ny1 = ny_array[0]
    ny2 = ny_array[1]
    
    # row = [] ; col = [] ; data = []
    # row = 0
    for k in range(nz):
        for j in range(ny1,ny2):
            for i in range(nx):
                # ta = time.time()
                
                if(dom[j,i,k] == 1):
                    eqn_n = Ntvx + c_dom[j,i,k]
                # Veins only
                    ij_sum = 0
                    if(dom[j,i-1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j,i-1,k], -tx])
                        
                        ij_sum = ij_sum + tx
                    if(dom[j,i+1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j,i+1,k], -tx])
                        
                    
                        ij_sum = ij_sum + tx
                    if(dom[j-1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j-1,i,k], -ty])
                        
                         
                        ij_sum = ij_sum + ty
                    if(dom[j+1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j+1,i,k], -ty])
                        
                        ij_sum = ij_sum + ty
                        
                    if(dom[j,i,k+1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j,i,k+1], -tz])
                        
                        ij_sum = ij_sum + tz
                    
                    if(dom[j,i,k-1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        
                        ws.writerow([eqn_n, Ntvx + c_dom[j,i,k-1], -tz])
                        
                        ij_sum = ij_sum + tz
                    
                    ij_sum = ij_sum + a*dx*dy*dz
                    
                    ws.writerow([eqn_n, Ntvx + c_dom[j,i,k], ij_sum])
                    ws.writerow([eqn_n, c_dom[j,i,k], -a*dx*dy*dz])
                    
                    
                    for b in range(len(v_outlets)):
                        if([j,i,k] in nbrhd_v[b]):
                            ele = v_element.loc[v_element['2']==v_outlets.iloc[b,0]]
                            k1 = np.pi*(ele.iloc[0,3]*dx)**4/(128*myu*(ele.iloc[0,4]*dx))
                            s = np.sqrt(((i-v_outlets.iloc[b,1])*dx)**2 + ((j-v_outlets.iloc[b,2])*dy)**2 + ((k-v_outlets.iloc[b,3])*dz)**2)
                            n_e_x = eta(s,nbrhd_v[b],Cv[b,0])
                            G = k1*n_e_x#*(1/e**3)*dx*dy*dz 
                            
                            ws.writerow([eqn_n, 2*Ntvx + nNa + ele.iloc[0,1], -G])
                            ws.writerow([eqn_n, 2*Ntvx + nNa + ele.iloc[0,2], G])
                            
                            
    # tb = time.time()
    # print('time taken to complete one iteration = ', tb-ta)                
    return(0)



if __name__ == "__main__":
    
    xl_title = 'new_method_flow_rcd/venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.csv'
    wb = open(xl_title, 'w')
    writer = csv.writer(wb)
    writer.writerow(['row', 'col', 'data'])
    
    # ny = 20
    test  = func([nx,[0,ny],nz,writer])
    
    wb.close()
    
    # rcd = []
    # dny = 100
    # pool_n = 100
    # ny_list = []
    # ny = int(ny/2)
    # for i in range(dny):
    #     ny_list.append([int(i*ny/dny),int((i+1)*ny/dny)])
    
    # p = Pool(pool_n)
    # input_value = [(nx,ny_list[x],nz) for x in range(dny)]
    # rcd.append(p.map(func,input_value))
    # p.close()
    # p.join()
    
    # row = []
    # col = []
    # data = []
    # N_unk = para.N_unk
    
    # for i in range(len(rcd[0])):
    #     row = row + rcd[0][i][0]
    #     col = col + rcd[0][i][1]
    #     data = data + rcd[0][i][2]
        
    # A = sp.csc_matrix((data,(row,col)), shape=(N_unk,N_unk))
    
    # plt.figure(figsize=(10,10))
    # plt.spy(A)
    # plt.savefig('venal_compartment_e_'+str(e)+'_dny_' + str(dny) + '_pool_' + str(pool_n) + '.png')
    
    
    # ra_title = 'new_method_flow_rcd/row_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_1_.npy'
    # ca_title = 'new_method_flow_rcd/col_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_1_.npy'
    # da_title = 'new_method_flow_rcd/data_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_1_.npy'
    
    # np.save(ra_title,row)
    # np.save(ca_title,col)
    # np.save(da_title,data)

    t2 = time.time()
    print('Time taken for the code = ', round((t2 - t1)/60,2), ' minutes')
    # print('Venal Compartment\n pool number = '+ str(pool_n)+' dny = '+str(dny)+' time = ',round((t2-t1)/60,3),' min')
    
    # print('max of row = ',max(row))
    # print('max of col = ',max(col))
    
    