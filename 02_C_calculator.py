import numpy as np
import pandas as pd
import neighbourhood_matrix as nbm
from concurrent.futures import ProcessPoolExecutor

old = False

a_element = pd.read_csv('arteries_element_database.csv')
v_element = pd.read_csv('veins_element_database.csv')

a_outlets = pd.read_csv('arteries_outlet_coordinates_3D_shifted.csv')
v_outlets = pd.read_csv('veins_outlet_coordinates_3D_shifted.csv')

e =  6.4e-05
dx = nbm.dx
dy = nbm.dy
dz = nbm.dz

if(old == True):
    title_nbr_a = 'nbrhd_matrices/'+str(e)+'/nbrhd_old_a_dx_dy_'+str(dx)+'_e_' + str(e) + '.npy'
    title_nbr_v = 'nbrhd_matrices/'+str(e)+'/nbrhd_old_v_dx_dy_'+str(dx)+'_e_' + str(e) + '.npy'
    
if(old == False):
    title_nbr_a = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_a_dx_dy_'+str(dx)+'_e_' + str(e) + '_new.npy'
    title_nbr_v = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_v_dx_dy_'+str(dx)+'_e_' + str(e) + '_new.npy'

nbrhd_a = np.load(title_nbr_a, allow_pickle=True).tolist()
nbrhd_v = np.load(title_nbr_v, allow_pickle=True).tolist()

def artery():
    global dx
    global dy
    global dz
    global e
    global a_outlets
    global nbrhd_a
    
    C_a = []
    
    for i in range(len(a_outlets)):
        c = []
        C_a.append(c)
    
    for i in range(len(a_outlets)):
        sum_exp = 0.0
        for j in range(len(nbrhd_a[i])):
            s = np.sqrt( ((nbrhd_a[i][j][1] - a_outlets.iloc[i,1])*dx)**2 + ((nbrhd_a[i][j][0] - a_outlets.iloc[i,2])*dy)**2 + ((nbrhd_a[i][j][2] - a_outlets.iloc[i,3])*dz)**2) 
            if(s/e < 1):
                exp = np.exp(1/(abs(s/e)**2 - 1))/(e**3) #;print(s,s/e)
                sum_exp = sum_exp + exp ; #print(exp)
        C_a[i].append(1/sum_exp)
    
    if(old == True):
        ca_title = 'constants/Ca_old_dx_dy_'+str(dx) +'_e_' + str(e) + '.npy' 
    if(old == False):
        ca_title = 'constants/Ca_3D_dx_dy_'+str(dx) +'_e_' + str(e) + '.npy' 
    np.save(ca_title,np.array(C_a))
    

def veins():
    global dx
    global dy
    global dz
    global e
    global v_outlets
    global nbrhd_v
    
    C_v = []
    
    for j in range(len(v_outlets)):
        c = []
        C_v.append(c)
        
    for i in range(len(v_outlets)):
        sum_exp = 0.0
        for j in range(len(nbrhd_v[i])):
            s = np.sqrt( ((nbrhd_v[i][j][1] - v_outlets.iloc[i,1])*dx)**2 + ((nbrhd_v[i][j][0] - v_outlets.iloc[i,2])*dy)**2 + ((nbrhd_v[i][j][2] - v_outlets.iloc[i,3])*dz)**2) 
            if(s/e < 1):
                exp = np.exp(1/(abs(s/e)**2 - 1))/(e**3) ;#print(s/e)
                sum_exp = sum_exp + exp ; #print(exp)
        C_v[i].append(1/sum_exp)
    
    if(old == True):
        cv_title = 'constants/Cv_old_dx_dy_'+str(dx) +'_e_' + str(e) + '.npy' 
    if(old == False):
        cv_title = 'constants/Cv_3D_dx_dy_'+str(dx) +'_e_' + str(e) + '.npy' 
    np.save(cv_title,np.array(C_v))

def main():
    with ProcessPoolExecutor(max_workers=3) as executor_p:
        executor_p.submit(artery)
        executor_p.submit(veins)
        
if __name__=='__main__':
    main()