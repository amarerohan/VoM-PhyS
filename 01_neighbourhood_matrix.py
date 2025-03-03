import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import os 


dx = dy = 6.4e-5
dz = 0.333e-3

e = 6.4e-5
old = False

a_element = pd.read_csv('arteries_element_database.csv')
v_element = pd.read_csv('veins_element_database.csv')

a_outlets = pd.read_csv('arteries_outlet_coordinates_3D_shifted.csv')
v_outlets = pd.read_csv('veins_outlet_coordinates_3D_shifted.csv')


dom = np.load('tongue_3D.npy')

ny, nx, nz = np.shape(dom)

# Number of elements in the arterial blood vessel network
nP = len(a_element)
nN = int(max((max(a_element.iloc[:,1]),max(a_element.iloc[:,2])))) + 1

# Number of elements in the venal blood vessel network
nPv = len(v_element)
nNv = int(max((max(v_element.iloc[:,1]),max(v_element.iloc[:,2])))) + 1


# Creating neighbourhood matrices
nbr_a = []
for n in range(len(a_outlets)):
    nbr = []
    nbr_a.append(nbr)
    
nbr_v = []
for n in range(len(v_outlets)):
    nbr = []
    nbr_v.append(nbr)


# Atempting multiththreading
def NBRHD_A(x,y,z,kn):
    global a_outlets
    global e
    global nbr_a
    
    s = np.sqrt( ((x - a_outlets.iloc[kn,1])*dx)**2 + ((y - a_outlets.iloc[kn,2])*dy)**2 +((z - a_outlets.iloc[kn,3])*dz)**2) 
    if(s < e):
        nbr_a[kn].append([y,x,z]) 
    
        
def NBRHD_V(x,y,z,kn):
    global v_outlets
    global e
    global nbrhd_v
    s = np.sqrt( ((x - v_outlets.iloc[kn,1])*dx)**2 + ((y - v_outlets.iloc[kn,2])*dy)**2 +((z - v_outlets.iloc[kn,3])*dz)**2) 
    if(s < e):
        nbr_v[kn].append([y,x,z])

def artery():
    global nz
    global ny
    global nx
    global dom
    global a_outlets
    
    executor = ThreadPoolExecutor(max_workers=100)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):   
                if(dom[j,i,k] == 1):
                    for kn in range(len(a_outlets)):
                        arr = [i,j,k,kn]
                        executor.submit(lambda p: NBRHD_A(*p),arr)
    artery_nbrhd = np.array(nbr_a, dtype = object)
    
    path = 'nbrhd_matrices/' + str(e)
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)
    
    if(old == True):
        title_a = 'nbrhd_matrices/'+str(e)+'/nbrhd_old_a_dx_dy_'+str(dx)+'_e_' + str(e) + '.npy'
    if(old == False):
        title_a = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_a_dx_dy_'+str(dx)+'_e_' + str(e) + '_new.npy'
    np.save(title_a,artery_nbrhd)
    print(e,np.shape(artery_nbrhd))

def veins():
    global nz
    global ny
    global nx
    global dom
    global v_outlets
    
    executor = ThreadPoolExecutor(max_workers=100)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):    
                if(dom[j,i,k] == 1):
                    for kn in range(len(v_outlets)):
                        arr = [i,j,k,kn] 
                        executor.submit(lambda p: NBRHD_V(*p),arr)
    vein_nbrhd = np.array(nbr_v, dtype = object)
    
    path = 'nbrhd_matrices/' + str(e)
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)
        
    if(old == True):
        title_v = 'nbrhd_matrices/'+str(e)+'/nbrhd_old_v_dx_dy_'+str(dx)+'_e_' + str(e) + '.npy'
    if(old == False):
        title_v = 'nbrhd_matrices/'+str(e)+'/nbrhd_3D_v_dx_dy_'+str(dx)+'_e_' + str(e) + '_new.npy'
    np.save(title_v,vein_nbrhd)
    print(e,np.shape(vein_nbrhd))
    
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=2) as executor_p:
        executor_p.submit(artery)
        executor_p.submit(veins)
    
    
    
    
    
    

