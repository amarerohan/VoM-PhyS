import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np
import Flow_solver_parameter_file_loader as para
import time
import os
import matplotlib.pyplot as plt
import pandas as pd

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
nNa = para.nNa
nNv = para.nNv
N_unk = para.N_unk


def eta(x,nbrhd,C):
    global e
    if(x/e < 1):
        eta_x = C*np.exp(1/((abs(x/e))**2 -1))#/(e**3)#*(dx*dy*dz) 
    else:
        eta_x = 0
    
    return eta_x 



if __name__ == "__main__":
    
    f = open('Time_analysis_' + sensitivity_para + '.txt','w')
    
    rcd_a = pd.read_csv('new_method_flow_rcd/arterial_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.csv')
    rcd_v = pd.read_csv('new_method_flow_rcd/venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.csv')
    rcd_o = pd.read_csv('new_method_flow_rcd/other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.csv')
    
    # np.save('A.npy',np.arange(0,100,1))
    
    r = rcd_a.iloc[:,0].tolist() + rcd_v.iloc[:,0].tolist() + rcd_o.iloc[:,0].tolist()
    c = rcd_a.iloc[:,1].tolist() + rcd_v.iloc[:,1].tolist() + rcd_o.iloc[:,1].tolist()
    d = rcd_a.iloc[:,2].tolist() + rcd_v.iloc[:,2].tolist() + rcd_o.iloc[:,2].tolist()
    
    # np.save('B.npy',np.arange(0,100,1))
    
    '''
    row1 = np.load('new_method_flow_rcd/row_arterial_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    col1 = np.load('new_method_flow_rcd/col_arterial_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    data1 = np.load('new_method_flow_rcd/data_arterial_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    
    row2 = np.load('new_method_flow_rcd/row_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    col2 = np.load('new_method_flow_rcd/col_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    data2 = np.load('new_method_flow_rcd/data_venal_compartment_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    
    row3 = np.load('new_method_flow_rcd/row_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    col3 = np.load('new_method_flow_rcd/col_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    data3 = np.load('new_method_flow_rcd/data_other_equations_e_'+str(e)+'_dx_'+str(dx)+'_sensitivity_'+sensitivity_para+'_.npy').tolist()
    
    r = row1 + row2 + row3
    c = col1 + col2 + col3
    d = data1 + data2 + data3
    '''
    path = 'new_method_flow_rcd/' + str(E)
    try:
        os.mkdir(path)
    except OSError as error:
        print('Directory already present')
        
    np.save(path + '/row_sensitivity_'+sensitivity_para+'.npy' ,r)
    np.save(path + '/col_sensitivty_'+sensitivity_para+'.npy' ,c)
    np.save(path + '/data_sensitivty_'+sensitivity_para+'.npy',d)
    
    ta = time.time()
    
    # r = np.load(path + '/row_sensitivity_'+sensitivity_para+'.npy')
    # c = np.load(path + '/col_sensitivity_'+sensitivity_para+'.npy')
    # d = np.load(path + '/data_sensitivity_'+sensitivity_para+'.npy')
    # np.save('C.npy',np.arange(0,100,1))
    print(max(r),N_unk)
    
    A = sp.csc_matrix((d,(r,c)), shape=(N_unk,N_unk))
    
    plt.figure(figsize=(10,10))
    plt.spy(A)
    plt.title('Flow Simulation ' + str(e))
    plt.savefig(path + '/sparse_matrix_'+sensitivity_para + '.png')
    
    B = np.zeros((N_unk,1), dtype = float)
    B[-4] = B[-3] = P_ba
    B[-2] = P_bv
    B[-1] = P_bv
    
    t2 = time.time()
    print('e = ',e,' time needed to load and create matrix = ',round((t2-t1)/60,3),' mins')
    f.write(' Time taken to load and create matrix A and b \n')
    f.write(str(round((t2-ta)/60,2)))
    
    e_ref = 0.001
    if(e > 0.001):
        r0 = np.load('new_method_flow_rcd/'+str(0.001)+'/row.npy')
        c0 = np.load('new_method_flow_rcd/'+str(0.001)+'/col.npy')
        d0 = np.load('new_method_flow_rcd/'+str(0.001)+'/data.npy')
        
        A0 = sp.csc_matrix((d0,(r0,c0)), shape=(N_unk,N_unk))
        LU = spla.splu(A0)
    
    elif(e <= 0.001):
        LU = spla.splu(A)
    
    X0 = np.load('new_method_flow_rcd/6.4e-05/flow_solution.npy')
    # X0 = np.copy(B) ; X[:] = P_bv
    
    M = spla.LinearOperator(np.shape(LU),LU.solve)
    
    tm = time.time()
    f.write(' Time taken to create M \n')
    f.write(str(round((tm-t2)/60,2)))
    
    # np.save('M_done.npy',np.arange(0,100,1))
    
    tol = 1e-8
    X = spla.gmres(A,B,M=M,x0=X0,tol=tol)
    print('error code = ',X[1])
    X = X[0]
    
    t3 = time.time()
    print('time required to solve the matrix = ',round((t3-t2)/60,3),' mins')
    
    f.write(' time required to solve the matrix \n')
    f.write(str(round((t3-tm)/60,2)))
    f.close()
    
    X_title = path + '/flow_solution_new_X0_6.4e-05_sensitivity_'+sensitivity_para+'.npy'
    np.save(X_title,X)
    
    print(np.max(X[:]),np.min(X[:]))
    
    np.save('X_done.npy',np.arange(0,100,1))
    
    # # #  Solver over Cross Verification and saving data begins # # # 
    
    t4 = time.time()
    
    a_p = np.zeros((ny,nx,nz),dtype = float)
    a_p[:,:,:] = 0
    counter = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(dom[j,i,k] == 1):
                    a_p[j,i,k] = X[counter]
                    counter = counter + 1

    aP_min = min(X[0:Ntvx])
    aP_max = max(X[0:Ntvx])


    v_p = np.zeros((ny,nx,nz),dtype = float)
    v_p[:,:,:] = 0
    counter = Ntvx
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(dom[j,i,k]==1):
                    v_p[j,i,k] = X[counter]
                    counter = counter + 1

    vP_max = max(X[Ntvx:2*Ntvx])
    vP_min = min(X[Ntvx:2*Ntvx])


    Pa = X[2*Ntvx:2*Ntvx + nNa]
    Pv = X[2*Ntvx + nNa:N_unk]


    # # # # Flow calculations # # # # 
    KA = []
    QA = []
    DPA = []
    for i in range(len(a_element)):
        phi = np.pi*(a_element.iloc[i,3]*dx)**4/(128*myu*a_element.iloc[i,4]*dx)
        KA.append(phi)
        delta_P = Pa[a_element.iloc[i,1]] - Pa[a_element.iloc[i,2]]
        DPA.append(delta_P)
        q = delta_P*phi
        QA.append(q)

    KV = []
    QV = []
    DPV = []
    for i in range(len(v_element)):
        phi = np.pi*(v_element.iloc[i,3]*dx)**4/(128*myu*v_element.iloc[i,4]*dx)
        KV.append(phi)
        delta_P = Pv[v_element.iloc[i,1]] - Pv[v_element.iloc[i,2]]
        DPV.append(delta_P)
        q = delta_P*phi
        QV.append(q)


    np.save(path + '/Arterial_Compartment_Pressure_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',a_p)
    np.save(path + '/Venal_Compartment_Pressure_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',v_p)
    np.save(path + '/Arterial_Nodal_Pressure_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',Pa)
    np.save(path + '/Venal_Nodal_Pressure_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',Pv)
    np.save(path + '/Arterial_Elemental_Flow_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',QA)
    np.save(path + '/Venal_Elemental_Flow_' +str(E)+'_sensitivity_'+sensitivity_para+'_.npy',QV)


    DP = a_p - v_p

    Qin = abs(QA[a_element.loc[a_element['1'] == a_inlet1].iloc[0,0]]) + abs(QA[a_element.loc[a_element['2'] == a_inlet2].iloc[0,0]])
    Qout = abs(QV[v_element.loc[v_element['1'] == v_inlet1].iloc[0,0]]) + abs(QV[v_element.loc[v_element['2'] == v_inlet2].iloc[0,0]])

    error_mass = Qin - Qout


    QA_flow = np.zeros((ny,nx,nz),dtype = float)
    QV_flow = np.copy(QA_flow)
    Q_perf  = np.copy(QA_flow)
    QA_diff = np.copy(QA_flow)
    QV_diff = np.copy(QA_flow)

    QA_percentage = QA[:]/Qin*100
    QV_percentage = QV[:]/Qin*100

    for j in range(len(nbrhd_a)):
        ele = a_element.loc[a_element['2'] == a_outlets.iloc[j,0]]
        for i in range(len(nbrhd_a[j])):
           y,x,z = nbrhd_a[j][i]
           s = np.sqrt(((x-a_outlets.iloc[j,1])*dx)**2 + ((y-a_outlets.iloc[j,2])*dy)**2  + ((a_outlets.iloc[j,3] - z)*dz)**2)
           n_e_x = eta(s,nbrhd_a[j],Ca[j,0])
           QA_flow[y,x,z] = QA_flow[y,x,z] + abs(QA[ele.iloc[0,0]])*n_e_x#*(1/e**3)*dx*dy*dz



    for j in range(len(nbrhd_v)):
        ele = v_element.loc[v_element['2'] == v_outlets.iloc[j,0]]
        for i in range(len(nbrhd_v[j])):
           y,x,z = nbrhd_v[j][i]
           s = np.sqrt(((x-v_outlets.iloc[j,1])*dx)**2 + ((y-v_outlets.iloc[j,2])*dy)**2  + ((v_outlets.iloc[j,3] - z)*dz)**2)
           n_e_x = eta(s,nbrhd_v[j],Cv[j,0])
           QV_flow[y,x,z] = QV_flow[y,x,z] + abs(QV[ele.iloc[0,0]])*n_e_x#*(1/e**3)*dx*dy*dz




    Q_perf = a*dx*dy*dz*(DP[:,:,:])



    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(dom[j,i,k]==1):
                    ij_sum = 0
                    if(dom[j,i-1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_a + dx/Lambda_a))
                        ij_sum = ij_sum + tx*(a_p[j,i-1,k] - a_p[j,i,k])

                    if(dom[j,i+1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_a + dx/Lambda_a))
                        ij_sum = ij_sum + tx*(a_p[j,i+1,k] - a_p[j,i,k])

                    if(dom[j-1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_a + dy/Lambda_a))
                        ij_sum = ij_sum + ty*(a_p[j-1,i,k] - a_p[j,i,k])

                    if(dom[j+1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_a + dy/Lambda_a))
                        ij_sum = ij_sum + ty*(a_p[j+1,i,k] - a_p[j,i,k])

                    if(dom[j,i,k+1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_a + dz/Lambda_a))
                        ij_sum = ij_sum + tz*(a_p[j,i,k+1] - a_p[j,i,k])

                    if(dom[j,i,k-1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_a + dz/Lambda_a))
                        ij_sum = ij_sum + tz*(a_p[j,i,k-1] - a_p[j,i,k])


                    QA_diff[j,i,k] = ij_sum

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if(dom[j,i,k]==1):                
                    ij_sum = 0
                    if(dom[j,i-1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        ij_sum = ij_sum + tx*(v_p[j,i-1,k] - v_p[j,i,k])

                    if(dom[j,i+1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        ij_sum = ij_sum + tx*(v_p[j,i+1,k] - v_p[j,i,k])

                    if(dom[j-1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        ij_sum = ij_sum + ty*(v_p[j-1,i,k] - v_p[j,i,k])

                    if(dom[j+1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        ij_sum = ij_sum + ty*(v_p[j+1,i,k] - v_p[j,i,k])

                    if(dom[j,i,k+1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        ij_sum = ij_sum + tz*(v_p[j,i,k+1] - v_p[j,i,k])

                    if(dom[j,i,k-1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        ij_sum = ij_sum + tz*(v_p[j,i,k-1] - v_p[j,i,k])

                    QV_diff[j,i,k] = ij_sum



    # Conservation of mass within arteries
    QA_out = 0.0
    for i in range(len(a_outlets)):
        ele = a_element.loc[a_element['2']==a_outlets.iloc[i,0]]
        QA_out = QA_out + abs(QA[ele.iloc[0,0]])
    QA_error = (abs(Qin) - abs(QA_out))/(min(abs(Qin),abs(QA_out)))*100

    QV_in = 0.0
    for i in range(len(v_outlets)):
        ele = v_element.loc[v_element['2'] == v_outlets.iloc[i,0]]
        QV_in = QV_in + abs(QV[ele.iloc[0,0]])
    QV_error = (abs(Qout) - abs(QV_in))/(min(abs(Qout),abs(QV_in)))*100

    Q_tissue_error = (abs(QA_out) - abs(QV_in))/(min(abs(QA_out),abs(QV_in)))*100


    QA_A_flow = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                QA_A_flow = QA_A_flow + QA_flow[j,i,k]

    QV_V_flow = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                QV_V_flow = QV_V_flow + QV_flow[j,i,k]

    Q_perfusion_net = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                Q_perfusion_net = Q_perfusion_net + Q_perf[j,i,k]

    print('net_error                    ',round(abs(error_mass)/min(Qin,Qout)*100,5),'%')
    print('error Q within arteries      ', round(QA_error,5),'%')
    print('error Q within veins         ', round(QV_error,5),'%')
    print('error Q QA_out QV_in         ', round(Q_tissue_error,5),'%')
    print('error Q QA_out Qtissue_in    ', round((QA_out - QA_A_flow)/(min(abs(QA_out),abs(QA_A_flow))),10),'%')
    print('\n')
    print('Q entering                   ', round(Qin,10))
    print('Q coming out of artery       ', round(QA_out,10))
    print('Q artery to tissue           ', round(QA_A_flow,10))
    print('Perf Q within tissue         ', round(Q_perfusion_net,10))
    print('Q tissue to vein             ', round(QV_V_flow,10))
    print('Q entering vein              ', round(QV_in,10))
    print('Q coming out of vein         ', round(Qout,10))

    # Q_A_compare = np.zeros((len(a_outlets),3),dtype = float)
    # for i in range(len(a_outlets)):
    #     ele = a_element.loc[a_element['2']==a_outlets.iloc[i,0]] #; print(ele.iloc[0,0])
    #     Q_A_compare[i,0] = abs(QA[ele.iloc[0,0]])
    #     Q_A_compare[i,1] = abs(QA_flow[nbrhd_a[i][0][0],nbrhd_a[i][0][1]])
    #     Q_A_compare[i,2] = abs(Q_perf[nbrhd_a[i][0][0],nbrhd_a[i][0][1]])


    # Plotting Contours
    
    cl = 100

    for i in range(1,4):
        plt.figure(figsize=(10,12))
        plt.contourf(np.flip(a_p[:,:,i],0),np.arange(aP_min,aP_max*(1.001),(aP_max*(1.001) - aP_min)/cl))
        plt.colorbar().ax.tick_params(labelsize = 18)
        plt.title('Arterial Compartment pressure'+ ' | dx = dy = ' + str(dx) + ' m ' + '| e = ' + str(E) + 'm' + ' layer ' + str(i))
        titlesave = path + '/Arterial_Compartment_Pressure_dx_'+str(dx)+'_e_'+str(E)+'_3D_layer_'+str(i)+'_sensitivity_'+sensitivity_para+'.png'
        plt.savefig(titlesave,dpi = 1000)

    for i in range(1,4):
        plt.figure(figsize=(10,12))
        plt.contourf(np.flip(v_p[:,:,i],0),np.arange(vP_min, vP_max*(1 + 1/cl),(vP_max*(1+1/cl) - vP_min)/cl))
        plt.colorbar().ax.tick_params(labelsize = 18)
        plt.title('Venal Compartment Pressure' + ' | dx = dy = ' + str(dx) + ' m ' + '| e = ' + str(e) + 'm'+ ' layer ' + str(i))
        titlesave = path + '/Venal_Compartment_Pressure_dx_'+str(dx)+'_e_'+str(E)+'_3D_layer_'+str(i)+'_sensitivity_'+sensitivity_para+'.png'
        plt.savefig(titlesave,dpi = 1000)
           
    print('\n PARAMETERS USED \n')
    print('Ka = ', Ka)
    print('Kv = ', Kv)
    print('alpha = ',a)
    q_df = pd.DataFrame(QV)
    q = q_df.loc[q_df[0] > 0]
    print('\n q \n', len(q), '\n', q)
    
    
    
    '''
    a_mass_cons = np.zeros((ny,nx,nz),dtype = float)
    v_mass_cons = np.copy(a_mass_cons)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                #print(j,i,k)
                if(dom[j,i,k] == 1):# and (dom[j,i,0] != 1 and dom[j,i,2] != 1)):
                # Arteries only
                    ij_sum = 0
                    if(dom[j,i-1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_a + dx/Lambda_a))
                        
                        ij_sum = ij_sum + tx*(a_p[j,i-1,k] - a_p[j,i,k])
                        
                    if(dom[j,i+1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_a + dx/Lambda_a))
                        
                        ij_sum = ij_sum + tx*(a_p[j,i+1,k] - a_p[j,i,k])
                        
                    if(dom[j-1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_a + dy/Lambda_a))
                        
                        ij_sum = ij_sum + ty*(a_p[j-1,i,k] - a_p[j,i,k])
                    
                    if(dom[j+1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_a + dy/Lambda_a))
                        
                        ij_sum = ij_sum + ty*(a_p[j+1,i,k] - a_p[j,i,k])
                    
                    if(dom[j,i,k+1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_a + dz/Lambda_a))
                        ij_sum = ij_sum + tz*(a_p[j,i,k+1] - a_p[j,i,k])
                    
                    if(dom[j,i,k-1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_a + dz/Lambda_a))
                        ij_sum = ij_sum + tz*(a_p[j,i,k-1] - a_p[j,i,k])
                    
                    ij_sum = ij_sum - a*dx*dy*dz*(a_p[j,i,k] - v_p[j,i,k])
                                        
                    for b in range(len(a_outlets)):
                        if([j,i,k] in nbrhd_a[b]):
                            # print(i,j,b)
                            ele = a_element.loc[a_element['2']==a_outlets.iloc[b,0]]
                            
                            k1 = np.pi*(ele.iloc[0,3]*dx)**4/(128*myu*(ele.iloc[0,4])*dx)
                            s = np.sqrt(((i-a_outlets.iloc[b,1])*dx)**2 + ((j-a_outlets.iloc[b,2])*dy)**2 + ((k-a_outlets.iloc[b,3])*dz)**2) #; print(s,[j,i,k])
                            n_e_x = eta(s,nbrhd_a[b],Ca[b,0]) 
                            G = abs(QA[ele.iloc[0,0]])*n_e_x
                            ij_sum = ij_sum + G
                            
                            
                    a_mass_cons[j,i,k] = ij_sum
    
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):  
                if(dom[j,i,k] == 1):
                    
                    ij_sum = 0
                    if(dom[j,i-1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        
                        ij_sum = ij_sum + tx*(v_p[j,i-1,k] - v_p[j,i,k])
                    if(dom[j,i+1,k] == 1):
                        tx = 2*dy*dz*(1/(dx/Lambda_v + dx/Lambda_v))
                        
                        ij_sum = ij_sum + tx*(v_p[j,i+1,k] - v_p[j,i,k])
                    if(dom[j-1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        
                        ij_sum = ij_sum + ty*(v_p[j-1,i,k] - v_p[j,i,k])
                    if(dom[j+1,i,k] == 1):
                        ty = 2*dx*dz*(1/(dy/Lambda_v + dy/Lambda_v))
                        
                        ij_sum = ij_sum + ty*(v_p[j+1,i,k] - v_p[j,i,k])
                        
                    if(dom[j,i,k+1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        
                        ij_sum = ij_sum + tz*(v_p[j,i,k+1] - v_p[j,i,k])
                    
                    if(dom[j,i,k-1] == 1):
                        tz = 2*dx*dy*(1/(dz/Lambda_v + dz/Lambda_v))
                        
                        ij_sum = ij_sum + tz*(v_p[j,i,k-1] - v_p[j,i,k])
                    
                    ij_sum = ij_sum + a*dx*dy*dz*(a_p[j,i,k] - v_p[j,i,k])
                    
                    
                    for b in range(len(v_outlets)):
                        if([j,i,k] in nbrhd_v[b]):
                            ele = v_element.loc[v_element['2']==v_outlets.iloc[b,0]]
                            k1 = np.pi*(ele.iloc[0,3]*dx)**4/(128*myu*(ele.iloc[0,4]*dx))
                            s = np.sqrt(((i-v_outlets.iloc[b,1])*dx)**2 + ((j-v_outlets.iloc[b,2])*dy)**2 + ((k-v_outlets.iloc[b,3])*dz)**2)
                            n_e_x = eta(s,nbrhd_v[b],Cv[b,0])
                            G = abs(QV[ele.iloc[0,0]])*n_e_x
                            
                            ij_sum = ij_sum - G
                    
                    v_mass_cons[j,i,k] = ij_sum
                    
    np.save(str(e)+'a_mass_cons.npy',a_mass_cons)
    np.save(str(e)+'v_mass_cons.npy',v_mass_cons)
    print('a_cons_error = ',np.max(np.abs(a_mass_cons[:,:,:])))
    print('v_cons_error = ',np.max(np.abs(v_mass_cons[:,:,:])))


    t5 = time.time()
    
    
    print('Time needed for saving data and cross_verification = ', round((t5-t4)/60,3),' mins')
    '''
    
    
    
    
    
    
    
    