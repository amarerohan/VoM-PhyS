# VoM-PhyS
VoM-PhyS is a simulation framework to simulate blood flow and heat transfer in biological domains. The repository shared here, consists of codes and example domain used for simulation.

Blood vessels segmented from the imaging data are modeled as 1D flow using Hagen-Poiseuille equation. The 3D tissue is modeled as porous domain, and 1D - 3D flow domain is coupled using Dirac distribution function.

The framework and code free for users to use and modify, but we request you to cite us. 

The paper associated with this work can be found at

https://www.nature.com/articles/s41598-022-18831-3


STEPS TO RUN THE CODES:

FLOW SOLVER 

1. Load the neighbourhood_matrix.py and specify the parameters of voxel resolution and the radius of SoI (e). The code generates two .npy files which will be required later.
2. Load and run the C_calculator.py. Ensure that the the parameters of SoI radius (e) and voxel resolution are same as the neighbourhood_matrix.py. The code will generate two .npy files which will be required later
3. Matrix Generation: This step is in 4 parts
A. Load Flow_solver_parameter_file_loader.py: For the frog tongue domain provided, the inlet and outlets will not change. The parameters like SoI radius (e) and voxel resolutions should match to the previous two steps. As this ensures correct .npy files are being loaded.
B. Flow_arterial_compartment_read_write.py - Run this code and it will generate an excel sheet that stores the row, col and data of the flow matrix related to the flow equations of the arterial compartment
C. Flow_venal_compartment_read_write.py - Run this code and it will generate an excel sheet that stores the row, col and data of the flow matrix related to the flow equations of the venous compartment
D. Flow_Arterial_Venal_Other_Equations_read_write.py - Run this code and it will generate an excel sheet that stores the row, col and data of the flow matrix related to the segmented arteris and veins
5. Flow_load_and_solve_new_method.py The final step is to load all the excel sheets and solve the flow matrix. The error should be 0 for mass conservation and depends on the iterations and tolerance used for the GMRES solver. In addition, a preconditioner may be needed depending on the SoI radius used. The code should plot the pressure maps when the solution is correctly solved.



HEAT SOLVER COUPLED WITH FLOW
This section will only work if the flow solver part has been succesfully executed.

1. Heat_solver_parameter_file_loader.py - Provide the parameters and ensure they are correct with respect to the flow solver section. This will ensure correct files are being loaded.
2. Matrix generation: This step is in two parts
   A. Heat_voxel_equations_read_write.py - Run this code and it will generate an excel sheet that stores the row, col and data of the heat equations for the tissue voxels.
   B. Heat_artery_vein_equations_read_write.py - Run this code and it will generate an excel sheet that stores the row, col and data of the heat equations for the segmented vessels.
3. Heat_load_and_solver_read_write.py - The final step to solve for the heat transfer equation.
