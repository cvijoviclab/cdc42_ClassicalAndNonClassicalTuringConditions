#==============================================================================================
# FEM-SIMULATION OF  Cdc42-POLARISATION
#==============================================================================================
# DATE:
# 2019-11-14
# WRITTEN BY:
# Johannes Borgqvist and Carl Lundholm
# DESCRIPTION:
# The program simulates the "bulk-surface" RD-model of cell polarisation where we have
#  a coupled system of three ODEs describing the dynamics of three states: GTP-bound or active
# Cdc42 u, GDP-bound inactive Cdc42 v and GDI-bound cytosolic Cdc42 V. The first two states live
#  on the membrane which is why we use the surface measure "ds" in the variational formulation
# while the cytosolic species V exists within the cell and therefore has the corresponding measure
# "dx". We solve the PDE by using the Finite Element Method (FEM) in space and
# Finite Differences (FD) in time. We use a "Implicit-Explicit" solution algorithm where the
# we use the explicit Euler algorithm in time (FD) and where the non-linearities in the reaction
# term are evaluated at the previous time point which boils down to solving a linear problem
# in space (FEM). The program starts with the various rate parameters as input and conducts
# the simulation in three steps.
# (1) Calculating the steady states
# (2) Checking the Turing conditions
# (3) If the Turing conditions are satisfied, solve the PDEs using the implicit-explicit FD-FEM algorithm.
#-----------------------------------------------------------------------------------------------------------------------------
# WE run the scripts twice, one for each parameter set:
# namely (high gamma, high d) and  (low gamma, low d)
#-----------------------------------------------------------------------------------------------------------------------------
#=================================================================================================
#=================================================================================================
#=================================================================================================
#=================================================================================================
# Two data sets, which we print to the user
print('===========================================================================================================================================================\n')
print('\tTwo data sets:\n\n\t\t(1) (high gamma,high d)\n\t\t(2)(low gamma,low d)\n')
print('===========================================================================================================================================================\n')
#=================================================================================================
#=================================================================================================
#=================================================================================================
#=================================================================================================
# The proportion of the surface are over the volume
a = 3
# The cytosolic diffusion is general for all data sets
D = 10000
# maximum concentration in the membrane    
cmax = 3
# Extract all parameters
 #       c2 = c2Vec[i] # Activation rate
c2 = 0.45
#c1 = c1Vec[i] # Influx from the cytosol
c1 = 0.05
#c_1 = c_1Vec[i] # Dissociation from the membrane
c_1 = 0.04
#gamma = gammaVec[i] # Relative strength of the reaction
gamma = 25
#=================================================================================================
# Variable for denoting the dataSet
#dataSet = 0
#=================================================================================================
#=================================================================================================
# Loop over the various data sets and solve the FEM problem
#for dataIter in range(8): #Dont do special cases as for now
for dataIter in range(1): #Dont do special cases as for now
    #dataSet = dataIter + 5
    dataSet = 4

    #--------------------------------------------------------------------------------------------------------------------------------------------------------
    # Now, we need to print this for the user
    print("\t Data set\t%d\tout of\t%d!\n\n\n"  %(dataSet, 4))
    # We now loop through the various values of c_1 to run the calculations
    #for i in range(c2Vec):
    for index in range(1):
        print('===========================================================================================================================================================\n')
        #print("\t Iteration %d out of %d\n\n"  %(i, len(c2Vec)))
        print('===========================================================================================================================================================\n')
        #--------------------------------------------------------------------------------------------------------------------------------------------------------
        # Extract all parameters
 #       c2 = c2Vec[i] # Activation rate
        c2 = 0.45
        #c1 = c1Vec[i] # Influx from the cytosol
        c1 = 0.05
        #c_1 = c_1Vec[i] # Dissociation from the membrane
        c_1 = 0.04
        #gamma = gammaVec[i] # Relative strength of the reaction
        gamma = 25
        i = 4
        #d = 30
        # VECTORS WITH STEADY STATES, HEY?
        dVec =  [5, 10, 15, 30, 50]        
        u0Vec =   [1.2263, 1.2263, 1.2263, 1.2263, 1.2263]
        v0Vec = [0.6276, 0.6276, 0.6276, 0.6276, 0.6276]
        V0Vec = [6.0, 6.0, 6.0, 6.0, 6.0]        
        #--------------------------------------------------------------------------------------------------------------------------------------------------------
        d = dVec[i] # Relative diffusion of active versus the inactive
        # Extract the initial conditions
        u0 = u0Vec[i] # Initial concentraiton active cdc42, u (i.e. cdc42-GTP)
        v0 = v0Vec[i] # Initial concentraiton inactive cdc42, v (i.e. cdc42-GDP)
        V0 = (V0Vec[i] - ( a * (u0+v0)) ) # Initial concentraiton cytosolic cdc42, V (i.e. cdc42-GDI)
        #=================================================================================================
        #=================================================================================================
        #=================================================================================================
        #=================================================================================================
        # Welcome prompt to user
        print('===========================================================================================================================================================\n')
        print('\tFEM-solution to the RD problem modelling Cdc42 activation\n')
        print('===========================================================================================================================================================\n')
        #=================================================================================================
        #=================================================================================================
        #=================================================================================================
        #=================================================================================================
        #-----------------------------------------------------------------------------------------------------------------------------
        # Calculating the FD-FEM solution of the PDE-problem
        #-----------------------------------------------------------------------------------------------------------------------------
        print('===========================================================================================================================================================\n')
        print('\tCalculating the FD-FEM solution of the PDE-problem!\n\n')
        #-----------------------------------------------------------------------------------------
        #import FEMandFD_solver_ImpExp_linear
        #FEMandFD_solver_ImpExp_linear.solvePDE(u0, v0, V0,c1,c_1,c2,cmax,gamma,d,D,a,dataSet,i)

        import FEMandFD_solver_ImpExp_linear_tAdaptive_PoleFinder
        #FEMandFD_solver_ImpExp_linear_tAdaptive.solvePDE(u0, v0, V0,c1,c_1,c2,cmax,gamma,d,D,a,dataSet,i)
        t_pole, ratio_pole, u_max, u_min, poleIndicator = FEMandFD_solver_ImpExp_linear_tAdaptive_PoleFinder.solvePDE(u0, v0, V0, c1,c_1,c2,cmax,gamma,d,D,dataSet,i)

        u_file = open("../u_Stats.txt", "w")
        u_file.write("\t\t d=\t%0.1f\n" % (d)) 
        u_file.write("\t\t time_Pole=\t%0.5f\n" % (t_pole)) 
        u_file.write("\t\t surfaceRatio=\t%0.5f\n" % (ratio_pole))
        u_file.write("\t\t u_max=\t%0.5f\n" % (u_max))
        u_file.write("\t\t u_min=\t%0.5f\n" % (u_min)) 
        u_file.write("\t\t poleIndicator=\t%d\n" % (poleIndicator))       
        u_file.close()
        
        #-----------------------------------------------------------------------------------------
        print('\t\tCalculations are finished!')
