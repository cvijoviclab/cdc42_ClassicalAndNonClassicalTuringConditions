#==============================================================================================
# FEM-SIMULATION OF  Cdc42-POLARISATION
#==============================================================================================
# DATE:
# 2020-03-09
# WRITTEN BY:
# Johannes Borgqvist, Adam Malik and mainly Carl Lundholm
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
# Nice package to create vectors, innit?
import numpy as np
import os
# The proportion of the surface are over the volume
a = 3
# The cytosolic diffusion is general for all data sets
D = 10000
# maximum concentration in the membrane    
cmax = 3
# Initial concentration V0
V0_init = 6.0
# The reaction strength parameter
d = 10
#=================================================================================================
# Variable for denoting the dataSet
#=================================================================================================
#=================================================================================================
# Loop over the various data sets and solve the FEM problem
for dataIter in range(2): #Dont do special cases as for now
    dataSet = dataIter + 1
    #dataSet = 2
    print("\t Dataset %d\n" %(dataSet))
    #-----------------------------------------------------------------------------------------------------------------------------
    # IMPORT PARAMETERS FROM THE PARAMETER FILE
    #-----------------------------------------------------------------------------------------------------------------------------
    # We have four parameter sets, hey?
    # The first two are the classic Turing and
    # the second two are the unclassic Turing.
    # The first and third are with high values of gamma, while
    # the second and fourth are with low values of gamma.
    #--------------------------------------------------------------------------------------------------------------------------------------------------------
    if dataSet == 1:# Classic, increasing d
        c1 =  0.05
        c_1 =  0.04
        c2 = 0.45
        u0 = 1.2263
        v0 = 0.6276
        strBase = '../Classic_tAdaptive'
        os.mkdir("../Classic_tAdaptive") 
    else: # Non-classic, increasing d
        c1 =  0.05
        c_1 = 0.03
        c2 = 0.15
        u0 = 0.9059
        v0 = 0.9332
        strBase = '../NonClassic_tAdaptive'
        os.mkdir("../NonClassic_tAdaptive")         
        
     # Initial concentraiton cytosolic cdc42, V (i.e. cdc42-GDI)
    V0 = (V0_init - ( a * (u0+v0)) )
    #--------------------------------------------------------------------------------------------------------------------------------------------------------

    # number of diffusion coefficients
    gammaVec = np.arange(10,163,10)
    #gammaVec = np.arange(10,23,10)
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------------------------------------------- 
    #--------------------------------------------------------------------------------------------------------------------------------------------------------    # Repat for stochasticity
    print('\tWe repeat three times!\n\n')
    nuOfRepeats = 3    
    # Files for saving the tex-files, hey? THE AVERAGE VALUES
    uMaxAvg_file = open("%s/uMaxAvg.tex" % (strBase), "w")
    uMinAvg_file = open("%s/uMinAvg.tex" % (strBase), "w")
    ratioPoleAvg_file = open("%s/ratioPoleAvg.tex" % (strBase), "w")
    tPoleAvg_file = open("%s/tPoleAvg.tex" % (strBase), "w")        
    # Write opening line for the plots, hey? 
    uMaxAvg_file.write("\\addplot[color=blue] coordinates {")
    uMinAvg_file.write("\\addplot[color=blue] coordinates {")
    ratioPoleAvg_file.write("\\addplot[color=blue] coordinates {")
    tPoleAvg_file.write("\\addplot[color=blue] coordinates {")
    # Define the average parameters
    u_maxAvg = np.zeros(len(gammaVec))
    u_minAvg = np.zeros(len(gammaVec))
    t_poleAvg = np.zeros(len(gammaVec))
    ratio_poleAvg = np.zeros(len(gammaVec))
    #-------------------------------------------------------------------------------------------------------------------------------------------------------     #-------------------------------------------------------------------------------------------------------------------------------------------------------      
    for innerIndex in range(nuOfRepeats):             
        #gammaVec = 200
        # Files for saving the tex-files, hey? THE INDIVIDUAL FILES
        uMax_file = open("%s/uMax_%d.tex" % (strBase, innerIndex), "w")
        uMin_file = open("%s/uMin_%d.tex" % (strBase, innerIndex), "w")
        ratioPole_file = open("%s/ratioPole_%d.tex" % (strBase, innerIndex), "w")
        tPole_file = open("%s/tPole_%d.tex" % (strBase, innerIndex), "w")        
        # Write opening line for the plots, hey? 
        uMax_file.write("\\addplot[color=blue] coordinates {")
        uMin_file.write("\\addplot[color=blue] coordinates {")
        ratioPole_file.write("\\addplot[color=blue] coordinates {")
        tPole_file.write("\\addplot[color=blue] coordinates {")            
        # We now loop through the various values of c_1 to run the calculations
        for i in range(len(gammaVec)):
            #--------------------------------------------------------------------------------------------------------------------------------------------------------
            #--------------------------------------------------------------------------------------------------------------------------------------------------------
            #gamma = 200
            gamma = gammaVec[i] # Relative diffusion of active versus the inactive
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
            import FEMandFD_solver_ImpExp_linear_tAdaptive_PoleFinder
            t_pole, ratio_pole, u_max, u_min, poleIndicator, strBase_2 = FEMandFD_solver_ImpExp_linear_tAdaptive_PoleFinder.solvePDE(u0, v0, V0, c1,c_1,c2,cmax,gamma,d,D,dataSet,i,innerIndex)
        
            uMax_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, u_max))
            uMin_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, u_min))
            tPole_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, t_pole))
            ratioPole_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, ratio_pole))
            u_maxAvg[i] += u_max
            u_minAvg[i] += u_min
            t_poleAvg[i] += t_pole
            ratio_poleAvg[i] += ratio_pole
            #-----------------------------------------------------------------------------------------
        print('\t\tCalculations are finished!')
        # Close the two tex-files which we write to as well
        uMax_file.write("};")
        uMin_file.write("};")
        tPole_file.write("};")
        ratioPole_file.write("};")                        
        uMax_file.close()
        uMin_file.close()
        tPole_file.close()
        ratioPole_file.close()                  
   
    for i in range(len(gammaVec)):
        # STOP REPATING=>CALCULATE AVERAGES!
        u_maxAvg[i] = u_maxAvg[i]/nuOfRepeats
        u_minAvg[i] = u_minAvg[i]/nuOfRepeats
        t_poleAvg[i] = t_poleAvg[i]/nuOfRepeats
        ratio_poleAvg[i] = ratio_poleAvg[i]/nuOfRepeats
        # SAVE TO FILES
        uMaxAvg_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, u_maxAvg[i]))
        uMinAvg_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, u_minAvg[i]))
        tPoleAvg_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, t_poleAvg[i]))
        ratioPoleAvg_file.write("\t\t (%d\t ,\t%0.5f\t ) \n" % (gamma, ratio_poleAvg[i]))        

    # Close the two tex-files which we write to as well
    uMaxAvg_file.write("};")
    uMinAvg_file.write("};")
    tPoleAvg_file.write("};")
    ratioPoleAvg_file.write("};")                        
    uMaxAvg_file.close()
    uMinAvg_file.close()
    tPoleAvg_file.close()
    ratioPoleAvg_file.close()                  
