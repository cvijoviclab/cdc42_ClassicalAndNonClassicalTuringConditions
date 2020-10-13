#==============================================================================================
# FEM-SIMULATION OF  Cdc42-POLARISATION
#==============================================================================================
# DATE:
# 2020-10-08
# WRITTEN BY:
# Johannes Borgqvist, Adam Malik and Carl Lundholm
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
import pandas as pd 
#import pymp
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
        strBase = "../../../Results/increasingGamma/Classical/21"        
        os.mkdir(strBase)
        
    else: # Non-classic, increasing d
        c1 =  0.05
        c_1 = 0.03
        c2 = 0.15
        u0 = 0.9059
        v0 = 0.9332
        strBase = "../../../Results/increasingGamma/NonClassical/21"
        os.mkdir(strBase)         
        
     # Initial concentraiton cytosolic cdc42, V (i.e. cdc42-GDI)
    V0 = (V0_init - ( a * (u0+v0)) )
    #--------------------------------------------------------------------------------------------------------------------------------------------------------
    # number of gamma values corresponding to the relative strength of the reaction
    # versus diffusion
    gammaVec = np.arange(10,163,10)
    # Save the gamma vector
    # convert array into dataframe 
    DF0 = pd.DataFrame(gammaVec) 
    # save the dataframe as a csv file 
    DF0.to_csv(strBase + "gammaVec.csv")    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------------------------------------------------- 
    #--------------------------------------------------------------------------------------------------------------------------------------------------------
    # Repeat for stochasticity
    nuOfRepeats = 1
    print('\tWe repeat %d times!\n\n' % (nuOfRepeats))    
    # Define the matrices in which we will save four parameters:
    # uMax, uMin, tPole and ratioPole
    u_maxAvg = np.zeros((len(gammaVec),nuOfRepeats),dtype=np.double)
    u_minAvg = np.zeros((len(gammaVec),nuOfRepeats),dtype=np.double)
    t_poleAvg = np.zeros((len(gammaVec),nuOfRepeats),dtype=np.double)
    ratio_poleAvg = np.zeros((len(gammaVec),nuOfRepeats),dtype=np.double)
    #-------------------------------------------------------------------------------------------------------------------------------------------------------     #-------------------------------------------------------------------------------------------------------------------------------------------------------      
    # We loop over the number of repetition (for stochasticity in the ICs)
    for innerIndex in range(nuOfRepeats):                     
        # We now loop through the various values of c_1 to run the calculations
        for i in range(len(gammaVec)):
            #--------------------------------------------------------------------------------------------------------------------------------------------------------
            #--------------------------------------------------------------------------------------------------------------------------------------------------------
            gamma = gammaVec[i] # Extract the reaction strength
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
            # Save the values to the matrix, hey?
            u_maxAvg[i,innerIndex] = u_max
            u_minAvg[i,innerIndex] = u_min
            t_poleAvg[i,innerIndex] = t_pole
            ratio_poleAvg[i,innerIndex] = ratio_pole
            #-----------------------------------------------------------------------------------------
            print('\t\tCalculations are finished!')

            
    # Finally, we save the matrices in csv files which will
    # be used for post-processing
    #---------------------------------------------------------                
    # u_min
    # convert array into dataframe 
    DF1 = pd.DataFrame(u_minAvg) 
    # save the dataframe as a csv file 
    DF1.to_csv(strBase + "uMin.csv")
    #---------------------------------------------------------                
    # u_max
    # convert array into dataframe 
    DF2 = pd.DataFrame(u_maxAvg) 
    # save the dataframe as a csv file 
    DF2.to_csv(strBase + "uMax.csv")
    #---------------------------------------------------------                
    # t_pole
    # convert array into dataframe 
    DF3 = pd.DataFrame(t_poleAvg) 
    # save the dataframe as a csv file 
    DF3.to_csv(strBase + "tPole.csv")
    #---------------------------------------------------------                
    # ratio_pole
    # convert array into dataframe 
    DF4 = pd.DataFrame(ratio_poleAvg) 
    # save the dataframe as a csv file 
    DF4.to_csv(strBase + "ratioPole.csv")    
    #---------------------------------------------------------                
