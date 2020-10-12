#==============================================================================================
# STEP 1: SETTING UP THE PROBLEM
#==============================================================================================
#--------------------------------------------------
# Import various packages
#--------------------------------------------------
from dolfin import *
import numpy as np
#--------------------------------------------------
# Python classes for initial conditions
#--------------------------------------------------
class ic(UserExpression):
    def __init__(self, *args, **kwargs):
        self.ic_u = kwargs.pop('ic_u')
        self.ic_v = kwargs.pop('ic_v')
        self.ic_V = kwargs.pop('ic_V')
        self.sigma = kwargs.pop('sigma')
        self.rc = kwargs.pop('rc')
        super(ic,self).__init__(*args, **kwargs)

    def eval(self, values, x):
        if (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] > self.rc ):
            values[0] = self.ic_u + np.random.normal(scale=self.sigma)
            values[1] = self.ic_v + np.random.normal(scale=self.sigma)
            #values[2] = 0
            values[2] = self.ic_V + np.random.normal(scale=self.sigma)
        elif (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] < self.rc ):
            values[0] = 0
            values[1] = 0
            values[2] = self.ic_V + np.random.normal(scale=self.sigma)

    def value_shape(self):
        return(3,)

#--------------------------------------------------
# Define function for solving the PDE
#--------------------------------------------------
def solvePDE(u0, v0, V0, c1,c_1,c2,cmax,gamma,d,D,dataSet,iteration):
    #--------------------------------------------------
    # STRING DEFINING DIRECTORY
    #--------------------------------------------------
    strOriginal = ''
    #--------------------------------------------------
    # DEFINE TIME STEPPING SCHEMES
    #--------------------------------------------------
    T = 5.0 # The final time
    dt = 1e-12 # time step size
    #--------------------------------------------------
    # REDEFINE PARAMETERS FOR THE VARIATIONAL FORMULATION
    #--------------------------------------------------
    # Import to the membrane: Cdc42-GDI->Cdc42-GTP
    c1_orig = c1
    c1 = Constant(c1)
    # Dissociation from the membrane: Cdc42-GTP->Cdc42-GDP
    c_1_orig = c_1
    c_1 = Constant(c_1)
    # Activity of feedback loop: Cdc42-GDP<->Cdc42-GTP
    c2_orig = c2
    c2 = Constant(c2)
    # Relative reaction force in the membrane: gamma
    gammaOrig = gamma
    gamma = Constant(gamma)
    # Relative diffusion of Cdc42-GTP in the membrane: d
    dOrig = d
    d = Constant(d)
    # Diffusion of Cdc42-GDI in the cytosol: D
    D_orig = D
    D = Constant(D)
    # Time step for the Euler forward scheme
    k = Constant(dt)
    # Maximum concentration in the membrane
    cmax_orig = cmax
    cmax= Constant(cmax)
    #--------------------------------------------------
    # READ THE MESH CREATED IN GMSH
    #--------------------------------------------------
    #mesh_cytosol_membrane = Mesh("./mesh_crude.xml")
    #mesh_cytosol_membrane = Mesh("./mesh_MiddleCoarse.xml")
    mesh_cytosol_membrane = Mesh("./mesh_ActuallyCoarse.xml")
    #--------------------------------------------------
    print("\n-----------------------------------------------------------------------------------------------------------------------\n","Setting up the Finite Element Method","\n-----------------------------------------------------------------------------------------------------------------------\n")
    # Let user know that we set up function space
    print("\n\tSetting up functions spaces for test functions\n")
    #==============================================================================================
    # STEP 2: FEM-function spaces, initial conditions and Variation Formulation (VF)
    #==============================================================================================
    # Define function space for system of concentrations
    # IMPORTANT WITH FEM: should be ''tetrahedon'' in three dimensions, ''triangle''
    # in two and ''line'' in one
    P1 = FiniteElement('P', tetrahedron, 1)
    # Define a mixed element since we want three different basis
    element = MixedElement([P1, P1, P1])
    # Collect everything in the function space
    H = FunctionSpace(mesh_cytosol_membrane, element)
    # Define test functions
    phi_1, phi_2, phi_3 = TestFunctions(H)
    # Define functions for concentrations of cdc42
    u, v, V = TrialFunctions(H) # Current time step in finite difference scheme
    #U = TrialFunction(H)
    uPrev = Function(H) # Previous time step in finite difference scheme
    # Let user know that it is done
    print("\n\t\tDONE!\n")
    #--------------------------------------------------
    # INITIAL CONDITIONS
    #--------------------------------------------------
    # Split system functions to access components
    # Let user know that we define initial conditions
    #print("\n\tDefining initial conditions")
    radius_cytosol = 1.0 - 0.1*mesh_cytosol_membrane.hmin()
    #ic = Expression(('pow(x[0],2) + pow(x[1],2) + pow(x[2],2) >pow(rc,2) ? ic_u*(1.0 +  (rand() / 2147483647.0 - 0.5) / 10.0): 0.0','pow(x[0],2) + pow(x[1],2) + pow(x[2],2) > pow(rc,2) ? ic_v*(1.0 + (rand() / 2147483647.0 - 0.5) / 10.0) : 0.0','pow(x[0],2) + pow(x[1],2) +  pow(x[2],2) <pow(rc,2) ? ic_V*(1.0 +  (rand() / 2147483647.0 - 0.5) / 10.0): 0.0'), degree=1, rc = radius_cytosol, ic_u=u0, ic_v=v0, ic_V=V0)
    if dataSet == 1:
        u0Exp = ic(ic_u=u0,ic_v=v0, ic_V=V0,sigma = 0.05,rc=radius_cytosol, element=H.ufl_element())
    else:
        u0Exp = ic(ic_u=u0,ic_v=v0, ic_V=V0,sigma = 0.05,rc=radius_cytosol, element=H.ufl_element())
    #u0Exp = Constant((u0, v0, V0))
    # Prepare initial condition
    uPrev.interpolate(u0Exp)
    #uPrev = interpolate(ic,H)
    # "interpolate" or "project"
    #---------------------------------------
    #--------------------------------------------------
    # SPLIT INTO THE VARIOUS COMPONENTS
    #--------------------------------------------------
    # Split system functions to access components
    u_prev, v_prev, V_prev = split(uPrev)
    # Let user know that it is done
    print("\n\t\tDONE!\n")
    #--------------------------------------------------
    #--------------------------------------------------
    # VARIATIONAL FORMULATION
    #--------------------------------------------------
    # Let user know that we define the variational formulation
    print("\n\tDefining variational formulation")
    # For the variation formulation we have two domains
    # and one interface in which the action occurs:
    # the membrane, the interface and the cytosol. So
    # the variation formulation consists of these three
    # parts.
    #-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    # Define the tangential gradient on the boundary
    n = FacetNormal(mesh_cytosol_membrane)
    def grad_G(w):
        return grad(w) - dot(grad(w), n)*n
    #_________________________________________________
    #PART 0: REACTION TERM AND BOUNDARY REACTION
    #_________________________________________________
    # REACTION FUNCTION
    #f = c2*v - u + u*u_prev*v_prev # maxed implicit for linear
    f = ( (c2 * v) - u + (u_prev * u_prev * v_prev) ) # mixed explicit-implicit
    #f = ( (c2 * v_prev) - u_prev + (u_prev * u_prev * v_prev) ) # fully explicit
    # ROBIN BOUNDARY FUNCTION
    #q = (c1 * V * (cmax - (u+v))) - (c_1 * v) # Alternative Original
    #q = c1*V_prev*(cmax - (u + v)) - c_1*v   # maxed implicit for linear
    #q = c1*V*(cmax - (u_prev + v_prev)) - c_1*v   # alternative implicit for linear
    #q = c1 * (Vinit - 3.0*(u_prev + v_prev)) * (cmax - (u + v)) - (c_1 * v)   # maxed implicit for linear (ODE-version)
    q = (c1 * V_prev * (cmax - (u_prev+v_prev))) - (c_1 * v)   # mixed explicit-implicit
    #q = (c1 * V_prev * (cmax - (u_prev+v_prev))) - (c_1 * v_prev)  # fully explicit
    #_________________________________________________
    #PART 1: MEMBRANE BOUND REACTIONS (dx(3))
    #_________________________________________________
    membraneActive_mass = (u-u_prev)*phi_1*ds
    membraneActive_stiffness = dot(grad_G(u), grad_G(phi_1))*ds
    membraneInactive_mass = (v-v_prev)*phi_2*ds
    membraneInactive_stiffness = d*dot(grad_G(v), grad_G(phi_2))*ds
    membrane_reaction = gamma*f*(phi_2 - phi_1)*ds
    #membrane = membraneActive_mass + membraneActive_stiffness
    #membrane += membraneInactive_mass + membraneInactive_stiffness
    #membrane += membrane_reaction
    #_________________________________________________
    #PART 2: REACTIONS OVER THE INTERFACE (dS(2))
    #_________________________________________________
    interface_reaction = gamma*q*(phi_3 - phi_2)*ds
    #interface = interface_reaction
    #_________________________________________________
    #PART 3: DIFFUSION IN THE CYTOSOL (dx(2))
    #_________________________________________________
    cytosol_mass = (V - V_prev)*phi_3*dx
    cytosol_stiffness = D*dot(grad(V), grad(phi_3))*dx
    # cytosol = cytosol_mass + cytosol_stiffness
    #_________________________________________________
    #PART 4: THE VARIATION FORMULATION
    #_________________________________________________

    Mass_form = membraneActive_mass + membraneInactive_mass + cytosol_mass
    Stiffness_form = membraneActive_stiffness + membraneInactive_stiffness + cytosol_stiffness
    #Stiffness_form = cytosol_stiffness
    Reaction_form = membrane_reaction + interface_reaction
    #F = membrane + interface + cytosol
    #--------------------------------------------------
    # FILES IN WHICH WE STORE DATA, HEY?
    #--------------------------------------------------
    # Let user know that we save the initial conditions
    # and that we set up the time stepping
    print("\n\tSave initial conditions and set up FD scheme")
    print("\t\tWe will save solutions in pvd-files in a folder called 'reaction_system'","\n\t\tIt contains three subfolders:\n","\t\t\t(1.) Active\n\t\t\t(2.) Inactive\n\t\t\t(3.) Cytosolic\n")
    #======================================================================================
    # Create VTK files for visualization output
    #======================================================================================    
    strOriginal = '../Results'
    # Create the folder structure
    if dataSet <5:
        strOriginal = strOriginal + '/Classic'
        # We have four subfolders, hey?
        if dataSet == 1:
            strOriginal = strOriginal + '/c_1VSc2'
        elif dataSet == 2:
            strOriginal = strOriginal + '/c1VSc2'            
        elif dataSet == 3:
            strOriginal = strOriginal + '/increasing_gamma'            
        else:
            strOriginal = strOriginal + '/increasing_d'
        # Now, we add the parameter number
        strOriginal = strOriginal + '/paramSet_' + str(int(iteration))
    elif dataSet>4 and dataSet <9:
        strOriginal = strOriginal + '/Unclassic'
        # We have four subfolders, hey?
        if dataSet == 5:
            strOriginal = strOriginal + '/c_1VSc2'
        elif dataSet == 6:
            strOriginal = strOriginal + '/c1VSc2'            
        elif dataSet == 7:
            strOriginal = strOriginal + '/increasing_gamma'            
        else:
            strOriginal = strOriginal + '/increasing_d'
        # Now, we add the parameter number
        strOriginal = strOriginal + '/paramSet_' + str(int(iteration))        
    else:
        strOriginal = strOriginal + '/SpecialCases'
        # We have four subfolders, hey?
        if iteration == 1:
            strOriginal = strOriginal + '/d_EQUALS_ONE'
        elif iteration == 2:
            strOriginal = strOriginal + '/noPattern_Classic'            
        else:
            strOriginal = strOriginal + '/noPattern_Unclassic'                                
        # We finish the job and define the final name of the directory
    strOriginal = strOriginal + '/ActiveCdc42_ConcentrationProfileOverTime/cdc42GTP.pvd'
    print('%s' %(strOriginal))
    # Save the final file
    vtkfile_u = File(strOriginal)
    #======================================================================================
    #======================================================================================        
    # Create progress bar
    #progress = Progress('Time-stepping',num_steps)
    #set_log_level(LogLevel.PROGRESS)

    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    # Save solution to file (VTK)
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    U = Function(H)
    u, v, V = split(U)

    U.assign(uPrev)
    t = 0.0
    _u, _v, _V = U.split()
    vtkfile_u << (_u, t)
    # Let user know that it is done
    print("\n\t\tDONE!\n")
    print("\n-----------------------------------------------------------------------------------------------------------------------\n","The Finite Element Method starts","\n-----------------------------------------------------------------------------------------------------------------------\n")
    #==============================================================================================
    # STEP 3: Finite differenses (FD) in time FEM in space(=> VF)
    #==============================================================================================

    ### Mass for stability check ###
    mass_u = u*ds
    mass_v = v*ds
    mass_V = V*dx

    membrane_area_int = Constant(1.0)*ds(domain=mesh_cytosol_membrane)
    cell_volume_int = Constant(1.0)*dx(domain=mesh_cytosol_membrane)
    membrane_area = assemble(membrane_area_int)
    cell_volume = assemble(cell_volume_int)

    ### Define residual for temporal adaptivity ###
    # Reaction factors for residual
    f_res = c2*v - u + u*u*v
    q_res = c1*V*(cmax - (u+v)) - c_1*v

    # Membrane terms for residual
    membraneActive_res = (u-u_prev)*phi_1*ds #+ dot(grad_G(u), grad_G(phi_1))*ds
    membraneInactive_res = (v-v_prev)*phi_2*ds # + d*dot(grad_G(v), grad_G(phi_2))*ds
    membraneReaction_res = gamma*f_res*(phi_2 - phi_1)*ds
    membrane_res = membraneActive_res + membraneInactive_res + membraneReaction_res

    # Interface terms for residual
    interface_res = gamma*q_res*(phi_3 - phi_2)*ds

    # Cytosol terms for residual
    cytosol_res = (V - V_prev)*phi_3*dx + D*dot(grad(V), grad(phi_3))*dx

    # Residual
    Residual_form = membrane_res + interface_res + cytosol_res

    ### Assemble parts of the linear system ###
    # Matrices on LHS
    Mass_matrix = assemble(lhs(Mass_form), keep_diagonal = True)
    Mass_matrix.ident_zeros()
    StiffnessReaction_matrix = assemble(lhs(Stiffness_form + Reaction_form), keep_diagonal = True)
    StiffnessReaction_matrix.ident_zeros()

    # Forms for vectors on RHS
    Mass_form_rhs = rhs(Mass_form)
    StiffnessReaction_form_rhs = rhs(Stiffness_form + Reaction_form)
    # Start time-stepping
    t_it = 0
    #TOL_tadapt = 1e-1*T
    TOL_tadapt = 1e-2*T
    dt_max = 0.1*T
    while t < T:
        # Update current time and iteration number
        t_it += 1
        t += dt
        k = Constant(dt)

        print("Step\t",t_it,"\t, dt = \t",dt,"\t, time\t",t,"\tof\t",T)

        # Assemble system matrix and vector
        A = Mass_matrix + dt*StiffnessReaction_matrix
        b = assemble(Mass_form_rhs + k*StiffnessReaction_form_rhs)

        # Solve linear variational problem for time step
        solve(A,  U.vector(), b)

        # Check mass conservation
        u_tot = assemble(mass_u)
        v_tot = assemble(mass_v)
        V_tot = assemble(mass_V)
        u_ode = u_tot/membrane_area
        v_ode = v_tot/membrane_area
        V_ode = V_tot/cell_volume
        mass_tot = u_tot + v_tot + V_tot
        print("Total mass = ", mass_tot)
        #print("u_tot = ", u_tot, ", v_tot = ", v_tot ,", V_tot = ", V_tot)
        print("u_avg = ", u_ode, ", v_avg = ", v_ode ,", V_avg = ", V_ode)
        #print("V0 - 3*(u + v) = ", Vinit - 3.0*(u_ode + v_ode))

        # Check parameter values
        #print('c1 = ', c1.values()[0], 'c_1 = ', c_1.values()[0],
        #      'c2 = ', c2.values()[0], 'Vinit = ', Vinit,
        #       'cmax = ', cmax.values()[0], 'gamma = ', gamma.values()[0])

        # Check q and f
        #q_val_onlyuv = c1_orig*(Vinit - 3.0*(u_ode + v_ode))*(cmax_orig - (u_ode+v_ode)) - c_1_orig*v_ode
        q_val = c1_orig*V_ode*(cmax_orig - (u_ode+v_ode)) - c_1_orig*v_ode
        f_val = c2_orig*v_ode - u_ode + u_ode*u_ode*v_ode
        print("q = ", q_val, ", f = ", f_val, "\n")
        if u_tot < 0.0 or  v_tot < 0.0 or V_tot < 0.0:
            print("ERROR: Negative amounts detected ==> System unstable")
            print("TERMINATING SIMULATION!!!")
            break

        # Compute next timestep adaptively with residual
        dt_old = dt
        R = assemble(Residual_form)
        dt = TOL_tadapt/norm(R, 'l2')
        dt = min(2.0*dt_old*dt/(dt_old + dt), dt_max)

        # Save solution to file (VTK)
        #-------------------------------------------------------------------------
        #-------------------------------------------------------------------------
        if dataSet == 1:
            if  (t_it % 15 ==0) or t + dt >= T:
                _u, _v, _V = U.split()
                vtkfile_u << (_u, t)
        else:
            if  (t_it % 15 ==0) or t + dt >= T:
                _u, _v, _V = U.split()
                vtkfile_u << (_u, t)

        # Update previous solution
        uPrev.assign(U)


