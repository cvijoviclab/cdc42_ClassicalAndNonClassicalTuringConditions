%% Parameter plot Turing space
% Date : 2020-01-20
% Author: Johannes Borgqvist and Adam Malik
% Description:
% We find a local steady state numerically, and check the condition
% for Turing instability. We save the parameters that satisfy criteria.
% We use a local search and least square minimisation. If we are within
% a tolerance value we say that where we end up is a steady state. In this
% case we plug in the steady state into the conditions and see whether or
% not if the parameters at hand actually satisfies the Turing instability
% criteria.
clear,clc

%% Print welcome prompt
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
fprintf('\n\n\t WELCOME TO THE "Turing parameter space"-SCRIPT!\n\n')
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
pause(2)
fprintf('\tDATE WHEN SCRIPT WAS WRITTEN:\n');
pause(1)
fprintf('\t\t2020-01-20\n');
pause(2)
fprintf('\tSCRIPT WAS WRITTEN BY:\n');
pause(2)
fprintf('\t\tJohannes Borgqvist\tand\tAdam Malik\n');
pause(2)
fprintf('\tDESCRIPTION:\n');
pause(1)
fprintf('\t\t The program calculates the steady states for given parameters\n');
fprintf('\t\t and evaluates if the system undergoes diffusion driven instability!\n\n');
pause(2)
fprintf('-------------------------------------------------------------');
fprintf('-------------------------------------------------------------\n');
pause(4)


%% DECIDE HOW MANY ITERATIONS WE SHOULD HAVE, i.e. the size of the image 
nuOfIter = 1500;
%nuOfIter = 200;
nuOfGuesses = 50; %For the numeric solver
dataSetVec = [1, 2];
% -------------------------------------------------------------------------
% We also have to define if we want to look at classic or unclassic
% Turing instability. So, we define a variable for this called
% "desiredSpace" which takes the value 1.0 for classic Turing instability
% and the value 0.5 for unclassic Turing instability
%desiredSpace = 1.0;
classic = 1.0;
nonClassic = 0.5;
% -------------------------------------------------------------------------
%% LOOP OVER DATA SET
% Loop over five data sets
for iter = 1:length(dataSetVec)
    %% Define Parameters
    dataSet = dataSetVec(iter);
    fprintf('\n\t\tDATA SET %d!\n\n',dataSet);
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    %% Define Parameters
    fprintf('\n\t\tALLOCATE PARAMETERS:\t');    
    if dataSet == 1
        %% Defining variables
        V0_init = 6.0; % Initial GDI-bound Cdc42
        a = 3.0; % Quotient between membrane-area and cytosol-volume
        cmax = 3.0; % Maximum amount of membrane bound cdc42
        c1_Fixed = 0.05; % Import f Cdc42-GDP to membrane from cytosol
        %c_1 = 0.05; % Import f Cdc42-GDP to membrane from cytosol
        %d = 5;% Diffusion "Cdc42GDP/Cdc42GTP" parameter: Reactions vs diffusion
        % -------------------------------------------------------------------------
        %% Define Parameter Vectors to loop over
        % Interval for c_1
        c_1_min = 0.01;
        c_1_max = 0.2;
        % Interval for c2
        c2_min = 0.01;
        c2_max = 2.00;
        % Create the vectors bsed on these limits
        c_1_Vec = linspace(c_1_min,c_1_max,nuOfIter);
        c2_Vec = linspace(c2_min,c2_max,nuOfIter);
        % Create a parameter mesh
        [X,Y] = meshgrid(c_1_Vec,c2_Vec);
    else
        %% Defining variables
        V0_init = 6.0; % Initial GDI-bound Cdc42
        a = 3.0; % Quotient between membrane-area and cytosol-volume
        cmax = 3.0; % Maximum amount of membrane bound cdc42
        %c1 = 0.05; % Import f Cdc42-GDP to membrane from cytosol
        c_1_Fixed = 0.05; % Import f Cdc42-GDP to membrane from cytosol
        %d = 5;% Diffusion "Cdc42GDP/Cdc42GTP" parameter: Reactions vs diffusion
        % -------------------------------------------------------------------------
        %% Define Parameter Vectors to loop over
        % Interval for c_1
        % Interval for c1
        c1_min = 0.01;
        c1_max = 10.00;
        % Interval for c2
        c2_min = 0.01;
        c2_max = 2.00;
        % Create vectors based on these limits
        c1_Vec = linspace(c1_min,c1_max,nuOfIter);
        c2_Vec = linspace(c2_min,c2_max,nuOfIter);
        % Create a mesh using the vectors
        [X,Y] = meshgrid(c1_Vec,c2_Vec);
    end    
    %%
    % Define the third thing for plotting
    TuringSpace = zeros(nuOfIter,nuOfIter);    
    % Prompt to the user
    pause(1)
    fprintf('\t\t\tDONE!\n\n');
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    pause(2)
    %%   LOOPING THROUGHT THE PARAMETERS SPACE
    tStart = 0;
    tEnd = 0;
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('\n\t\tLOOPING THROUGH THE PARAMETER SPACE\n\n');
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('\n\t\t\tWe will conduct %d iterations\n',nuOfIter);
    fprintf('\t\t\tThis will take some time...\n\n');
    pause(1.5)
    fprintf('\t\t\tGo ahead and have a cup of coffee!\n\n\n');
    pause(1)
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    pause(2)
    % Loop through the parameter space
    parfor i = 1:nuOfIter
        tStart = tic;
        fprintf('\t\t\tIteration %d out of %d\n',i,nuOfIter);
        for j = 1:nuOfIter
            %-------------------------------------------------------
            % Set up the parameters and the critical values
            %-------------------------------------------------------
            % The activation rate
            c2 = Y(i,j);
            if mod(dataSet,2)~=0                 
                % The dissociation rate
                c_1 = X(i,j);
                % The cytosolic flux
                c1 = c1_Fixed;
            else
                % The dissociation rate
                c_1 = c_1_Fixed;
                % The cytosolic flux
                c1 = X(i,j);
            end            
            %---------------------------------------------------------------------------------
            % Calculate the steady states
            [uStar,vStar,VStar] = steadyStateCalculator(c1,c_1,c2,a,V0_init,cmax,nuOfGuesses);
            %[uStar,vStar, VStar] = SymbolicSS(c1,c_1,c2,a,V0_init,cmax);
            %---------------------------------------------------------------------------------
            % Checking the Turing conditions: d=5
            [indicator1,value1,u_SS,v_SS,V_SS] = TuringCond(c1, c_1, c2, 5, a, V0_init, cmax,uStar,vStar);
            % Checking the Turing conditions: d=10
            [indicator2,value2,u_SS,v_SS,V_SS] = TuringCond(c1, c_1, c2, 10, a, V0_init, cmax,uStar,vStar);
            % Checking the Turing conditions: d=30
            [indicator3,value3,u_SS,v_SS,V_SS] = TuringCond(c1, c_1, c2, 30, a, V0_init, cmax,uStar,vStar);            
            %---------------------------------------------------------------------------------
            % Save parameters if we have the desired instability for the
            % given parameters
            if indicator1 && (value1 == nonClassic)
                TuringSpace(i,j) = value1;                            
            else
                if indicator1 && (value1 == classic) 
                    TuringSpace(i,j) = 1.0;
                else
                    if indicator2 && (value2 == classic) && (TuringSpace(i,j)==0)
                        TuringSpace(i,j) = 1.5;
                    else
                        if indicator3 && (value3 == classic) && (TuringSpace(i,j)==0)
                            TuringSpace(i,j) = 2.0;
                        end
                    end
                end
            end
        end
        % Stop taking time
        tEnd = toc(tStart);
        fprintf('\t\t\t\tElapsed Time:\t%8.2f\tseconds\n', round(tEnd,2));
    end
    % Refine the Turing space to remove NANs
    M = isnan(TuringSpace);
    [r,c] = find(M==1);
    TuringSpace(r,c) = 0;    
    % Prompt to the user
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------\n');
    fprintf('\tCalculations are done!\n\tLet us plot the space!\n')    
    pause(2)        
    %% Plot parameter space and save the plot
    fig = figure('units','normalized','outerposition',[0 0 1 1],...
        'PaperPositionMode','auto');
    clf
    contourf(X,Y,TuringSpace)
    %----------------------------------------------------------------------------
    Viridis = viridis(5);
    colormap(Viridis);
    set(gca,'Ticklabelinterpreter','latex','Fontsize',40);
    %----------------------------------------------------------------------------
    hold off
    caxis([0,2.0]);
    if (dataSet == 1)
        print(fig, '-depsc', '-r0', './Figures/paramPlot/c_1VSc2.eps');
    else
        print(fig, '-depsc', '-r0', './Figures/paramPlot/c1VSc2.eps');
    end
    close(fig);    
    TuringSpace = zeros(nuOfIter,nuOfIter);
end
