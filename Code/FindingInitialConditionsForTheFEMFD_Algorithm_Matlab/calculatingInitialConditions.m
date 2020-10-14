%% Testing single steady states
% Date : 2019-11-14
% Author: Johannes Borgqvist and Adam Malik
% Description: 
% Given various parameters, we calculate the 
% the steady states and check the Turing
% conditions. We also print this information. 
% The steady states (u0,v0,V0) are later used
% to set the initial conditions in the FEM-FD
% simulations. There is an output called
% ``value'' from the function ``TuringCond''
% which has the value 1.0 for the classic case
% and it has the value 0.5 for the non-classic
% case
clear,clc
%% Define logical variable depending on which intersection of parameter space
alternative = true; % "c_1 vs c2"-plane
%alternative = false; % "c1 vs c2"-plane

%% Defining variables
%------------------------------------------------
% Unclassic
%V0_init = 10.0; % Initial GDI-bound Cdc42
% Classic
V0_init = 6.0; % Initial GDI-bound Cdc42
%------------------------------------------------
a = 3.0; % Quotient between membrane-area and cytosol-volume
cmax = 3.0; % Maximum amount of membrane bound cdc42
%------------------------------------------------
% Unclassic
%d = 1;% Diffusion "Cdc42GDP/Cdc42GTP" parameter: Reactions vs diffusion
% Classic low diffusion
%d = 10;% Diffusion "Cdc42GDP/Cdc42GTP" parameter: Reactions vs diffusion
% Classic high diffusion
d = 5;% Diffusion "Cdc42GDP/Cdc42GTP" parameter: Reactions vs diffusion
%------------------------------------------------
c2 = 0.4; % Activation rate
c_1 = 0.02; % Dissociation rate of Cdc42-GDP from membrane to cytosol
c1 = 0.05; % Import of Cdc42-GDP to membrane from cytosol
%------------------------------------------------
% Alternatives depending on parameter space
if alternative %"c_1 vs c2"
    c_1Max = c_1;
    h = ((c_1Max-c_1)/(10));
    %c_1Vec = c_1:h:c_1Max;
    c_1Vec = c_1;
    nuOfEle = length(c_1Vec);    
else %"c1 vs c2"
    c1Max = c1;
    h = ((c1Max-c_1)/(10));
    %c1Vec = c_1:h:c1Max;    
    c1Vec = c1;
    nuOfEle = length(c1Vec);
end
% Since we use a Newton solver locally around a start guess
% we do multiple start guesses u0 for the active state u
% in the interval $u\in (\sqrt{c_2},\min\left\{c_{\max},V_0/a\right\})$.
nuOfGuesses = 50; % Number of guesses for the steady state
% -------------------------------------------------------------------------
% Loop over the start guesses and display the steady states
for i = 1:nuOfEle
    if alternative
        c_1 = c_1Vec(i);
    else
       c1 = c1Vec(i); 
    end
    %% Calculate the steady states and check if we had Turing or not
    % Calculate the steady states
    [uStar,vStar,VStar] = steadyStateCalculator(c1,c_1,c2,a,V0_init,cmax,nuOfGuesses);
    % Checking the Turing conditions: Classic => value=1.0 and Non-classic=> value=0.5
    [indicator,value,u_SS,v_SS,V_SS] = TuringCond(c1, c_1, c2, d, a, V0_init, cmax,uStar,vStar);
    % Print the information
    if indicator 
        if alternative
            % Print the information
            fprintf('------------------------------------------------------\n');
            fprintf('c_1\tc2\tTuring\tu0\tv0\tV0\n');
            fprintf('------------------------------------------------------\n');
            fprintf('%0.3f\t%0.3f\t%0.1f\t%0.3f\t%0.3f\t%0.5f\n',c_1,c2,value,uStar(1),vStar(1),VStar(1));
            fprintf('------------------------------------------------------\n');
            break;
        else
            % Print the information
            fprintf('------------------------------------------------------\n');
            fprintf('c1\tc2\tTuring\tu0\tv0\tV0\n');
            fprintf('------------------------------------------------------\n');
            fprintf('%0.1f\t%0.3f\t%0.1f\t%0.3f\t%0.3f\t%0.5f\n',c1,c2,value,uStar(1),vStar(1),VStar(1));
            fprintf('------------------------------------------------------\n');            
        end
    else
        % We do not have a steady state for the proposed parameters
        fprintf('Sorry, we do not have a steady state for these parameters...\n');
    end
end
