function [uStar,vStar,VStar] = steadyStateCalculator(c1,c_1,c2,a,V0_init,cmax,nuOfGuesses)
% Define the temporary parameter vector for the optimisation
theta = [c1, c_1, c2, V0_init, a, cmax];
% Define the critical point
uCrit = sqrt(c2);
%uCrit = 0;
% Define the value of m
m = (V0_init/a);
% Define the maximum value
uMax = min(cmax,m);
% Define the step in the u direction
uStep = ( (uMax - uCrit) / (nuOfGuesses) );
% The vector f steady states
uGuesses = (uCrit:uStep:uMax)';
% Outputs
uStar = zeros(length(uGuesses),1);
vStar = zeros(length(uGuesses),1);
VStar = zeros(length(uGuesses),1);
% We do not want to display everything for fzero
options = optimset('Display','off');
% Loop through the startguess
for i = 1:length(uGuesses)
    %---------------------------------
    % Steady state calculation
    %---------------------------------
    % We start by solving the upward curve and
    % then the downward curve
    % Define a start guess
    u0 = uGuesses(i,1);
    % Solve the equation
    uTemp = fzero(@(u)(v_q_down(theta,u)-v_f(c2,u)),u0,options);
    % Find the v-value
    vTemp = v_f(c2,uTemp);
    % Check the big V
    VTemp = V0_init - (a * (uTemp+vTemp) );
    % Check if this is a steady state?
    if ( abs( q(uTemp,vTemp,V0_init,a,cmax,c1,c_1)) < 0.001) && ( abs( f(uTemp,vTemp,c2)) < 0.001) && VTemp>0 && ((uTemp+vTemp)<min(cmax,m))
        % Add the u-value
        uStar(i,1) = uTemp;
        % Add the v-value
        vStar(i,1) = vTemp;
        % Add the V-value
        VStar(i,1) = VTemp;
    else
        % Solve the equation
        uTemp = fzero(@(u)(v_q_up(theta,u)-v_f(c2,u)),u0,options);
        % Find the v-value
        vTemp = v_f(c2,uTemp);
        % Check the big V
        VTemp = V0_init - (a * (uTemp+vTemp) );        
        % Check if this is a steady state?
        if ( abs( q(uTemp,vTemp,V0_init,a,cmax,c1,c_1)) < 0.001) && ( abs( f(uTemp,vTemp,c2)) < 0.001) && VTemp>0 && ((uTemp+vTemp)<min(cmax,m))
            % Add the u-value
            uStar(i,1) = uTemp;
            % Add the v-value
            vStar(i,1) = vTemp;
            % Add the V-value
            VStar(i,1) = VTemp;        
        end
    end
end
% Remove the zero-elements, hey?
[r,c] = find(uStar == 0);
uStar(r) = [];
vStar(r) = [];
VStar(r) = [];
end