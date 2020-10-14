function [indicator,value,u_SS,v_SS,V_SS] = TuringCond(c1, c_1, c2, d, a, V0_init, cmax,uStar,vStar)
% Set an indicator to check whether we have Turing
% or not
indicator = false;
value = 0;
% Define new poarameter m used for the integration
m = V0_init/a;
% Define the partial derivatives
f_u =@(u,v) ( (2*u*v) - 1);
f_v =@(u,v)  ( c2 + (u*u) );
%q_u =@(u,v)(-c1*a) * (cmax + m - ( 2* (u+v) ) ); % the old way
q_u =@(u,v)(-c1*a) * (m - (u+v) ); % Total derivative
q_v =@(u,v) q_u(u,v) - c_1;
q_V =@(u,v)  c1 * (cmax - (u+v) );
VPrime=-a;
% Allocate memory for the successful steady state
u_SS = 0;
v_SS = 0;
V_SS = 0;
% Check that we even have a steady state
if ~(isempty(uStar))
    % Loop through the steady states and
    % see if we have Turing conditions or not
    for i = 1:length(uStar)
        % Extract the steady states used  to
        % check the conditions
        u1 =  uStar(i);
        v1 =  vStar(i);
        % Now we calculate a bunch of conditions
        % Negative trace of homogeneous system
        cond1 = ( ( f_u(u1,v1) -  f_v(u1,v1) +  q_v(u1,v1) + (q_V(u1,v1) * VPrime)) < 0);
        % Positive determinant of homogeneous system
        Term1 = f_u(u1,v1) * ( q_v(u1,v1) + (q_V(u1,v1) * VPrime) );
        Term2 = f_v(u1,v1) * ( q_u(u1,v1) + (q_V(u1,v1) * VPrime) );
        cond2 = ( (Term1 - Term2 ) > 0);
        % Classic Turing:
        cond3 = ( ( ( f_u(u1,v1) * q_v(u1,v1) ) - ( f_v(u1,v1) * q_u(u1,v1) ) ) >= 0);
        cond4 = ( ( (d*f_u(u1,v1)) -f_v(u1,v1)+q_v(u1,v1) ) > 0);
        Term3 = ( (d*f_u(u1,v1)) -f_v(u1,v1)+q_v(u1,v1) );
        Term4 = ( ( f_u(u1,v1) * q_v(u1,v1) ) - ( f_v(u1,v1) * q_u(u1,v1) ) );
        Q = ( (Term3 * Term3) - (4*d*Term4));
        cond5 = ( ( (Term3 * Term3) - (4*d*Term4)) >= 0 );
        % Unclassic Turing:
        cond7 = ( Term4 < 0 );
        cond8 =  ( ( (1/(2*d)) * (Term3 + sqrt(Q)) ) > 0 );
        if cond1 && cond2 && cond3 && cond4 && cond5
            % We had a classic Turing steady state
            % Change the indicator
            indicator = true;
            % Set the value to 1
            value = 1.0;
            % We save the steady state values, hey?
            u_SS = u1;
            v_SS = v1;
            V_SS = V0_init-(a*(u1+v1));
            % We break the loop because we are satisfied
            break
        elseif cond1 && cond2 && cond7 && cond8
            % We had a unclassic Turing steady state
            % Change the indicator
            indicator = true;
            % Set the value to 1
            value = 0.5;
            % We save the steady state values, hey?
            u_SS = u1;
            v_SS = v1;
            V_SS = V0_init-(a*(u1+v1));            
            % We break the loop because we are satisfied
            break
        end
    end
end
end