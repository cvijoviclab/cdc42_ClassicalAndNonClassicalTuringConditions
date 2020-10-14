function v = v_q_down(theta,u)
    % We have the parameters in a vector theta and they are as follows:
    c1 = theta(1); % import from cytosol to membrane
    c_1 = theta(2); % dissociation from membrane to cytosol
    V0_init = theta(4); % Initial concentration of Cdc42-GDI
    a = theta(5); % Quotient between membrane area and cytosol volume
    cmax = theta(6); % Maximum concentration be
    %% Optimisation bit: Define the function for finding the zeros
    A = (cmax - u);
    B = (V0_init - (a * u) );
    K = ( (1/a) * ( (a*A + B) + ( (c_1) / (c1) ) ) );
    L = ( (A * B)/(a));    
    v = 0.5 * (K - sqrt(K.^2 - (4 .* L)));
end